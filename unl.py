import pandas as pd
import numpy as np
from math import cos, sin, exp
import matplotlib.pyplot as plt

class Controller(object):
    def __init__(self, zd=np.array([[.0],[2.2]])):
        self.a0 = 2.5
        self.a1 = .1
        self.k1 = .4
        self.k2 = .4
        self.zd = zd
        self.J  = np.array([[0.,-1.],[1.,0.]])
        self.internalTime = 0.
	
    def derivs(self,zd, term2):
        return -self.a1*self.zd + term2 * np.dot(self.J,self.zd)
    
    def calculate(self,errorGlobal, currentPose, desiredVelocity, dt):
        x,y,theta = 0,1,2
        #print errorGlobal[x], errorGlobal[y], errorGlobal[theta]
        
        # -------------------------------------------------
        # Transform errorGlobal to [w,z]
        w = (-errorGlobal[theta]*cos(currentPose[theta]) + 2.*sin(currentPose[theta]))*errorGlobal[x]
        + (-errorGlobal[theta]*sin(currentPose[theta]) - 2.*cos(currentPose[theta]))*errorGlobal[y]
        z = np.array([
            [errorGlobal[theta]], 
            [cos(currentPose[theta])*errorGlobal[x] + sin(currentPose[theta])*errorGlobal[y]]
            ])

        # -------------------------------------------------
        # local definitions
        f      = 2.*(desiredVelocity[1]*z[1] - desiredVelocity[0]*np.sin(z[0]))
        delta  = self.a0 * np.exp(-self.a1*self.internalTime)
        delta2 = delta*delta
        Om_1   = self.k2 - self.a1 + w*( (self.k1*w + f) / delta2 )
        term2  = (self.k1*w + f) / delta2 + w*Om_1

        # ------------------------------------------------
        # Euler step
        zd_dot  = self.derivs(self.zd,term2)
        self.zd = self.zd + zd_dot*dt
        #print np.linalg.norm(self.zd[:][0])
        

        # ------------------------------------------------
        # auxiliary control
        ua = ((self.k1*w + f) / delta2)*np.dot(self.J,self.zd) + Om_1*self.zd
        u  = ua - self.k2*z

        # -------------------------------------------------
        # transform to desired feedback law
        T = np.array([[(errorGlobal[x]*sin(currentPose[theta])-errorGlobal[y]*cos(currentPose[theta])), 1.],[1.,0.]])
        b = np.array([
            [desiredVelocity[0]*cos(errorGlobal[theta]) + desiredVelocity[1]*(errorGlobal[x]*sin(currentPose[theta])-errorGlobal[y]*cos(currentPose[theta]))],
            [desiredVelocity[1]]
            ])
        
        vfb = T.dot(u) + b

        self.internalTime += dt 

        return vfb[0][0], vfb[1][0], w, z



def normalizeAngle(angle):
    return np.mod((angle + np.pi), (2 * np.pi)) - np.pi


def angleDiff(angleA, angleB):
    """ Returns the angle needed to go from angleA to angleB, i.e. (angleB - angleA) """

    a = np.array([np.cos(angleA), np.sin(angleA)])
    b = np.array([np.cos(angleB), np.sin(angleB)])

    return np.arctan2(np.cross(a, b), np.dot(a, b))


if __name__ == "__main__":
    # execute only if run as a script
    dt = 0.02

    # desired speeds (vd will be increased within loop)
    vd, wd = 0., 0.

    # current robot position and heading
    xr, yr, thr = -.1, .1, 1.4
    v, w = 0., 0.

    UNL = Controller(zd= np.array([[1.4],[0.2]]))

    # first point along trajectory (zero velocity)
    xd, yd, thd = 0., 0., np.pi/2.

    # acceleration
    a = .05

    # initialize time
    time = 0.

    timeMax = 40.
    vMax = 0.5
    r = .25
    wMax = vMax / r

    T  = []
    Xr = []
    Yr = []
    Tr = []
    Xd = []
    Yd = []
    Td = []
    Vd = []
    Vrc = []
    Wrc = []
    Wd  = []
    controlW = []
    controlZ0 = []
    controlZ1 = []
    controlw = 0
    controlz = [0, 0]

    # main loop
    while time < timeMax:
    	T.append(time)
    	Xr.append(xr)
    	Yr.append(yr)
    	Tr.append(thr)
    	Vrc.append(v)
    	Wrc.append(w)
    	Xd.append(xd)
    	Yd.append(yd)
    	Td.append(thd)
    	Vd.append(vd)
        Wd.append(wd)
        controlW.append(controlw)
        controlZ0.append(controlz[0])
        controlZ1.append(controlz[1])
        '''
        # plot stuff
        plt.figure('')
        robot_circle = plt.Circle((xr, yr), r, color='b')
        goal_circle = plt.Circle((xd, yd), .1, color='r')
        plt.gca().add_artist(robot_circle)
        plt.gca().add_artist(goal_circle)
        # plt.hold(True)
        plt.plot(Xr, Yr, c='b', label='act')
        plt.plot(Xd, Yd, c='r', label='des')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.axis([-10, 10, 0, 20])
        plt.legend()
        plt.show()
        plt.pause(.0001)
        '''
    	if vd > vMax:
    		a = 0.

    	x, y, theta = 0, 1, 2
        thd = normalizeAngle(thd)
        thr = normalizeAngle(thr)
        errorGlobal = [
            xr - xd,
            yr - yd,
            # angleDiff(currentPose[theta], desiredPose[theta])
            normalizeAngle(angleDiff(thd, thr))
        ]
    	currentPose = [xr,yr,thr]
    	desiredVelocity = [vd,wd]
    	v,w, controlw, controlz = UNL.calculate(errorGlobal, currentPose, desiredVelocity, dt)
    	if v > vMax:
    		v = vMax
        # elif v <= 0:
        #     v = 0
    	elif v < -vMax:
    		v = -vMax

    	if w > wMax:
    		w = wMax
    	elif w < -wMax:
    		w = -wMax

    	# setup derivative of robot state at current velocity
    	# NOTE: this assumes infinite actuator bandwidth and no bias/scaling
    	# essentially, this means command = achieved
    	q = np.array([[xr],[yr],[thr]])
        S = np.array([[np.cos(thr), 0.],[np.sin(thr), 0.],[0.,1.]])
        V = np.array([[v],[w]])
        qdot = np.dot(S,V)

        # heun's step
        qpred = q + qdot*dt
        qpred = qpred.ravel()
        Spred = np.array([[np.cos(qpred[2]), 0.],[np.sin(qpred[2]), 0.],[0.,1.]])
        qdotpred = np.dot(Spred,V)
        q = q.ravel() + dt/2.*(qdot.ravel()+qdotpred.ravel())

        xr,yr,thr = q[0],q[1],q[2]

    	# increment to next path point and velocity at that point
    	yd   += vd*dt
    	vd   += a*dt

    	#print xd,yd

    	# increment time
    	time += dt

    # plot stuff
    plt.figure('path')
    # plt.hold(True)
    plt.plot(Xr, Yr, '--', label='act')
    plt.plot(Xd, Yd, label='des')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim([-2, 2])
    plt.legend()

    plt.figure('control output')
    plt.plot(T,Vrc,'b',label='vact')
    plt.plot(T,Vd,'b--',label='vdes')
    plt.plot(T,Wrc,'r',label='wact')
    plt.plot(T,Wd,'r--',label='wdes')
    plt.xlabel('Time')
    plt.ylabel('Velocities')
    plt.legend()

    '''
    plt.figure('error')
    ax = plt.subplot(3, 1, 1)
    plt.plot(T, Xr, label='x real')
    plt.plot(T, Xd, label='x des')
    plt.legend()
    ax = plt.subplot(3, 1, 2)
    plt.plot(T, Yr, label='y real')
    plt.plot(T, Yd, label='y des')
    plt.legend()
    ax = plt.subplot(3, 1, 3)
    plt.plot(T, Tr, label='theta real')
    plt.plot(T, Td, label='theta des')
    plt.legend()

    plt.figure('w/z error')
    ax = plt.subplot(3, 1, 1)
    plt.plot(T, controlW, label='w')
    plt.legend()
    ax = plt.subplot(3, 1, 2)
    plt.plot(T, controlZ0, label='z0')
    plt.legend()
    ax = plt.subplot(3, 1, 3)
    plt.plot(T, controlZ1, label='z1')
    plt.legend()
    '''
    plt.show()