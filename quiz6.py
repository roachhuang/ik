# cartesion space
import numpy as np
import pandas as pd

np.set_printoptions(precision=2, suppress=True)

col_names = ['ti', 'xi', 'yi', 'qi']
# 4 points has 3 linear lines and 4個2次段
traj = np.array([[0, -4, 0, 120], [2, -5, 5, 45], [4, 2, 3, 30], [9, 5, -3, 0]])
# t = traj[:, 0]
t = np.diff(traj,axis=0)
v1s =np.array([])
v2s =np.array([])
v3s =np.array([])
# duration
#print (t)
TRAJ = pd.DataFrame(traj, columns=col_names)
print(TRAJ)
# 每段parabolic func 區間長0.5s
durationOfPara = 0.5
# vel = ds/dt
for col in range(1, 4):
    v1s = np.append(v1s, t[0,col] / (t[0, 0] - durationOfPara / 2))
    v2s = np.append(v2s, t[1,col] / (t[1, 0]))
    v3s = np.append(v3s, t[2,col] / (t[2, 0] - durationOfPara / 2))
v = np.array([[0,0,0],v1s,v2s,v3s,[0,0,0]])
# np.asarray(v)
col_names = ['x', 'y', 'q (deg/s)']
row_names=['v0','v1','v2','v3','vf']
V = pd.DataFrame(v, columns=col_names,  index=row_names)
print(V)
a = np.diff(v,axis=0)/durationOfPara
row_names=['a0','a1','a2','af']
A = pd.DataFrame(a, columns=col_names,  index=row_names)
print(A)
'''
linear: q(t)=qj+q'jk*dt=q'jk*(t-tqj)
para:   q(t)=qj+q'jk*dt1+1/2*q"k*dt2**2=qj+q'jk*(t-tqj)+1/2*q"k(t-tqk+1/2tk)**2
xeq1(t)=x0+v0dt+1/2*a0*dt**2 = -4+0(t-0)+1/2*4.57(t-0)**2               (t in 0, 0.5)
xeq2(t)=x0+v1*dt = -4+2.29(t-0.25)                                      (t in 0.5, 1.75)

xeq3(t)=x0+v1*dt1+1/2*a1*dt2**2 = -4+2.29(t-0.25)+1/2*(-1.57)(t-1.75)**2 (t in 1.75, 2.25)
xeq4(t)=x1+v2*dt = 0+1.5(t-2)                                           (t in 2.25, 3.75)

xeq5(t)=x1+v2*dt1+1/2*a2*dt2**2 = 0+1.5(t-2)+1/2*(-2.27)(t-3.75)**2     (t in 3.75, 4.25)
xeq6(t)=x2+v3*dt = 3+0.35(t-4)                                          (t in 4.25, 6.5)

xeq7(t)=x2+v3dt1+1/2*a3*dt2**2 = 3+0.36(t-4)+1/2*(-0.73)(t-6.5)**2      (t in 6.5, 7)
pick some points to run ik
'''
#
def eq1(t):
    dt=t-0
    for i in range(3):
        v0=v[0,i]
        a0=a[0,i]
        pos=traj[0,i+1]+v0*dt+1/2*a0*dt**2
        print (pos)

def eq2(t):
    dt=t-0.25
    for i in range(3):
        v1=v[1,i]
        pos=traj[0,i+1]+v1*dt
        print (pos)

def eq3(t):
    for i in range(3):
        v1=v[1,i]
        a1=a[1,i]
        dt1=t-0.25
        dt2=t-1.75
        pos=traj[0,i+1]+v1*dt1+1/2*a1*dt2**2
        print (pos)

def eq4(t):
    dt=t-2
    for i in range(3):
        v2=v[2,i]
        pos=traj[1,i+1]+v2*dt
        print (pos)

def eq5(t):
    dt1=t-2
    dt2=t-3.75
    for i in range(3):
        v2=v[2,i]
        a2=a[2,i]
        pos = traj[1,i+1]+v2*dt1+1/2*a2*dt2**2
        print ('t=4', pos)

def eq6(t):
    dt=t-4
    for i in range(3):
        v3=v[3,i]
        pos=traj[2,i+1]+v3*dt
        print (pos)

def eq7(t):
    dt1=t-4
    dt2=t-8.5
    for i in range(3):
        v3=v[3,i]
        a3=a[3,i]
        pos = traj[2,i+1]+v3*dt1+1/2*a3*dt2**2
        print (pos)


#eq1(0)
#eq2(1)
#eq3(2)
#eq4(3)
eq5(4)
#eq6(8)
#eq7(9)

'''
for t in 9s: every 0.25s
    switch(t):
        case (t in 0, 0.5)
            set a0 to motor
        case (t in 0.5, 1.75)
            set v1 to motor
        case (t in 1.75, 2.25)
            set a1 to motor
        case (t in 2.25, 3.75)
            set v2 to motor

import time
def setInterval(func, sec):
    time.sleep(sec)
    func()
    setInterval(func(), sec)

def call():
    print('hello!')

setInterval(call, 3)
'''