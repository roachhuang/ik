# cartesion space
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
totalPoints, num_cols = traj.shape
segs = totalPoints - 1
# substract 頭, 尾
viaPoints = totalPoints - 2
ts=traj[:, 0]
# in 0, 0.5s
def eq1(t, col):
    dt=t-0
    v0=v[0, col]
    a0=a[0, col]
    return traj[0,col+1]+v0*dt+1/2*a0*dt**2

# in 0.5, 1.75
def eq2(t, col):
    dt=t-0.25
    v1=v[1,col]
    return traj[0, col+1]+v1*dt

# in 1.7, 2.25
def eq3(t, col):
    v1=v[1,col]
    a1=a[1,col]
    dt1=t-0.25
    dt2=t-(ts[1]-0.25)
    return traj[0,col+1]+v1*dt1+1/2*a1*dt2**2

# in 2.25, 3.75
def eq4(t, col):
    dt=t-ts[1]
    v2=v[2,col]
    return traj[1,col+1]+v2*dt

# 3.75, 4.25
def eq5(t, col):
    dt1=t-ts[1]
    dt2=t-(ts[2]-0.25)
    v2=v[2,col]
    a2=a[2,col]
    return traj[1,col+1]+v2*dt1+1/2*a2*dt2**2

# 4.25, 6.4
def eq6(t, col):
    dt=t-ts[2]
    v3=v[3,col]
    return traj[2,col+1]+v3*dt

# col - 0:x, 1:y, 2:theta
def eq7(t, col):
    dt1=t-ts[2]
    dt2=t-(ts[totalPoints-1]-0.5)
    v3=v[3,col]
    a3=a[3,col]
    return traj[2,col+1]+v3*dt1+1/2*a3*dt2**2

# 0 ~ 9s
timeAxis = np.arange(0.0, 9.0, 0.1)
inputPoints=[]
inputPoints.append([])
inputPoints.append([])
inputPoints.append([])
plt.xlabel('Time')
# col - 0~2, denote x, y or theta data
for col in range(3):
    plt.ylabel(str(col))
    for t in timeAxis:
        if t >= ts[0] and t <= ts[0]+0.5:
            inputPoints[col].append(eq1(t, col))
        elif t > ts[0]+0.5 and t <= ts[1]-0.25:
            inputPoints[col].append(eq2(t, col))
        elif t > ts[1]-0.25 and t <= ts[1]+0.25:
            inputPoints[col].append(eq3(t, col))
        elif t > ts[1]+0.25 and t <= ts[2]-0.25:
            inputPoints[col].append(eq4(t, col))
        elif t > ts[2]-0.25 and t <= ts[2]+0.25:
            inputPoints[col].append(eq5(t, col))
        elif t > ts[2]+0.25 and t <= ts[totalPoints-1]-0.5:
            inputPoints[col].append(eq6(t, col))
        elif t > ts[totalPoints-1]-0.5 and t <= ts[totalPoints-1]:
            inputPoints[col].append(eq7(t, col))

    #array_x = [x for x, y in inputPoints]   # time

    # plt.plot(array_x, array_y, marker='o', color='crimson', linestyle='-')
    # this fig has 1 row, 3 col in one page
    plt.subplot(1, 3, col+1)
    plt.plot(timeAxis, inputPoints[col], 'r')
    plt.grid()

plt.show()
    # plt.pause(5)

ax = plt.axes(projection='3d')
ax.plot3D(inputPoints[0], inputPoints[1], timeAxis, 'r')
plt.show()

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