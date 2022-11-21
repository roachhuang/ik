from math import atan2, sqrt
import craig as cg
import numpy as np
import pieper as pp
import pandas as pd
import matplotlib.pyplot as plt

def main():

    np.set_printoptions(precision=4, suppress=True)
    dh_tbl = np.array([[0, 0, 0], [np.deg2rad(-90), -30, 0], [0, 340, 0],
                       [np.deg2rad(-90), -40, 338], [np.deg2rad(90), 0, 0],
                       [np.deg2rad(-90), 0, 0]])

    cg.setDhTbl(dh_tbl)

    # 從機械手臂的Frame {0}座標系來看，杯子的中心（Frame {C}原點）在不同時間點的位置及姿態分別在下表列出。
    # time, x, y, z, tx, ty, tz(wrt world frame)
    p = np.array([
        [0, 550, 270, 19.5, 0, 0, 35],
        [2, 550, 270, 79.5, 0, 0, 35],  # viapoint1
        [6, 330, 372, 367, 0, -60, 0],  # viapoint2
        [9, 330, 472, 367, 0, -60, 0],
    ])

    # duration
    #print (t)
    col_names = ['ti', 'xi', 'yi', 'zi', 'qx', 'qy', 'qz']
    row_names = ['p0', 'p1', 'p2', 'pf']
    P = pd.DataFrame(p, columns=col_names, index=row_names)
    print(P)
    print(' ')

    Tcup_6 = np.array([[0, 0, 1, 0], [0, -1, 0, 0], [1, 0, 0, 206],
                       [0, 0, 0, 1]])
    tc_0 = []
    # pg 6, compute each point's tc_0 transformation matrix based on cartesion space. range: 0~totalPoints-1
    tf = np.empty(shape=[4, 4])  #p0 to pf
    t6_0 = []
    totalPoints, num_cols = p.shape
    segs = totalPoints - 1
    # substract 頭, 尾
    viaPoints = totalPoints - 2
    DOF = 6

    for i in range(totalPoints):
        tx, ty, tz = np.deg2rad(p[i, 4:7])
        # combine rot and translation vector into transformation matrix
        tf[:3, :3] = cg.Rx(tx) @ cg.Ry(ty) @ cg.Rz(tz)  # rotation matrix
        tf[:3, 3] = p[i, 1:4]  # x,y,z
        tf[3, :] = [0, 0, 0, 1]
        tc_0.append(tf)
        # get t6_0 for each point
        t6_0 = tc_0[i] @ np.linalg.inv(Tcup_6)
        # replace p with ik's result - thetas
        p[i, 1:7] = np.rad2deg(pp.pieper(t6_0))

    col_names = ['ti', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6']
    P = pd.DataFrame(p, columns=col_names, index=row_names)
    print(P.round(2))
    print(' ')

    # 開始規劃 trajectory
    t = np.diff(p, axis=0)
    v1s = np.array([])
    v2s = np.array([])
    v3s = np.array([])

    # 對P6_0 在各點的pos and 姿態, compute vel, 每段parabolic func 區間長0.5s
    durationOfPara = 0.5
    for col in range(1, 7):
        v1s = np.append(v1s, t[0, col] / (t[0, 0] - durationOfPara / 2))
        v2s = np.append(v2s, t[1, col] / (t[1, 0]))
        v3s = np.append(v3s, t[2, col] / (t[2, 0] - durationOfPara / 2))

    v = np.array([[0, 0, 0, 0, 0, 0], v1s, v2s, v3s, [0, 0, 0, 0, 0, 0]])
    # np.asarray(v)
    row_names = ['v0', 'v1', 'v2', 'v3', 'vf']
    col_names = ['q1', 'q2', 'q3', 'q4', 'q5', 'q6']
    V = pd.DataFrame(v, columns=col_names, index=row_names)
    print(V.round(2))

    a = np.diff(v, axis=0) / durationOfPara
    row_names = ['a0', 'a1', 'a2', 'af']
    A = pd.DataFrame(a, columns=col_names, index=row_names)
    print(A)

    ts=p[:, 0]

    # in 0, 0.5s. col[0~2]: x, y, z
    def eq1(t, col):
        dt=t-0
        v0=v[0, col]
        a0=a[0, col]
        return p[0,col+1]+v0*dt+1/2*a0*dt**2

    # in 0.5, 1.75
    def eq2(t, col):
        dt=t-0.25
        v1=v[1,col]
        return p[0, col+1]+v1*dt

    # in 1.7, 2.25
    def eq3(t, col):
        v1=v[1,col]
        a1=a[1,col]
        dt1=t-0.25
        dt2=t-(ts[1]-0.25)
        return p[0,col+1]+v1*dt1+1/2*a1*dt2**2

    # in 2.25, 5.75
    def eq4(t, col):
        dt=t-ts[1]
        v2=v[2,col]
        return p[1,col+1]+v2*dt

    # 5.75, 6.25
    def eq5(t, col):
        dt1=t-ts[1]
        dt2=t-(ts[2]-0.25)
        v2=v[2,col]
        a2=a[2,col]
        return p[1,col+1]+v2*dt1+1/2*a2*dt2**2

    # 6.25, 8.5
    def eq6(t, col):
        dt=t-ts[2]
        v3=v[3,col]
        return p[2,col+1]+v3*dt

    # col - 0:x, 1:y, 2:theta
    # 8.5 ~9s
    def eq7(t, col):
        dt1=t-ts[2]
        dt2=t-(ts[totalPoints-1]-0.5)
        v3=v[3,col]
        a3=a[3,col]
        return p[2,col+1]+v3*dt1+1/2*a3*dt2**2

    # 0 ~ 9s
    timeAxis = np.arange(0.0, 9.0, 0.1)
    inputPoints=[]
    plt.xlabel('Time')
    # col - 0~2, denote x, y or theta data
    arrT = np.arange(0, 9, 0.1, dtype=int)
    for col in range(3):
        inputPoints=[]
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

        # plt.plot(array_x, array_y, marker='o', color='crimson', linestyle='-')
        ax = plt.axes(projection='3d')
        ax.plot3D(inputPoints[0], inputPoints[1], inputPoints[2], 'gray')
        plt.show()

    # plot 建立並繪出各DOF 在每個時間區段軌跡
    # linear/parabolic 共7段 （每段parabolic curve 時間設定為0.5s）

    # ik p0 ~ pf 所有的點
    # FK to xyz space
    # plt simulation

if __name__ == "__main__":
    main()
