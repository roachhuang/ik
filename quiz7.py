import craig as cg
import numpy as np
import pieper as pp
import pandas as pd

# SPACE='cartesion'
SPACE = 'joint'
# https://arduinogetstarted.com/faq/how-to-control-speed-of-servo-motor

def main():
    #np.set_printoptions(precision=2, suppress=True)
    dh_tbl = np.array([[0, 0, 0], [np.deg2rad(-90), -30, 0], [0, 340, 0],
                       [np.deg2rad(-90), -40, 338], [np.deg2rad(90), 0, 0],
                       [np.deg2rad(-90), 0, 0]])

    cg.setDhTbl(dh_tbl)

    # 從機械手臂的Frame {0}座標系來看，杯子的中心（Frame {C}原點）在不同時間點的位置及姿態分別在下表列出。
    # time, x, y, z, tx, ty, tz
    p = np.array([[0, 630, 364, 20, 0, 0, 0],
                  [3, 630, 304, 220, 60, 0, 0],
                  [7, 630, 220, 24, 180, 0, 0]],
                 dtype=float)

    col_names = ['ti', 'xi', 'yi', 'zi', 'qx', 'qy', 'qz']
    row_names = ['p0', 'p1', 'p2']
    P = pd.DataFrame(p, columns=col_names, index=row_names)
    print(P)
    print('-------------------------------------------------')

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
        if SPACE == 'cartesion':
            col_names = ['ti', 'xi', 'yi', 'zi', 'qx', 'qy', 'qz']
            p[i, 4:7] = cg.rotationMatrixToEulerAngles(t6_0)
            p[i, 1:4] = t6_0[0:3, 3].T
        else:
            col_names = ['ti', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6']
            p[i, 1:7] = np.rad2deg(pp.pieper(t6_0))

    P = pd.DataFrame(p, columns=col_names, index=row_names)
    print(P)
    print('-------------------------------------- ')

    # 開始規劃 trajectory
    t = np.diff(p, axis=0)
    # segs + 2(head and tail) = no. of v
    v1s = np.array([])
    v2s = np.array([])
    v3s = np.array([])

    # 對P6_0 在各點的pos and 姿態, compute vel, 每段parabolic func 區間長0.5s
    durationOfPara = 0.5
    for col in range(1, 7):
        v1s = np.append(v1s, t[0, col] / (t[0, 0] - durationOfPara / 2))
        v2s = np.append(v2s, t[1, col] / (t[1, 0] - durationOfPara / 2))
        #v3s = np.append(v3s, t[2, col] / (t[2, 0] - durationOfPara / 2))
    v = np.array([[0, 0, 0, 0, 0, 0], v1s, v2s, [0, 0, 0, 0, 0, 0]])
    # np.asarray(v)
    if SPACE == 'cartesion':
        col_names = ['xi', 'yi', 'zi', 'qx', 'qy', 'qz']
    else:
        col_names = ['q1', 'q2', 'q3', 'q4', 'q5', 'q6']

    # time 0, 3, 7s
    row_names = ['v0', 'v1', 'v2', 'vf']
    V = pd.DataFrame(v, columns=col_names, index=row_names)
    # v1: 0.5~2.75, v2: 3.75~6.5. linear segs
    print(V)

    # parabolic segs
    a = np.diff(v, axis=0) / durationOfPara
    row_names = ['a0', 'a1', 'af']
    A = pd.DataFrame(a, columns=col_names, index=row_names)
    # a0: 0~0.5s, a1:2.75~3.25, af: 6.5~7s
    print(A)

    # plot 建立並繪出各DOF 在每個時間區段軌跡, x,y,z to time
    # plot thetas to time for the 6 axes （以此theats 對 time 的關係來control motors)
    # linear/parabolic 共7段 （每段parabolic curve 時間設定為0.5s）

    # ik p0 ~ pf 所有的點
    # FK to xyz space to verify
    # plt simulation

    # xeq4(t)=x1+v2*dt = 0+1.5(t-2)
    dt = 5 - 4
    for i in range(3):
        v2 = v[2, i]
        pos = p[1, i + 1] + v2 * dt
        print(pos)

if __name__ == "__main__":
    main()
