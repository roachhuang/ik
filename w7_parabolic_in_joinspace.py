from math import atan2, sqrt
import craig as cg
import numpy as np
import pieper as pp
import pandas as pd
import plan_traj as pt

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
        tf[:3, :3] = cg.Rot('x', tx) @ cg.Rot('y', ty) @ cg.Rot('z', tz)  # rotation matrix
        tf[:3, 3] = p[i, 1:4]  # x,y,z
        tf[3, :] = [0, 0, 0, 1]
        tc_0.append(tf)
        # get t6_0 for each point
        t6_0 = tc_0[i] @ np.linalg.inv(Tcup_6)
        # replace p with ik's result - thetas
        p[i, 1:7] = np.rad2deg(pp.pieper(t6_0))
        # fk to see if q1-q6 computed by ik are correct
        # +0.0 to fix -0 in array, decimals=1 for fixing allclose returns false if decimals=2
        fk_t6_0 = np.around(cg.fk_6axes(np.deg2rad(p[i,1:7])), decimals=1)+0.0
        print (f'fk_t6_0: {fk_t6_0}')
        assert np.allclose(np.around(t6_0, decimals=1), fk_t6_0)

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

    pt.planTraj(p)

    # plot 建立並繪出各DOF 在每個時間區段軌跡
    # linear/parabolic 共7段 （每段parabolic curve 時間設定為0.5s）

    # ik p0 ~ pf 所有的點
    # FK to xyz space
    # plt simulation


if __name__ == "__main__":
    main()
