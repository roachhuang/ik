import craig as cg
from math import radians
import numpy as np
import pieper as pp
from sympy import Matrix, Symbol, Eq, solve, sin, cos

def main():
    dh_tbl=np.array([[0,0,0],
            [np.deg2rad(-90), -30, 0],
            [0, 340, 0],
            [np.deg2rad(-90), -40, 338],
            [np.deg2rad(90), 0, 0],
            [np.deg2rad(-90),0,0]])

    cg.setDhTbl(dh_tbl)

    viaPoint = 2
    # time, x, y, z, tx, ty, tz
    C = np.array([
        [0, 550, 270, 19.5, 0, 0, 35],
        [2, 550, 270, 79.5, 0, 0, 35],
        [6, 330, 372, 367, 0, -60, 0],
        [9, 330, 472, 367, 0, -60, 0],
    ])
    Tcup_6 = np.array([[0, 0, 1, 0], [0, -1, 0, 0], [1, 0, 0, 206], [0, 0, 0, 1]])

    dt1 = C[1, 1] - C[0, 1]
    dt2 = C[2, 1] - C[1, 1]
    dt3 = C[3, 1] - C[2, 1]
    T = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, dt1, dt1**2, dt1**3, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, dt2, dt2**2, dt2**3, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, dt3, dt3**2, dt3**3],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2 * dt3, 3 * (dt3**2)],
                [0, 1, 2 * dt1, 3 * (dt1**2), 0, -1, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 6 * dt1, 0, 0, -2, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 2 * dt2, 3 * (dt2**2), 0, -1, 0, 0],
                [0, 0, 0, 0, 0, 0, 2, 6 * dt2, 0, 0, -2, 0]])

    if (np.linalg.det(T)):
        T_inv = np.linalg.inv(T)
    else:
        print('error det !=0; no inverse')
    #print(T_inv)
    #print(np.matmul(T_inv, C_of_X))

    #C_of_Q = np.ndarray(shape=(0, 3), dtype=float)
    Q12x3 = []
    num_rows, num_cols=C.shape
    ty = C[1, 6]
    tz = C[3, 5]
    (x, y, z) = C[0, 1:4]
    tc_0_0 = np.array([[cos(tz), -sin(tz), 0, x], [sin(tz), cos(tz), 0, y],
                    [0, 0, 0, z], [0, 0, 0, 1]])

    # @2s rotate in z aix
    (x, y, z) = C[1, 1:4]
    tc_0_2 = np.array([[cos(tz), -sin(tz), 0, x], [sin(tz), cos(tz), 0, y],
                    [0, 0, 0, z], [0, 0, 0, 1]])

    # @6s rotate in y aix
    (x, y, z) = C[2, 1:4]
    tc_0_6 = np.array([[cos(ty), 0, sin(ty), x], [0, 1, 0, y],
                    [-sin(ty), 0, cos(ty), z], [0, 0, 0, 1]])
    # @9s
    (x, y, z) = C[3, 1:4]
    tc_0_9 = np.array([[cos(ty), 0, sin(ty), x], [0, 1, 0, y],
                    [-sin(ty), 0, cos(ty), z], [0, 0, 0, 1]])

    # convert cartesian space to joint space: 3 segments 0-2s 2-4s 4-9s
    t6_0 = tc_0_0 @ np.linalg.inv(Tcup_6)
    pp.pieper(t6_0)
    # Q12x3 = np.append(Q12x3, q123)
    t6_0 = tc_0_2 @ np.linalg.inv(Tcup_6)
    pp.pieper(t6_0)
    t6_0 = tc_0_6 @ np.linalg.inv(Tcup_6)
    pp.pieper(t6_0)
    t6_0 = tc_0_9 @ np.linalg.inv(Tcup_6)
    pp.pieper(t6_0)

    # print (Q12x1)
    #padding 0 for last 6 rows
    # col:theta1~theat3 of the arm
    Q12x3.resize((12, 3), refcheck=False)
    print('theta:', Q12x3.real)

    # 2 dicimal points accroding to 5-6 exm2
    np.set_printoptions(precision=3, suppress=True)
    A12x3 = (np.round(T_inv @ Q12x3, 2))

    # make a0~a3 in one row, 0~1s,2~4s,4~9s in one group(3x4), each group represent one theta
    # 3xseg, 3xtheta, 4xai
    A12x3 = np.transpose(A12x3).reshape(3, 3, 4)
    #A = np.reshape(A, (3, 3, 4))
    # A is now a 3-D array
    # row: theta, col:aij,
    print('A12x3:', A12x3.real)

if __name__ == "__main__":
    main()
    