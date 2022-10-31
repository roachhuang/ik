
# aix 2&3 decide heigh of z and sqt(x**2+y**2)
# q4,q5, q6
# refer to must read UCLA ik.pdf for r6_3 )
# https://univ.deltamoocx.net/courses/course-v1:AT+AT_010_1102+2022_02_01/courseware/a3e573de127b85f1dcb23ea797cd253f/dc947a72e470ca516e9270c3bb4424e1/?child=first
from cmath import atan, cos, sin, sqrt
from math import atan2
import numpy as np
import craig as cg

def ver456(r6_3, q4s, q5s, q6s):
    for t4 in q4s.real:
        for t5 in q5s.real:
            for t6 in q6s.real:
                t4_3 = cg.get_ti2i_1(4, t4)
                t5_4 = cg.get_ti2i_1(5, t5)
                t6_5 = cg.get_ti2i_1(6, t6)
                t6_3 = t4_3 @ t5_4 @ t6_5
                R6_3 = t6_3[0:3, 0:3]
                #print('r6_3:', r6_3)
                #print('R6_3:', R6_3)
                if np.allclose(r6_3.astype('float64'), R6_3.astype('float64')):
                    print('q4, q5, q6:', t4 * 180 / np.pi, t5 * 180 / np.pi,
                          t6 * 180 / np.pi)
                    return (np.array([t4, t5, t6], dtype=np.float64))

def ik456(r6_0, q1, q2, q3):
    q4s = []
    q5s = []
    q6s = []
    """
    q4,q5,q6=euler angles phi, theta, psi. see unit5 part3
    R6-3=Rz,phi Ry,theat Rz,psi = Rz,q4 Ry,q5 Rz,q6
    (q1,q2,q3)->R3-0->R6-3=(R3-0)t * R6-0->q4,q5,q6
    use standard DH table to compute R3-0
    1) R6-3=R4-3(q4)R5-4(q5)R6-5(q6) see unit5 part3
    2) also R6-3(q4,q5,q6) = (R3-0)t * R6-0
    for (1) = (2)

    """
    #after compute q1-3, we can know R3_0
    #T3_0=T1_0@T2_1@T3_2
    """
    r3_0=np.array([[c1c23, -c1s23, s1],
     [s1c23,    -s1s23,  -c1 ],
     [s23,  c23,  0   ])
    """
    t1 = q1
    t2 = q2
    t3 = q3
    # r3_0=X(alp0)Z(t1)X(alp1)Z(t2)X(alp2)Z(t3
    # = X(0)Z(58.61)X(-90)Z(-64.46)X(0)Z(-11.98)
    # this r3_0 uses standard dh table.
    r3_0 = np.array(
        [[cos(t1) * cos(t2 + t3), -cos(t1) * sin(t2 + t3), -sin(t1)],
         [sin(t1) * cos(t2 + t3), -sin(t1) * sin(t2 + t3),
          cos(t1)], [-sin(t2 + t3), -cos(t2 + t3), 0]])

    R3_0 = np.array([[0.6006, 0.7082, -0.3710], [0.24, 0.2830, 0.9286],
                     [0.7627, -0.6468, 0]])
    #print('r3_0', r3_0)
    r6_3 = np.transpose(r3_0) @ r6_0
    #print('r6-3', r6_3)

    # Elur angle zyz, NTU
    # firstable rotate r4_3 in x axis to be in line w/ euler,
    alp3 = cg.dh_tbl[3, 0]
    r3prime_0 = r3_0 @ cg.Rx(alp3)
    r6_3prime = r3prime_0.T @ r6_0
    #print('r6_3prime:', r6_3prime)
    # r6_3prime = r4_3prime_z_y_z(alp, beta, gama)
    r13 = r6_3prime[0, 2]
    r23 = r6_3prime[1, 2]
    r31 = r6_3prime[2, 0]
    r32 = r6_3prime[2, 1]
    r33 = r6_3prime[2, 2]
    beta = atan2(sqrt(r31**2 + r32**2).real, r33)
    alpa = atan2(r23 / sin(beta), r13 / sin(beta))
    gama = atan2(r32 / sin(beta), -r31 / sin(beta))

    # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/
    R = r6_3prime
    eul1 = atan2(R.item(1,2),R.item(0,2))
    sp = sin(eul1)
    cp = cos(eul1)
    eul2 = atan2(cp*R.item(0,2)+sp*R.item(1,2), R.item(2,2))
    eul3 = atan2(-sp*R.item(0,0)+cp*R.item(1,0),-sp*R.item(0,1)+cp*R.item(1,1))
    print("phi =", eul1)
    print("theta =", eul2)
    print("psi =", eul3)
    '''
    r6_3=np.array([[c4c5c6-s4s6, -c4c5s6-s4c6, -c4s5],
                [s5c6, -s5s6, -c5],
                [-s4c5c6-c4s6, s4c5s6-c4c6, s4s5]])
    '''
    c4s5 = -r6_3[0, 2]
    s4s5 = r6_3[2, 2]
    c5 = -r6_3[1, 2]
    s5c6 = r6_3[1, 0]
    s5s6 = -r6_3[1, 1]
    print('c5', c5)

    # s5!=0, i.e q5 !=0 or pi, s5=+-sqrt(1-c5**2)
    q5_1 = atan(sqrt(1 - c5**2) / c5)
    q5s = np.append(q5s, q5_1)
    #print('q5_1:', q5_1 * 180 / pi)

    q4_1 = atan2(s4s5, c4s5)
    q4s = np.append(q4s, q4_1)
    #print('q4_1:', np.rad2deg(q4_1))

    q6_1 = atan2(s5s6, s5c6)
    q6s = np.append(q6s, q6_1)
    #print('q6_1', np.rad2deg(q6_1))

    # s5 < 0, 0 ~ -pi
    q5_2 = atan(-sqrt(1 - c5**2) / c5)
    q5s = np.append(q5s, q5_2)
    #print('q5_2', q5_2 * 180 / pi)

    # coz s5 is neg
    q4_2 = atan2(-s4s5, -c4s5)
    q4s = np.append(q4s, q4_2)
    #print('q4_2', np.rad2deg(q4_2))

    q6_2 = atan2(-s5s6, -s5c6)
    q6s = np.append(q6s, q6_2)
    #print('q6_2', np.rad2deg(q6_2))

    #DH6_3 = np.array([[c(q4-q6), sin(q4-q6),0, 0],
    # sin(q4-q6), -cos(q4-q6),0, 0],
    # [0,0,-1, d4],
    # [0,0,0,1])
    # R6_3=np.transpose(R3_0)@R6_0
    return ver456(r6_3, q4s, q5s, q6s)
