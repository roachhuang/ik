# traj
#https://www.coursera.org/learn/robotics1/lecture/bNQfV/4-3-piepers-solution-1
from cmath import atan, isclose, pi, sqrt
from math import atan2
from sympy import Eq, Symbol, init_printing, solve, solveset, sin, cos, simplify, symbols, trigsimp

import numpy as np
import sympy as sp
import craig as cg

init_printing(use_unicode=True, use_latex='mathjax')
np.set_printoptions(precision=3, suppress=True)


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
                    print('q4, q5, q6:', t4 * 180 / pi, t5 * 180 / pi,
                          t6 * 180 / pi)
                    return ([t4, t5, t6])


'''
    for t1 in q4s:
        for t2 in q5s:
            for t3 in q6s:
                #myX = expr_x.subs([(q1, t1), (q2, t2), (q3, t3)])
                #myY = expr_y.subs([(q1, t1), (q2, t2), (q3, t3)])
                #myZ = expr_z.subs([(q2, t2), (q3, t3)])
                #if isclose(myX, x) and isclose(myY, y) and isclose(myZ, z):
                    #print('q1, q2, q3:', t1 * 180 / pi, t2 * 180 / pi,
                    #        t3 * 180 / pi)
                    #qs.append((t1, t2, t3))
'''


# aix 2&3 decide heigh of z and sqt(x**2+y**2)
# q4,q5, q6
# refer to must read UCLA ik.pdf for r6_3 )
# https://univ.deltamoocx.net/courses/course-v1:AT+AT_010_1102+2022_02_01/courseware/a3e573de127b85f1dcb23ea797cd253f/dc947a72e470ca516e9270c3bb4424e1/?child=first
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

    # Elur angle zyz
    # firstable rotate r4_3 in x axis to be in line w/ euler,
    alp3 = cg.dh_tbl[3, 0]
    r3prime_0 = r3_0 @ cg.Rx(alp3)
    r6_3prime = np.transpose(r3prime_0) @ r6_0
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


def pieper(t6_0):
    q1s = []
    q2s = []
    q3s = []
    #solved 3 axies in tuple
    qs = []

    u, q1, q2, q3 = sp.symbols('u,q1,q2,q3')

    (alp4, a4, d5) = cg.dh_tbl[4, :]
    (alp3, a3, d4) = cg.dh_tbl[3, :]
    (alp2, a2, d3) = cg.dh_tbl[2, :]
    (alp1, a1, d2) = cg.dh_tbl[1, :]

    print('t6_0', t6_0)
    #p4_0_org=P6_0_org
    (x, y, z) = t6_0[0:3, 3]
    print('4-0-org:', x, y, z)

    #(alp3, a3, d4)=cg.dh_tbl[3, :]
    #p4_3=np.array([a3, -d4*sin(alp3), d4*round(cos(alp3)),1])
    #print('p4_3:', p4_3)
    '''
    [x,y,z,1]t=P4_0_org=t1_0@t2_1@t3_2@p4_3_org
    =
    '''
    t4_3 = cg.get_ti2i_1(4)
    p4_3 = t4_3[:, 3]
    print('p4_3:', p4_3)
    t1_0 = cg.get_ti2i_1(1)
    t2_1 = cg.get_ti2i_1(2)
    t3_2 = cg.get_ti2i_1(3)

    #f is a func of theta3
    f = trigsimp(t3_2 @ p4_3)
    f1 = f[0]
    f2 = f[1]
    f3 = f[2]
    # g is a func of theta2 and theta3
    g = trigsimp(t2_1 @ f)
    g1 = trigsimp(g[0])
    g2 = trigsimp(g[1])
    g3 = trigsimp(g[2])
    #print('g1:', g1)
    #print('g2:', g2)
    #print('g3:', g3)
    '''
    p4_0=[x,y,z,1] = t1_0@g = [c1g1-s1g2, s1g1+c1g2, g3, 1]
    x=p4_0orgx=c1g1-s1g2, y=s1g1+c1g2
    frame0 to wrist 的長度平方 so as to eliminate q1
    c1g1**2-s1g2**2+s1g1**2+c1g2**2+g3**2=> g1**2(c1**2+s1**2)+g2**2(s1**2+c1**2)+g3**2
    hence, r=x**2+y**2+z**2=g1**2+g2**2+g3**2
    r is func of only theta2 & 3
    z = g3 (z is func of theta 3)
    '''

    r = x**2 + y**2 + z**2
    k1 = f1
    k2 = -f2
    k3 = trigsimp(f1**2 + f2**2 + f3**2 + a1**2 + d2**2 + 2 * d2 * f3)
    k4 = f3 * cos(alp1) + d2 * cos(alp1)

    #
    #r=trigsimp((k1*cos(q2)+k2*sin(q2))*2*a1 + k3)
    #z=trigsimp((k1*sin(q2)-k2*cos(q2))*sin(alp1) + k4)

    # watch out for the case when a+c=0, theta = 180
    if (a1 == 0.0):
        '''
        the distance from the origins of the 0 and 1 frames to the center of the spherical wrist
        (which is the origin of frames 4, 5, and 6) is a function of θ3 only;
        r=k3
        '''
        #q3s=solve(Eq(k3, x**2+y**2+z**2), q3)
        q3s = solve(Eq(r, k3), q3)
        rExpr = (k1 * sin(q2) - k2 * cos(q2)) * sin(alp1) + k4
        lExpr = z
    # general case when a1 != 0 & alpha1 != 0
    # combine k3 and g3
    elif (alp1 == 0.0):
        # z=k4
        q3s = solve(Eq(z, k4), q3)
        rExpr = (k1 * cos(q2) + k2 * sin(q2)) * 2 * a1 + k3
        lExpr = r
    else:
        '''
        r=(k1c2+k2s2)*2a1+k3
        z=(k1s2-k2c2)*sin(alp1)+k4
        r**2+z**2 to eliminate theta2, you can see from above 2 equations
        move k3&k4 to the left side, r and z 取平方和 to eliminate q2
        (r-k3/2*a1)**2 + (z-k4/sin(alp1))**2 = k1**2+k2**2
        =>(r-k3)**2/4*a1**2 + (z-k4)**2/sin(alp1)**2 = k1**2+k2**2
        '''
        #u=tan(q3/2)
        c3 = (1 - u**2) / (1 + u**2)
        s3 = (2 * u) / (1 + u**2)
        # very tricky - lexpr needs () for all expression!!!
        lExpr = (r - k3)**2 / (4 * (a1**2)) + ((z - k4)**2) / (sin(alp1))**2
        lExpr = sp.expand_trig(lExpr)
        #print('lExpr:', lExpr)
        lExpr = lExpr.subs([(sin(q3), s3), (cos(q3), c3)])
        #lExpr = lExpr.subs([(sin(q3), simplify('2*u/1+u**2')), (cos(q3), simplify('1-u**2/1+u**2'))])
        lExpr = lExpr.expand()
        #print('lExpr:', lExpr)
        rExpr = k1**2 + k2**2
        rExpr = sp.expand_trig(rExpr)
        #print('rExpr:', rExpr)
        rExpr = rExpr.subs([(sin(q3), s3), (cos(q3), c3)])
        rExpr = rExpr.expand()
        #print('rexpr:', rExpr)
        roots = solveset(Eq(lExpr, rExpr), u)
        # print('u:', roots)
        for root in roots:
            #q3s.append(2 * atan2(root, 1))
            q3s = np.append(q3s, 2 * atan2(root, 1))
            #rad = 2*atan(root)

        # this is for computing q2
        rExpr = (k1 * cos(q2) + k2 * sin(q2)) * 2 * a1 + k3
        lExpr = r

    #print('lExpr=r', lExpr, r)
    # remove duplicates from the list
    # solve q2: r=(k1*c2+k2*s2)*2*a1+k3
    q3s = list(dict.fromkeys(q3s))
    for t3 in q3s:
        #print('@t3:=', t3 * 180 / pi)
        #print('cos(t3)= {}'.format(cos(t3)))
        rExpr = rExpr.subs(q3, t3)
        tmp = solve(Eq(lExpr, rExpr), q2)
        q2s.extend(tmp)

    q2s = list(dict.fromkeys(q2s))
    #for t2 in q2s:
        #print('@t2:=', t2 * 180 / pi)

        # solve q1: x=c1*g1(q2,q3)-s1*g2(q2,q3)
    for t3 in q3s:
        #print('@t3:=', t3 * 180 / pi)
        for t2 in q2s:
            #print('@t2:=', t2 * 180 / pi)
            rExpr = cos(q1) * g1 - sin(q1) * g2
            rExpr = rExpr.subs([(q2, t2), (q3, t3)])
            tmp = solve(Eq(x, rExpr), q1)
            q1s.extend(tmp)

    q1s = list(dict.fromkeys(q1s))

    #verify ik
    #x=c1g1-s1g2, y=s1g1+c1g2
    expr_x = sp.expand_trig(cos(q1) * g1 - sin(q1) * g2)
    expr_y = sp.expand_trig(sin(q1) * g1 + cos(q1) * g2)
    expr_z = g3
    for t1 in q1s:
        for t2 in q2s:
            for t3 in q3s:
                myX = expr_x.subs([(q1, t1), (q2, t2), (q3, t3)])
                myY = expr_y.subs([(q1, t1), (q2, t2), (q3, t3)])
                myZ = expr_z.subs([(q2, t2), (q3, t3)])
                if isclose(myX, x) and isclose(myY, y) and isclose(myZ, z):
                    print('q1, q2, q3:', t1 * 180 / pi, t2 * 180 / pi,
                          t3 * 180 / pi)
                    #qs.append((t1, t2, t3))
                    #q1-3=np.append(qs, [t1,t2,t3])
                    # qs array contains verified q1~q3
                    qs = np.append(qs, [t1, t2, t3])
                    qs = np.append(qs, ik456(t6_0[0:3, 0:3], t1, t2, t3))
                    # get one verified q1-3 is enough
                    break
    print('q1-6:', qs)
    return qs

