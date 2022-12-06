# traj
#https://www.coursera.org/learn/robotics1/lecture/bNQfV/4-3-piepers-solution-1
# from cmath import isclose
from math import atan2, pi, isclose
from sympy import Eq, Symbol, init_printing, solve, solveset, sin, cos, simplify, symbols, trigsimp

import ik456 as ik
import numpy as np
import sympy as sp
import craig as cg

init_printing(use_unicode=True, use_latex='mathjax')

def pieper(t6_0):
    np.set_printoptions(precision=4, suppress=True)
    q1s = []
    q2s = []
    q3s = []
    #solved 3 axies in tuple
    qs = []

    u, q1, q2, q3 = symbols('u,q1,q2,q3')

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
    '''
    t4_3 = cg.get_ti2i_1(4)
    p4_3 = t4_3[:, 3]
    print('p4_3:', p4_3)
    t1_0 = cg.get_ti2i_1(1)
    t2_1 = cg.get_ti2i_1(2)
    t3_2 = cg.get_ti2i_1(3)

    #f is a func of theta3, trigsimp her may not need, remove it later
    # f = trigsimp(t3_2 @ p4_3)
    f = t3_2 @ p4_3

    f1, f2, f3 = f[0:3]
    # g is a func of theta2 and theta3, trigsimp in her is essential!!!
    g = trigsimp(t2_1 @ f)
    g1, g2, g3 = g[0:3]

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
    q3s = [*set(q3s)]
    for t3 in q3s:
        #print('@t3:=', t3 * 180 / pi)
        #print('cos(t3)= {}'.format(cos(t3)))
        r = rExpr.subs(q3, t3)
        tmp = solve(Eq(lExpr, r), q2)
        q2s.extend(tmp)

    #q2s = list(dict.fromkeys(q2s))
    q2s = [*set(q2s)]
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

    # q1s = list(dict.fromkeys(q1s))
    q1s = [*set(q1s)]

    # q1 is simply atan2(y,x) according to ntu
    # q1s.append(atan2(y,x))

    #verify ik
    #x=c1g1-s1g2, y=s1g1+c1g2
    # expr_x = sp.expand_trig(cos(q1) * g1 - sin(q1) * g2)
    # expr_y = sp.expand_trig(sin(q1) * g1 + cos(q1) * g2)
    expr_x = cos(q1) * g1 - sin(q1) * g2
    expr_y = sin(q1) * g1 + cos(q1) * g2
    expr_z = g3
    for t1 in q1s:
        for t2 in q2s:
            for t3 in q3s:
                myX = expr_x.subs([(q1, t1), (q2, t2), (q3, t3)])
                myY = expr_y.subs([(q1, t1), (q2, t2), (q3, t3)])
                myZ = expr_z.subs([(q2, t2), (q3, t3)])
                if isclose(myX, x) and isclose(myY, y) and isclose(myZ, z):
                    print(f'q1: {t1 * 180 / pi}, q2: {t2 * 180 / pi}, q3: {t3 * 180 / pi}')
                    #qs.append((t1, t2, t3))
                    #q1-3=np.append(qs, [t1,t2,t3])
                    # qs array contains verified q1~q3
                    #qs = np.append(qs, [t1, t2, t3])
                    q123 = np.array([t1, t2, t3],dtype=np.float)
                    q456 = ik.ik456(t6_0[0:3, 0:3], t1, t2, t3)
                    # q456 = ik.ik4_5_6(t6_0[0:3, 0:3], t1, t2, t3)
                    qs = np.concatenate((q123, q456))
                    #qs = np.append(qs, ik.ik456(t6_0[0:3, 0:3], t1, t2, t3))
                    # get one verified q1-3 is enough
                    print(f'q1-6: {np.rad2deg(qs)}')
                    return qs
    return 'no match'
