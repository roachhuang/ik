from cmath import acos, atan, pi, sqrt
from math import radians
import craig as cg
import sympy as sp
from sympy import Symbol, init_printing, solve, symbols, trigsimp, simplify
import numpy as np

print (np.version.version)

np.set_printoptions(precision=3, suppress=True)
init_printing(use_unicode=True)
# init_printing( use_latex='mathjax' )  # use pretty math output
q1, q2, q3=sp.symbols('q1,q2,q3')

dh_tbl=np.array([[0,0,0],
                    [np.deg2rad(-90), 0, 220],
                    [0, 430, -90],
                    [np.deg2rad(-90), 0, 430],
                    [np.deg2rad(90), 0, 0],
                    [np.deg2rad(-90),0,0]])

cg.setDhTbl(dh_tbl)
t1=15
t2=-40
t3=-50
t4=30
t5=70
t6=25

def quiz4():
    t1_0 = cg.get_ti2i_1(1, np.deg2rad(t1))
    t2_1 = cg.get_ti2i_1(2)
    t2_1 = t2_1.subs(q2, t2)
    print('t21 subs:', sp.sympify(t2_1, evaluate=True))
    t3_2 = cg.get_ti2i_1(3, np.deg2rad(t3))
    t4_3 = cg.get_ti2i_1(4, np.deg2rad(t4))
    t5_4 = cg.get_ti2i_1(5, np.deg2rad(t5))
    t6_5 = cg.get_ti2i_1(6, np.deg2rad(t6))

    t4_0 = sp.trigsimp(t1_0*t2_1*t3_2*t4_3)
    print ('t4_0:', t4_0)
    #t40 = np.round(t4_0.astype(np.double), decimals=2)
    p4_0 = t4_0[:, 3]
    print('p4_0:', p4_0)

    t6_0=t4_0@t5_4@t6_5
    print ('t6_0:', t6_0)
    #print ('t6_0:', np.round(t6_0.astype(np.double),2))
    #ans: 0//220//0.64//-0.77

#ans to quiz 2: 0//430//-90//90//0//0

def quiz3():
    t2_1 = cg.get_ti2i_1(2, radians(t2))
    print (type(t2_1))
    #print ('t2_1:', t2_1)
    #print ('t2_1:', np.round(t2_1.astype(np.double),2))
    #ans: 0//220//0.64//-0.77

quiz3()

quiz4()