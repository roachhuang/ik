from cmath import acos
from math import atan2, pi, sqrt
# the reason for importing cos and sin from sympy is Rot case (can't import them from math, plz note!!!)
from sympy import trigsimp, Symbol, init_printing, sin, cos, symbols, simplify
import numpy as np
#import sympy as sp
from math import log10, floor
import pandas as pd

np.set_printoptions(precision=2, suppress=True)

# 取有效數字 sig 位, 科學記號表示時就先取有效數字兩位,最後乘上科學記號的冪次
def my_sigfig(x, sig=2):
    return round(x, sig - int(floor(log10(abs(x)))) - 1)

#dh for quiz4
dh_tbl = []

def setDhTbl(dh):
    global dh_tbl
    dh_tbl = dh
    col_names = ['alphai-1', 'ai-1', 'di']
    DH = pd.DataFrame(dh, columns=col_names)
    print(DH)
    print('------------------------------------ ')


def Rot(axis, rad):
    theta = Symbol('theta')
    theta = rad

    if axis == 'x':
        return np.array([[1, 0, 0], [0, cos(theta), -sin(theta)],
                         [0, sin(theta), cos(theta)]])
    elif axis == 'y':
        return np.array([[cos(theta), 0, sin(theta)], [0, 1, 0],
                         [-sin(theta), 0, cos(theta)]])
    elif axis == 'z':
        return np.array([[cos(theta), -sin(theta), 0],
                         [sin(theta), cos(theta), 0], [0, 0, 1]])

# input: a rotation matrix, output: theta x, y, z in deg
# https://learnopencv.com/rotation-matrix-to-euler-angles/
def rotationMatrixToEulerAngles(R):
    # assert(isRotationMatrix(R))
    sy = sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
    singular = sy < 1e-6
    if not singular:
        x = atan2(R[2, 1], R[2, 2])
        y = atan2(-R[2, 0], sy)
        z = atan2(R[1, 0], R[0, 0])
    else:
        x = atan2(-R[1, 2], R[1, 1])
        y = atan2(-R[2, 0], sy)
        z = 0
    return np.rad2deg(np.array([x, y, z]))


# ti_i-1
def get_ti2i_1(i, theta=None):
    # init_printing(use_unicode=True)  # use pretty math output
    # np.set_printoptions(precision=2, suppress=True)

    # fill in dh tbl wrt robot arms' dh params

    # array idx starts frm 0, so i-1
    alp, ai, di, t = symbols('alp, ai, di, t')
    (alp, ai, di) = dh_tbl[i - 1, :]
    if theta is None:
        t = 'q' + str(i)
    else:
        t = theta
    #ci = Symbol('cos'+str(i))
    #si = Symbol('sin'+str(i))
    """
    m = np.array([[cos(t), -sin(t), 0, ai],
                  [
                      sin(t) * round(cos(alp), 2),
                      cos(t) * round(cos(alp), 2),
                      round(-sin(alp),2),
                      round(-sin(alp),2) * di
                  ],
                  [
                      sin(t) * round(sin(alp), 2),
                      cos(t) * round(sin(alp), 2),
                      round(cos(alp),2),
                      round(cos(alp),2) * di
                  ], [0, 0, 0, 1]])

    """
    m = np.array(
        [[cos(t), -sin(t), 0, ai],
         [sin(t) * cos(alp),
          cos(t) * cos(alp), -sin(alp), -sin(alp) * di],
         [sin(t) * sin(alp),
          cos(t) * sin(alp),
          cos(alp),
          cos(alp) * di], [0, 0, 0, 1]])

    if theta is None:
        #t=t.evalf(2)
        #t = np.round(t.astype(np.double), 2)
        print(f't{i}-{i-1}:', m)
        return m
    else:
        print(f't{i}-{i-1}:', m)
        #print (f't{i}-{i-1}:', np.round(t.astype(np.double),2))
        #return (np.format_float_scientific(m))
        return m.astype(float)


'''
ntu:
acos(x)+bsin(x)=c
u=tan(x/2)
cosx=1-u**2/1+u**2
sinx=2*u/1+u**2
subs cosx and sinx and multiple 1+u**2 on both sides
a(1-u**2)+2bu=c(1+u**2) -> (a+c)u^2-2bu+(c-a)=0
一元2次方程式: u=b+/-srqt(b^2+a^2-c^2)/a+c
x=2atan2(b+/-srqt(b^2+a^2-c^2), a+c)
'''


def trig_equ(a, b, c):
    np.set_printoptions(precision=3, suppress=True)
    r = np.sqrt(a**2 + b**2)
    alp = atan2(b, a)
    #r*cos(q+alp)=c
    # or 360-
    qNalp1 = acos(c / r)
    qNalp2 = 2 * pi - qNalp1
    q_1 = (qNalp1 - alp).real
    q_2 = (qNalp2 - alp).real
    print('q3:', q_1 * 180 / pi, q_2 * 180 / pi)
    return (q_1, q_2)


def is_negative_number_digit(n: str) -> bool:
    try:
        int(n)
        return True
    except ValueError:
        return False


def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def extract_num(inp_str):
    # inp_str=inp_str.strip()
    # inp_str= inp_str.replace("+ ", "+")
    # inp_str= inp_str.replace("- ", "-")
    newG = ''
    num = 0

    print("Original String : ", inp_str)
    s = '+'
    for c in inp_str.split():
        if c == '+' or c == '-':
            s = c  # keep the sign of nxt numeric
            newG += c
        #elif c.isnumeric():
        elif isfloat(c):
            num = num - float(c) if s == '-' else num + float(c)
            newG = newG[:-1]
            print("Extracted numbers from the list : ", num)
        else:
            newG += c
    # conver string to sympy expressions
    t = (trigsimp(newG), num)
    return t


def fk_3axes(l1, l2, l3, q1, q2, q3):
    x = l1 * cos(q1) + l2 * cos(q1 + q2) + l3 * cos(q1 + q2 + q3)
    y = l1 * sin(q1) + l2 * sin(q1 + q2) + l3 * sin(q1 + q2 + q3)
    print('x, y:', x, y)
    return (x, y)


def fk_6axes(q):
    return get_ti2i_1(1, q[0]) @ get_ti2i_1(2, q[1]) @ get_ti2i_1(
        3, q[2]) @ get_ti2i_1(4, q[3]) @ get_ti2i_1(5, q[4]) @ get_ti2i_1(
            6, q[5])


#l3=206
# endEffector # tuple (x,y)
def verify_ik(pc_0, l1, l2, l3, q1s, q2s, q3s):
    endEffector = pc_0
    for q1 in q1s:
        for q2 in q2s:
            for q3 in q3s:
                if (fk_3axes(l1, l2, l3, q1, q2, q3) == endEffector):
                    print(q2, q3)
