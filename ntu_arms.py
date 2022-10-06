
import pieper as pp
import craig as cg
import numpy as np
from sympy import Matrix, Symbol, Eq, solve, sin, cos

def ntu_pieper():
    #warnings.filterwarnings('ignore')
    #from inverse_kinematic import trig_equ
    dh_tbl=np.array([[0,0,0],
                        [np.deg2rad(-90), -30, 0],
                        [0, 340, 0],
                        [np.deg2rad(-90), -40, 338],
                        [np.deg2rad(90), 0, 0],
                        [np.deg2rad(-90),0,0]])

    cg.setDhTbl(dh_tbl)
    #getcontext().prec=2
    #getcontext().rounding = ROUND_UP

    ty = np.deg2rad(-60) # rotate y axis
    tcup_0_2s = np.array([[cos(ty), 0, sin(ty), 330], [0, 1, 0, 372],
                            [-sin(ty), 0, cos(ty), 367], [0, 0, 0, 1]])
    Tcup_6 = np.array([[0, 0, 1, 0], [0, -1, 0, 0], [1, 0, 0, 206],
                        [0, 0, 0, 1]])
    tcup_0_2s = np.array([[cos(ty), 0, sin(ty), 330], [0, 1, 0, 372],
                            [-sin(ty), 0, cos(ty), 367], [0, 0, 0, 1]])

    t6_0 = tcup_0_2s @ np.linalg.inv(Tcup_6)
    pp.pieper(t6_0)

ntu_pieper()

