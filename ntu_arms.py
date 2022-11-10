
import pieper as pp
import craig as cg
import numpy as np
from sympy import sin, cos

def ntu_pieper():
    #warnings.filterwarnings('ignore')
    #from inverse_kinematic import trig_equ
    dh_tbl=np.array([[0,0,0],
                        [np.deg2rad(-90), -30, 0],
                        [0, 340, 0],
                        [np.deg2rad(-90), -40, 338],
                        [np.deg2rad(90), 0, 0],
                        [np.deg2rad(-90),0,0]], dtype=float)

    cg.setDhTbl(dh_tbl)

    ty = np.deg2rad(-60) # rotate y axis
    tcup_0_2s = np.array([[cos(ty), 0, sin(ty), 330], [0, 1, 0, 372],
                            [-sin(ty), 0, cos(ty), 367], [0, 0, 0, 1]], dtype=float)
    Tcup_6 = np.array([[0, 0, 1, 0], [0, -1, 0, 0], [1, 0, 0, 206],
                        [0, 0, 0, 1]], dtype=int)

    t6_0 = tcup_0_2s @ np.linalg.inv(Tcup_6)
    # t6_0=  np.round(t6_0.astype(np.double), 3)
    pp.pieper(t6_0)

ntu_pieper()
"""
q3: -11.98, 178.48
q2: -64.46, 20.37
q1: 58.61
q4: 25.30, -154.70
q5: -87.13, 87.13
q6: -56.19, 123.81
"""


