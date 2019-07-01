from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals
path.__delitem__(0)

from numpy import array, pi, empty, abs, sqrt, savez
from matplotlib import pyplot as plt

if __name__== "__main__":
    t_0 = 0
    t_end = 100
    h = 1e-3


    r, v, m = empty((1, 2)), empty((1, 2)), empty(2)
    r[0], v[0], m[0] = get_start_vals('Mercury', -1, pi)
    m[-1] = 1
    t, x, y, v_x, v_y =RK4_exp_static_sun(1, t_0, t_end, r, v, m, h, der_r, der_v)

    plt.figure('Test')
    ax = plt.subplot(111)
    ax.plot(0, 0, 'yo', ms=20)
    ax.plot(x[0], y[0], 'b-', linewidth=0.5)
    ax.set_aspect('equal')
    plt.show()
