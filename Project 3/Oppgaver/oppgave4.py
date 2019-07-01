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

if __name__=="__main__":
    N_OBJECTS = 1

    t_0, t_end = 0, 50
    h = 1e-3


    m = empty(3)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = get_start_vals('Earth', 1, 0)
    #r[1], v[1], m[1] = get_start_vals('Mars', -1, 0)
    #r[2, 0], r[2, 1] = 0, 0
    #v[2, 0], v[2, 1] = 0, 0
    m[-1] = 1

    #m[0] = 1e-3

    print (m)
    print (r)
    print (v)

    t, x, y, v_x, v_y = LG4_imp_static_sun(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v, der_v_der_primary, der_v_der_secondary)


    plt.figure('r')
    #plt.plot(x[2], x[2], 'y-', ms=20)
    plt.plot(x[0], y[0], 'b-')
    #plt.plot(x[1], y[1], 'r-')
    plt.show()