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
    t_0, t_end = 0, 500
    h = 1e-3
    N_OBJECTS = 9

    m = empty(N_OBJECTS)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = array([0, 0]), array([1.7e-3, -9e-4]), 1
    #r[0], v[0], m[0] = array([0, 0]), array([0, 0]), 1
    r[1], v[1], m[1] = get_start_vals('Mercury', 0, 0)
    r[2], v[2], m[2] = get_start_vals('Venus', -1, 0)
    r[3], v[3], m[3] = get_start_vals('Earth', 1, 0)
    r[4], v[4], m[4] = get_start_vals('Mars', 0, 0)
    r[5], v[5], m[5] = get_start_vals('Jupiter', 0.5, 0)
    r[6], v[6], m[6] = get_start_vals('Saturn', -0.5, 0)
    r[7], v[7], m[7] = get_start_vals('Uranus', 1, 0)
    r[8], v[8], m[8] = get_start_vals('Neptune', -1, 0)
    r[4, 1], v[4, 0] = r[4, 1]*-1, v[4, 0]*-1
    r[6, 1], v[6, 0] = r[6, 1]*-1, v[6, 0]*-1



    #m[3] = 317.89/333480
    #print (m[3])

    t, x, y, v_x, v_y = RK4_exp(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v)#, der_v_der_primary, der_v_der_secondary)

    print (max(v_x[0]))
    print ("Non-static sun")



    plt.figure('1')
    ax = plt.subplot(111)
    ax.set_xlim(-30, 30)
    ax.set_ylim(-30, 30)
    ax.set_aspect('equal')

    ax.plot(x[0], y[0], 'yo', ms=10)
    ax.plot(x[1], y[1], 'k-')
    ax.plot(x[2], y[2], 'g-')
    ax.plot(x[3], y[3], 'b-')
    ax.plot(x[4], y[4], 'r-')
    ax.plot(x[5], y[5], 'c-')
    ax.plot(x[6], y[6], 'm-')
    ax.plot(x[7], y[7], 'k-')
    ax.plot(x[8], y[8], 'g-')
    plt.show()