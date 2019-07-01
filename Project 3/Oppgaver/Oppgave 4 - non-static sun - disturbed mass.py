from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals
path.__delitem__(0)

from numpy import array, empty
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd


if __name__=="__main__":


    t_0, t_end = 0, 500
    h = 1e-4
    N_OBJECTS = 5

    m = empty(N_OBJECTS)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = array([-0.00190482984382/2, 0]), array([-2.8e-7, -5.865e-3]), 1
    r[1], v[1], m[1] = get_start_vals('Mercury', 0, 0)
    r[2], v[2], m[2] = get_start_vals('Venus', -1, 0)
    r[3], v[3], m[3] = get_start_vals('Earth', 1, 0)
    r[4], v[4], m[4] = get_start_vals('Mars', 0, 0)
    r[4, 1], v[4, 0] = r[4, 1]*-1, v[4, 0]*-1

    m[3] = 317.89/333480
    print (m[3])

    t, x, y, v_x, v_y = RK4_exp(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v)#, der_v_der_primary, der_v_der_secondary)

    print (max(v_x[0]))
    print ("Non-static sun")

    new_shape = x[0, ::10].shape[0]
    x_new, y_new = empty((N_OBJECTS, new_shape)), empty((N_OBJECTS, new_shape))
    for i in range(N_OBJECTS):
        x_new[i] = x[i, ::10]
        y_new[i] = y[i, ::10]


    fig1 = plt.figure('Figur {}'.format(v[0, 1]))
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)
    ax = fig1.add_subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.set_xlabel(r'$x$   \   Au', fontsize=28)
    ax.set_ylabel(r'$y$   \   Au', fontsize=28)
    ax.set_ylim(-1.8, 1.8)
    ax.set_yticks(array([-1, 0, 1]))
    ax.set_xlim(-1.8, 1.8)
    ax.set_xticks(array([-1, 0, 1]))
    ax.plot(x_new[0], y_new[0], 'yo', ms=10, label = 'Sun')
    ax.plot(x_new[1], y_new[1], 'k-', linewidth=1.5, label='Mercury')
    ax.plot(x_new[2], y_new[2], 'g-', linewidth=1.5, label='Venus')
    ax.plot(x_new[3], y_new[3], 'b-', linewidth=1.5, label='Earth')
    ax.plot(x_new[4], y_new[4], 'r-', linewidth=1.5, label='Mars')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
              borderaxespad=0., fontsize=22,
              bbox_transform=ax.transAxes, markerscale=0.7)

    plt.show()