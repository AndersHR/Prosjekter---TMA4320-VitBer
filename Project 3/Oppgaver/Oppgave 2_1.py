from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_static_sun
from Gravitational_force import der_r, der_v, get_start_vals

path.__delitem__(0)

import matplotlib.gridspec as grd
from numpy import pi, empty, sqrt
from matplotlib import pyplot as plt

G = 4*pi**2

def squared(x, y):
    return sqrt(x**2+y**2)

if __name__== "__main__":
    t_0 = 0
    t_end = 2*.2408
    h = 1e-3


    r, v, m = empty((1, 2)), empty((1, 2)), empty(2)
    r[0], v[0], m[0] = get_start_vals('Mercury', -1, 0)
    r[0, 0], v[0, 1] = r[0, 0]*-1, v[0, 1]*-1
    print (r[0], v[0])
    m[-1] = 1
    t, x, y, v_x, v_y =RK4_exp_static_sun(1, t_0, t_end, r, v, m, h, der_r, der_v)

    fig1 = plt.figure('Figur r')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)

    ax = fig1.add_subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.set_xlabel(r'$x$   \   Au', fontsize=28)
    ax.set_ylabel(r'$y$   \   Au', fontsize=28)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(-0.5, 0.5)
    ax.plot(0, 0, 'yo', ms=40, label='Sun')
    ax.plot(x[0], y[0], 'b-', linewidth=1.5, label='Mercury')
    ax.set_aspect('equal')
    lgnd = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand",
              borderaxespad=0., fontsize=22,
              bbox_transform=ax.transAxes, markerscale=0.7)

    plt.show()


    v = squared(v_x[0], v_y[0])
    r = squared(x[0], y[0])

    V = - G * m[0] / r * 333480
    K = (1 / 2) * m[0] * (v**2) * 333480
    E = V + K


    fig1 = plt.figure('Figur V')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(11, 1)

    ax1 = fig1.add_subplot(gs[:5, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_ylabel(r'$v$   \   $\frac{\mathrm{Au}}{\mathrm{year}}$', fontsize=28)
    ax1.tick_params(
        axis='x',
        which='both',
        labelbottom='off'
    )
    ax1.plot(t, v_x[0], 'b-', linewidth=1.5, label='$V_x(t)$')
    ax1.plot(t, v_y[0], 'r-', linewidth=1.5, label='$V_y(t)$')
    ax1.plot(t, v, 'k-', linewidth=1.5, label='$V(t)$')
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax1.transAxes, markerscale=0.7)

    ax2 = fig1.add_subplot(gs[6:, 0])
    ax2.set_ylim(-8, 6)
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.set_xlabel(r'$t$   \   year', fontsize=28)
    ax2.set_ylabel(r'$E$   \   $\frac{\mathrm{m}_\mathrm{E}\mathrm{Au}^2}{\mathrm{year}^2}$', fontsize=28)
    ax2.plot(t, K, 'b-', linewidth=1.5, label='$K(t)$')
    ax2.plot(t, V, 'r-', linewidth=1.5, label='$V(t)$')
    ax2.plot(t, E, 'k-', linewidth=1.5, label='$E(t)$')
    ax2.legend(bbox_to_anchor=(0.1, 1.02, .9, .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax2.transAxes, markerscale=0.7)

    plt.show()
