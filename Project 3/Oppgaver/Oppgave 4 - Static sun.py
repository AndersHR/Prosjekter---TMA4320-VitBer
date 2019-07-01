from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals
path.__delitem__(0)

from numpy import array, pi, empty, abs, sqrt, savez
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd



G = 4*pi**2

def squared(x, y):
    return sqrt(x**2+y**2)

if __name__=="__main__":
    t_0, t_end = 0, 500
    h = 1e-3
    N_OBJECTS = 2

    m = empty(N_OBJECTS+1)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = get_start_vals('Earth', 1, 0)
    r[1], v[1], m[1] = get_start_vals('Mars', -1, 0)

    m[-1] = 1
        #317.89/333480

    t, x, y, v_x, v_y = RK4_exp_static_sun(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v)#, der_v_der_primary, der_v_der_secondary)

    fig1 = plt.figure('Figur r')
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
    ax.plot(0, 0, 'yo', ms=40, label='Sun')
    ax.plot(x[0], y[0], 'b-', linewidth=1.5, label='Earth')
    ax.plot(x[1], y[1], 'r-', linewidth=1.5, label='Mars')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax.transAxes, markerscale=0.7)

    plt.show()

    r_p = squared(x[0]-x[1], y[0]-y[1])
    Rs = empty((2, t.shape[0]))
    for i in range(N_OBJECTS):
        v = squared(v_x[i], v_y[i])
        r_s = squared(x[i], y[i])
        V = - G*333480*(m[i]*m[-1]/r_s+ m[0]*m[1]/r_p)
        K = (1 / 2) * m[i] * (v ** 2) * 333480
        E = V + K
        Rs[i] = r_s

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
        ax1.plot(t, v_x[i], 'b-', linewidth=1.5, label='$V_x(t)$')
        ax1.plot(t, v_y[i], 'r-', linewidth=1.5, label='$V_y(t)$')
        ax1.plot(t, v, 'k-', linewidth=1.5, label='$V(t)$')
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                   borderaxespad=0., fontsize=22,
                   bbox_transform=ax1.transAxes, markerscale=0.7)

        ax2 = fig1.add_subplot(gs[6:, 0])
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
    for i in range(1, t.shape[0]):
        Rs[0, i] = abs(Rs[0, i] - Rs[0, 0])
        Rs[1, i] = abs(Rs[1, i] - Rs[1, 0])
    Rs[0, 0] = abs(Rs[0, 0] - Rs[0, 0])
    Rs[1, 0] = abs(Rs[1, 0] - Rs[1, 0])

    fig1 = plt.figure('Figur V')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(141, 1)

    ax1 = fig1.add_subplot(gs[:70, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_ylabel(r'$\|r_s-r_0\|$    \   Au', fontsize=28)
    ax1.tick_params(
        axis='x',
        which='both',
        labelbottom='off'
    )
    mars_E, = ax1.semilogy(t, Rs[1], 'r-', linewidth=1.5, label='$Mars$')

    ax2 = fig1.add_subplot(gs[71:, 0])
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.set_xlabel(r'$t$   \   years', fontsize=28)
    ax2.set_ylabel(r'$\|r_s-r_0\|$    \   Au', fontsize=28)
    earth_E, = ax2.semilogy(t, Rs[0], 'b-', linewidth=1.5, label='$Earth$')
    ax2.legend(handles = [mars_E, earth_E], bbox_to_anchor=(0.1, 1.02, .9, .102), loc=3, ncol=3, mode="expand",
               borderaxespad=0., fontsize=22,
               bbox_transform=ax1.transAxes, markerscale=0.7)

    plt.show()
