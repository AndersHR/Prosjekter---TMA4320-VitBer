from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_static_sun
from Gravitational_force import der_r, der_v_SI
path.__delitem__(0)

from numpy import empty, array, sqrt
from scipy.constants import G
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd

EARTH_RADIUS = 6.3781e6
R_0 = 322e3 + EARTH_RADIUS


def squared(x, y):
    return sqrt(x**2+y**2)

if __name__=="__main__":

    t_0, t_end = 0, 10000
    h = 1


    m =array([720, 5.972*(10**24)])
    r = array([[EARTH_RADIUS+322e3, 0]])
    v = array([[0, sqrt(G*m[-1]/r[0, 0])]])

    h_S = array([100, 10, 5, 1, 5e-1, 1e-1])
    ESS = empty(h_S.shape[0])
    for i in range (h_S.shape[0]):
        t, x, y, v_x, v_y = RK4_exp_static_sun(1, t_0, t_end, r, v, m, h_S[i], der_r, der_v_SI)
        R = squared(x[0], y[0])
        print (i)
        print (max(R), min(R))
        print (R[0])
        ESS[i] = (max(R)-min(R))

    print (ESS)
    plt.figure('TEST')
    plt.loglog(h_S, ESS, 'r-')
    plt.show()



    fig1 = plt.figure('Figur r')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)

    ax = fig1.add_subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.set_xlabel(r'$x$   \   km', fontsize=28)
    ax.set_ylabel(r'$y$   \   km', fontsize=28)
    plt.plot(0, 0, 'bo', label='Earth', ms=28)
    circle = plt.Circle((0, 0), EARTH_RADIUS*1e-3, color='b', label='Earth')
    ax.add_artist(circle)
    line, = ax.plot(x[0]*1e-3, y[0]*1e-3, 'k-', linewidth=1.5, label='Satelite')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand",
              borderaxespad=0., fontsize=22,
              bbox_transform=ax.transAxes, markerscale=0.7)
    ax.set_axisbelow(True)

    plt.show()

    v = squared(v_x[0], v_y[0])
    r = squared(x[0], y[0])

    V = - G *m[-1]* m[0] / r
    K = (1 / 2) * m[0] * (v**2)
    E = V + K

    t = t/3600

    fig1 = plt.figure('Figur V')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(11, 1)

    ax1 = fig1.add_subplot(gs[:5, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_ylabel(r'$v$   \   $\frac{\mathrm{km}}{\mathrm{s}}$', fontsize=28)
    ax1.tick_params(
        axis='x',
        which='both',
        labelbottom='off'
    )
    ax1.plot(t, v_x[0]*1e-3, 'b-', linewidth=1.5, label='$V_x(t)$')
    ax1.plot(t, v_y[0]*1e-3, 'r-', linewidth=1.5, label='$V_y(t)$')
    ax1.plot(t, v*1e-3, 'k-', linewidth=1.5, label='$V(t)$')
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax1.transAxes, markerscale=0.7)

    ax2 = fig1.add_subplot(gs[6:, 0])
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.set_xlabel(r'$t$   \   hour', fontsize=28)
    ax2.set_ylabel(r'$E$   \   J', fontsize=28)
    ax2.plot(t, K, 'b-', linewidth=1.5, label='$K(t)$')
    ax2.plot(t, V, 'r-', linewidth=1.5, label='$V(t)$')
    ax2.plot(t, E, 'k-', linewidth=1.5, label='$E(t)$')
    ax2.legend(bbox_to_anchor=(0.1, 1.02, .9, .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0., fontsize=22,
                     bbox_transform=ax2.transAxes, markerscale=0.7)

    plt.show()

