from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, der_v_einstein, get_start_vals

path.__delitem__(0)

N = 10
TRAIL_LENGTH = 200

from numpy import array, arange, empty
from scipy.constants import pi
import matplotlib.pyplot as plt
from matplotlib import animation

def init():
    ax.set_ylim(-2, 2)
    ax.set_xlim(-2, 2)
    ax.set_aspect('equal')
    ax.set_title('Animation of TBP', fontsize=17)
    ax.grid(color='lightgrey')
    planet1.set_data([], [])
    trail1.set_data([], [])
    planet2.set_data([], [])
    trail2.set_data([], [])
    time_text.set_text('')
    return planet1, trail1, planet2, trail2, time_text,


def animate(i):
    if (i<TRAIL_LENGTH):
        trail1.set_data(x[0, :N*i+1:10], y[0, :N*i+1:10])
        trail2.set_data(x[1, :N*i+1:10], y[1, :N*i+1:10])
    else:
        trail1.set_data(x[0, N*(i-TRAIL_LENGTH):N * i + 1:10], y[0, N*(i-TRAIL_LENGTH):N * i + 1:10])
        trail2.set_data(x[1, N*(i-TRAIL_LENGTH):N * i + 1:10], y[1, N*(i-TRAIL_LENGTH):N * i + 1:10])
    planet1.set_data(x[0, N*i], y[0, N*i])
    planet2.set_data(x[1, N*i], y[1, N*i])
    time_text.set_text('{:.1f}'.format(i*1e-3/N))
    return planet1, trail1, planet2, trail2, time_text


def get_der_v_einstein(alfa):
    return lambda x, y: der_v_einstein(x, y, alfa)

if __name__ == "__main__":
    t_0, t_end = 0, 1
    h = 1e-4
    N_OBJECTS = 2

    m = empty(N_OBJECTS + 1)
    r = empty((N_OBJECTS, 2))
    v = empty((N_OBJECTS, 2))
    r[0], v[0], m[0] = get_start_vals('Earth', 1, 0)
    r[1], v[1], m[1] = get_start_vals('Mars', 0.1, 0)

    m[-1] = 1
    # m[1] = 317.89 / 333480
    # m[0] = 317.89 / 333480
    m[1] = m[-1]

    t, x, y, v_x, v_y = LG4_imp_static_sun(N_OBJECTS, t_0, t_end, r, v, m, h, der_r, der_v, der_v_der_primary,
                                           der_v_der_secondary)
    # Set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    planet1, = ax.plot([], [], 'bo', markersize=12)
    trail1,  = ax.plot([], [], 'b-', linewidth=1)
    planet2, = ax.plot([], [], 'ro', ms=12)
    trail2,  = ax.plot([], [], 'r-', linewidth=1)
    ax.plot([0], [0], 'ro', markersize=20)
    time_text = ax.text(0.03, 0.05, '', transform=ax.transAxes, fontsize=14,
                        bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=int(t_end/(N*h)), interval=10, blit=True)

    plt.show()