from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals
path.__delitem__(0)

from numpy import array, arange
from scipy.constants import pi
import matplotlib.pyplot as plt
from matplotlib import animation

"""
g = 4*pi**2

def der_r(v):
    return v

def der_v(primary_coord_diference, secondary_cood_diff):
    return (g*(primary_coord_diference)/((primary_coord_diference)**2+(secondary_cood_diff)**2)**(3/2))

def der_v_der_primary(primary_coord_diference, secondary_cood_diff):
    return g*((1-3*(primary_coord_diference)**2)/((primary_coord_diference)**2+(secondary_cood_diff)**2)) / \
           (((primary_coord_diference)**2+(secondary_cood_diff)**2)**(3/2))

def der_v_der_secondary(primary_coord_diference, secondary_cood_diff):
    return -g*(3*(primary_coord_diference*secondary_cood_diff))/(((primary_coord_diference)**2+(secondary_cood_diff)**2)**(5/2))
"""
def init():
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.set_title('Animation of TBP', fontsize=17)
    ax.grid(color='lightgrey')
    planet.set_data([], [])
    planet2.set_data([], [])
    sun.set_data([], [])
    time_text.set_text('')
    return planet, planet2, sun, time_text,


def animate(i):
    planet.set_data(x[1, 4*i], y[1, 4*i])
    planet2.set_data(x[2, 4*i], y[2, 4*i])
    sun.set_data(x[0, 4*i], y[0, 4*i])
    print (x[0, 4*i], y[0, 4*i])
    time_text.set_text('{:.1f}'.format(i*1e-3*4))
    return planet, planet2, sun, time_text


if __name__ == "__main__":
    t_0, t_end = 0, 10
    r = array([[0, 0], [0.46668776, 0], [0.72821844, 0], ])
    v = array([[0, 0], [0, 8.197591649], [0, 7.337825638]])
    m = [1, 0.0553/333480, 0.8150/333480]
    m_static_sun = [0.0553/333480, 0.8150/333480, 1]
    r_static_sun = array([[0.46668776, 0], [0.72821844, 0], ])
    v_static_sun = array([[0, 8.197591649], [0, 7.337825638]])
    h = 1e-3

    t, x, y, v_x, v_y = EC1_exp(3, t_0, t_end, r, v, m, h, der_r, der_v)#, der_v_der_primary, der_v_der_secondary)
    # Set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(x[1], y[1], 'b-', linewidth=1)
    ax.plot(x[2], y[2], 'g-', linewidth=1)
    planet, = ax.plot([], [], 'bo', markersize=12)
    planet2, = ax.plot([], [], 'go', markersize=12)
    sun, = ax.plot([0], [0], 'ro', markersize=20)
    time_text = ax.text(0.03, 0.05, '', transform=ax.transAxes, fontsize=14,
                        bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})
    print ("TEST")
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=int(t_end/(4*h)), interval=10, blit=True)

    plt.show()