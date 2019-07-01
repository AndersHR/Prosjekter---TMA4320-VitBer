from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Protein_class')
from Protein_class import Protein
path.__delitem__(0)

from numpy import linspace, exp, empty, arange
from matplotlib import pyplot as plt
from matplotlib import animation
from numpy.random import seed

def d(T, d_max, s):
    return d_max*exp(-s*T)

def init():
    ax.set_ylim(90, 250)
    ax.set_xlim(80, 220)
    ax.set_aspect('equal')
    ax.set_title('Animation of protein', fontsize=17)
    ax.grid(color='lightgrey')
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text,

def animate(i):
    line.set_data(x[i], y[i])
    time_text.set_text('{}'.format(i))
    return line, time_text

if __name__ == "__main__":
    seed(5476)
    protein = Protein(300, 100)
    it = 0
    x, y = empty((100, 300)), empty((100, 300))
    while (it < 71):
        print (it)
        if protein.bend_protein():
            x[it], y[it] = protein.get_pos_as_lists()
            it += 1

    # Set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    line, = ax.plot([], [], 'ro-', lw=2)
    time_text = ax.text(0.03, 0.05, '', transform=ax.transAxes, fontsize=14,
                        bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5})

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=arange(30, 72, 1), interval=1000, blit=True)

    plt.show()


