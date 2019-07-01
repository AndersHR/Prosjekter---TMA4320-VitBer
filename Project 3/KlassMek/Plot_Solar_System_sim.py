from numpy import load, arange, empty
from matplotlib.pyplot import figure, subplot2grid, plot, text, xlim, ylim, xlabel,ylabel, rc, legend, show, subplot, title
from matplotlib import animation

if __name__ == "__main__":
    solar_system_coords = load('Solar_System_coords.npz')

    x_coords = solar_system_coords['x']
    y_coords = solar_system_coords['y']

    num_planets = len(x_coords)

    fig = figure(dpi=200)
    #fig = figure()

    ax = subplot2grid((6, 7), (0, 0), colspan=6, rowspan=7)
    title('Numerical Simulation of the Solar System')
    #ax = subplot(111)
    ax.plot(0, 0, 'oy', MarkerSize=5)
    colors = ['red', 'yellow', 'cyan', 'black', 'red', 'orange', 'green', 'blue', 'purple', 'gray', 'indigo', '.r',
              '.r']
    labels = ['Mercury', 'Venus', 'Earth', 'Eros (asteroid)', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto',
              '''Halley's comet''', 'test', 'test']

    time_box = dict(boxstyle='square,pad=0.3', fc='white')
    #ax = text(22, -45, ('Time lapsed = %d' % end) +
    #              ' years', fontsize=6, bbox=time_box)

    limit = 40
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel('$x$ (AU)', fontsize=8)
    ax.set_ylabel('$y$ (AU)', fontsize=8)
    #ax = rc('xtick', labelsize=8)
    #ax = rc('ytick', labelsize=8)

    planet_size = 2

    planet_1, = ax.plot([], [], 'o', color=colors[0], label=labels[0], Markersize=planet_size)
    planet_2, = ax.plot([], [], 'o', color=colors[1], label=labels[1], Markersize=planet_size)
    planet_3, = ax.plot([], [], 'o', color=colors[2], label=labels[2], Markersize=planet_size)
    planet_4, = ax.plot([], [], 'o', color=colors[3], label=labels[3], Markersize=planet_size)
    planet_5, = ax.plot([], [], 'o', color=colors[4], label=labels[4], Markersize=planet_size)
    planet_6, = ax.plot([], [], 'o', color=colors[5], label=labels[5], Markersize=planet_size)
    planet_7, = ax.plot([], [], 'o', color=colors[6], label=labels[6], Markersize=planet_size)
    planet_8, = ax.plot([], [], 'o', color=colors[7], label=labels[7], Markersize=planet_size)
    planet_9, = ax.plot([], [], 'o', color=colors[8], label=labels[8], Markersize=planet_size)
    planet_10, = ax.plot([], [], 'o', color=colors[9], label=labels[9], Markersize=planet_size)
    planet_11, = ax.plot([], [], 'o', color=colors[10], label=labels[10], Markersize=planet_size)

    tail_width = 0.5

    tail_1, = ax.plot([], [], '-', color=colors[0], linewidth=tail_width)
    tail_2, = ax.plot([], [], '-', color=colors[1], linewidth=tail_width)
    tail_3, = ax.plot([], [], '-', color=colors[2], linewidth=tail_width)
    tail_4, = ax.plot([], [], '-', color=colors[3], linewidth=tail_width)
    tail_5, = ax.plot([], [], '-', color=colors[4], linewidth=tail_width)
    tail_6, = ax.plot([], [], '-', color=colors[5], linewidth=tail_width)
    tail_7, = ax.plot([], [], '-', color=colors[6], linewidth=tail_width)
    tail_8, = ax.plot([], [], '-', color=colors[7], linewidth=tail_width)
    tail_9, = ax.plot([], [], '-', color=colors[8], linewidth=tail_width)
    tail_10, = ax.plot([], [], '-', color=colors[9], linewidth=tail_width)
    tail_11, = ax.plot([], [], '-', color=colors[10], linewidth=tail_width)

    if False:
        for planet in range(num_planets):
            x = x_coords[planet]
            y = y_coords[planet]

            ax.plot(x, y, linestyle='-', color=colors[planet])
            ax.plot(x[-1], y[-1], linestyle='-', color=colors[planet], linewidth=0.3)


    def init():
        """
        planet_obj = empty(num_planets)
        tail_obj = empty(num_planets)
        for j in range(num_planets):
            planet_obj[j], = ax.plot([],[],'o',color=colors[j],label=labels[j],Markersize=3)
            tail_obj[j], = ax.plot([], [], '-', color=colors[j], label=labels[j],linewidth=1)
        return planet_obj,tail_obj

        """
        planet_1, = ax.plot([], [], 'o', color=colors[0], label=labels[0], Markersize=planet_size)
        planet_2, = ax.plot([], [], 'o', color=colors[1], label=labels[1], Markersize=planet_size)
        planet_3, = ax.plot([], [], 'o', color=colors[2], label=labels[2], Markersize=planet_size)
        planet_4, = ax.plot([], [], 'o', color=colors[3], label=labels[3], Markersize=planet_size)
        planet_5, = ax.plot([], [], 'o', color=colors[4], label=labels[4], Markersize=planet_size)
        planet_6, = ax.plot([], [], 'o', color=colors[5], label=labels[5], Markersize=planet_size)
        planet_7, = ax.plot([], [], 'o', color=colors[6], label=labels[6], Markersize=planet_size)
        planet_8, = ax.plot([], [], 'o', color=colors[7], label=labels[7], Markersize=planet_size)
        planet_9, = ax.plot([], [], 'o', color=colors[8], label=labels[8], Markersize=planet_size)
        planet_10, = ax.plot([], [], 'o', color=colors[9], label=labels[9], Markersize=planet_size)
        planet_11, = ax.plot([], [], 'o', color=colors[10], label=labels[10], Markersize=planet_size)

        """
        tail_1, = ax.plot([], [], '-', color=colors[0], linewidth=tail_width)
        tail_2, = ax.plot([], [], '-', color=colors[1], linewidth=tail_width)
        tail_3, = ax.plot([], [], '-', color=colors[2], linewidth=tail_width)
        tail_4, = ax.plot([], [], '-', color=colors[3],  linewidth=tail_width)
        tail_5, = ax.plot([], [], '-', color=colors[4],  linewidth=tail_width)
        tail_6, = ax.plot([], [], '-', color=colors[5],  linewidth=tail_width)
        tail_7, = ax.plot([], [], '-', color=colors[6],  linewidth=tail_width)
        tail_8, = ax.plot([], [], '-', color=colors[7],  linewidth=tail_width)
        tail_9, = ax.plot([], [], '-', color=colors[8],  linewidth=tail_width)
        tail_10, = ax.plot([], [], '-', color=colors[9],  linewidth=tail_width)
        tail_11, = ax.plot([], [], '-', color=colors[10],  linewidth=tail_width)
        """

        return planet_1,planet_2,planet_3,planet_4,planet_5,planet_6,planet_7,planet_8,planet_9,planet_10,planet_11#,tail_1,tail_2,tail_3,tail_4,tail_5,tail_6,tail_7,tail_8,tail_9,tail_10,tail_11



    def animate(i):
        """
        for j in range(num_planets):
            planet_obj[j].set_data(x_coords[j][i],y_coords[j][i])
            tail_obj[j].set_data(x_coords[j][0:i],y_coords[j][0:i])
        return planet_obj[0],planet_obj[1],planet_obj[2],planet_obj[3],planet_obj[4],planet_obj[5],planet_obj[6],tail_obj[0],tail_obj[5]
        """
        planet_1.set_data(x_coords[0][i], y_coords[0][i])
        planet_2.set_data(x_coords[1][i], y_coords[1][i])
        planet_3.set_data(x_coords[2][i], y_coords[2][i])
        planet_4.set_data(x_coords[3][i], y_coords[3][i])
        planet_5.set_data(x_coords[4][i], y_coords[4][i])
        planet_6.set_data(x_coords[5][i], y_coords[5][i])
        planet_7.set_data(x_coords[6][i], y_coords[6][i])
        planet_8.set_data(x_coords[7][i], y_coords[7][i])
        planet_9.set_data(x_coords[8][i], y_coords[8][i])
        planet_10.set_data(x_coords[9][i], y_coords[9][i])
        planet_11.set_data(x_coords[10][i], y_coords[10][i])

        max_length = 5000000
        if i < max_length:
            tail_length = i
        else:
            tail_length = max_length
        """
        tail_1.set_data(x_coords[0][i - tail_length:i], y_coords[0][i - tail_length:i])
        tail_2.set_data(x_coords[1][i - tail_length:i], y_coords[1][i - tail_length:i])
        tail_3.set_data(x_coords[2][i - tail_length:i], y_coords[2][i - tail_length:i])
        tail_4.set_data(x_coords[3][i - tail_length:i], y_coords[3][i - tail_length:i])
        tail_5.set_data(x_coords[4][i - tail_length:i], y_coords[4][i - tail_length:i])
        tail_6.set_data(x_coords[5][i - tail_length:i], y_coords[5][i - tail_length:i])
        tail_7.set_data(x_coords[6][i - tail_length:i], y_coords[6][i - tail_length:i])
        tail_8.set_data(x_coords[7][i - tail_length:i], y_coords[7][i - tail_length:i])
        tail_9.set_data(x_coords[8][i - tail_length:i], y_coords[8][i - tail_length:i])
        tail_10.set_data(x_coords[9][i - tail_length:i], y_coords[9][i - tail_length:i])
        tail_11.set_data(x_coords[10][i - tail_length:i], y_coords[10][i - tail_length:i])
        """
        return planet_1, planet_2, planet_3, planet_4, planet_5, planet_6, planet_7, planet_8, planet_9, planet_10, planet_11#, tail_1, tail_2, tail_3, tail_4, tail_5, tail_6, tail_7, tail_8, tail_9, tail_10, tail_11


    ax.legend(bbox_to_anchor=(1.26, 1), fontsize=6)

    anim = animation.FuncAnimation(fig,animate,arange(1000,x_coords[0].size,100),init_func=init,interval=2)

    show()

    if True:

        limit = 50

        fig = figure(dpi=200)

        ax = subplot2grid((6, 7), (0, 0), colspan=6, rowspan=7)
        title('Numerical Simulation of the Solar System')

        ax.plot(0, 0, 'oy', MarkerSize=5)

        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xlabel('$x$ (AU)', fontsize=8)
        ax.set_ylabel('$y$ (AU)', fontsize=8)

        for planet in range(num_planets):
            x = x_coords[planet]
            y = y_coords[planet]


            ax.plot(x, y, linestyle='-', color=colors[planet], label=labels[planet])
            ax.plot(x[-1], y[-1], linestyle='-', color=colors[planet], linewidth=0.5)
            ax.legend(bbox_to_anchor=(1.26, 1), fontsize=6)

        show()