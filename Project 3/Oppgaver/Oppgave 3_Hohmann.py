from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_step_static_sun
from Gravitational_force import der_r, der_v_SI
path.__delitem__(0)

from numpy import arange, empty, array, sqrt, savez, cos, sin, linspace, arctan2
from scipy.constants import G, pi
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd


EARTH_RADIUS = 6.3781e6
R_0 = 322e3 + EARTH_RADIUS
R_max = 35680e3 + EARTH_RADIUS

def get_angle(y, x):
    deg = arctan2(y, x)
    if (deg <0):
        deg += 2*pi
    return deg

def RK4_exp_static_earth_transfer(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v, r_transfer, t_first_transfer):
    sum_delta_v = 0

    t = arange(t_0, t_end + h / 2, h)
    v_x = empty((n_objects, t.shape[0]))
    v_y = empty((n_objects, t.shape[0]))
    x = empty((n_objects, t.shape[0]))
    y = empty((n_objects, t.shape[0]))
    for i in range(n_objects):
        v_x[i, 0] = v_0[i, 0]
        v_y[i, 0] = v_0[i, 1]
        x[i, 0] = r_0[i, 0]
        y[i, 0] = r_0[i, 1]
    i = 1
    while (i<t.shape[0]):
        print (i*h)
        if (t[i] < t_first_transfer):
            v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                             x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
            i+=1
        else:
            print ("SHIT", i)
            break
    for i_r in range(r_transfer.shape[0]):
        r_initial = sqrt(x[0, i-1]**2+y[0, i-1]**2)

        delta_v = sqrt(G*m[-1]/r_initial)*(sqrt(2*r_transfer[i_r]/(r_initial + r_transfer[i_r]))-1)
        sum_delta_v += delta_v

        angle = get_angle(y[0, i-1], x[0, i-1])
        v_x [0, i-1] -= delta_v*sin(angle)
        v_y [0, i-1] += delta_v*cos(angle)
        angle += pi
        angle_round = 0

        new_angle = get_angle(y[0, i-1], x[0, i-1])
        while (new_angle<angle):
            print(i * h)
            v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                                        x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
            new_angle = get_angle(y[0, i], x[0, i])+angle_round*2*pi
            if (new_angle<(angle-pi)):
                angle_round += 1
                new_angle += 2*pi
            i += 1
            if (i==t.shape[0]):
                break

        time_end = t[i]

        delta_v = sqrt(G*m[-1]/r_transfer[i_r]) * (1-sqrt(2*r_initial/(r_initial+r_transfer[i_r])))
        sum_delta_v += delta_v
        angle = arctan2(v_y[0, i - 1], v_x[0, i - 1])
        v_x[0, i - 1] = v_x[0, i - 1] + delta_v * cos(angle)
        v_y[0, i - 1] = v_y[0, i - 1] + delta_v * sin(angle)
    while (i < t.shape[0]):
        print (i*h)
        v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                         x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
        i += 1


    return t, x, y, v_x, v_y, sum_delta_v, time_end


if __name__=="__main__":
    n_transfers = 10

    t_0 = 0
    t_first_transfer = 15000
    h = 1


    m =array([720, 5.972*(10**24)])
    r = array([[EARTH_RADIUS+322e3, 0]])
    v = array([[0, sqrt(G*m[-1]/r[0, 0])]])
    t_ends = array([1.5e5, 1.5e5, 2e5, 2e5, 2e5, 2.5e5, 2.5e5, 3e5, 3e5, 4e5])

    sum_delta_vs = empty(n_transfers)

    for i in range(n_transfers):
        r_transfer = linspace(R_0, R_max, i+2)[1:]


        t, x, y, v_x, v_y, sum_delta_v, time_end = RK4_exp_static_earth_transfer(1, t_0, t_ends[i], r, v, m, h, der_r, der_v_SI, r_transfer, t_first_transfer)
        sum_delta_vs[i] = sum_delta_v

        #print (time_end-t_first_transfer)
        #print (2*pi*sqrt(((r[0, 0]+R_max)/2)**3/(G*(m[-1])))/2)


        if (True):
            fig1 = plt.figure('Figur r')
            plt.rc('xtick', labelsize=22)
            plt.rc('ytick', labelsize=22)
            gs = grd.GridSpec(1, 1)

            ax = fig1.add_subplot(gs[0, 0])
            ax.grid(color="lightgrey", linestyle='dashed')
            ax.set_xlabel(r'$x$   \   km', fontsize=28)
            ax.set_ylabel(r'$y$   \   km', fontsize=28)
            plt.plot(0, 0, 'bo', label='Earth', ms=28)
            circle = plt.Circle((0, 0), EARTH_RADIUS * 1e-3, color='b', label='Earth')
            ax.add_artist(circle)
            line, = ax.plot(x[0] * 1e-3, y[0] * 1e-3, 'k-', linewidth=1.5, label='Satelite')
            ax.set_aspect('equal')
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand",
                      borderaxespad=0., fontsize=22,
                      bbox_transform=ax.transAxes, markerscale=0.7)
            ax.set_axisbelow(True)
            #plt.show()


    file = abs_path + r'/Data_files/Oppgave_3_2_Hohmann.npz'
    savez(file, n_transfers=array([n_transfers]), sum_delta_vs=sum_delta_vs)

