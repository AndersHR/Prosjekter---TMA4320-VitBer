from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers_static_sun import RK4_exp_step_static_sun
from Gravitational_force import der_r, der_v_SI
path.__delitem__(0)

from numpy import arange, empty, zeros, array, sqrt, arctan, cos, sin, arctan2, geomspace, savez, flip
from scipy.constants import G, pi
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd


EARTH_RADIUS = 6.3781e6

def get_angle(x, y):
    if x<0:
        deg = arctan(y/x)+pi
    elif y<0:
        deg = arctan(y/x)+2*pi
    else:
        deg = arctan(y/x)
    return deg

def get_r(x, y):
    return sqrt(x**2+y**2)

def RK4_exp_static_earth_transfer(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v, delta_v):
    sum_delta_v = 0
    min_t = 5000
    delta_v = delta_v*h
    R_max = 35680e3+EARTH_RADIUS

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

    check = True

    for i in range(1, t.shape[0]):
        print (t[i])
        v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                         x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
        if (sqrt(x[0, i]**2+y[0, i]**2)<EARTH_RADIUS):
            print ("CRASH")
            break
        if (check):
            if (t[i]>min_t and sqrt(x[0, i]**2+y[0, i]**2)<R_max):
                print ("Accelerating")
                angle = arctan2(v_y[0, i], v_x[0, i])
                sum_delta_v += delta_v
                v_x[0, i] += delta_v * cos(angle)
                v_y[0, i] += delta_v * sin(angle)
            elif (sqrt(x[0, i]**2+y[0, i]**2)>=R_max):
                new_v = sqrt(G*m[-1]/sqrt(x[0, i]**2+y[0, i]**2))
                angle = arctan2(y[0, i], x[0, i])
                new_vx = -new_v*sin(angle)
                new_vy = new_v*cos(angle)
                sum_delta_v += sqrt((v_x[0, i]-new_vx)**2+(v_y[0, i]-new_vy)**2)
                v_x[0, i], v_y[0, i] = new_vx, new_vy
                check=False

    return t, x, y, v_x, v_y, sum_delta_v

def RK4_exp_static_earth_transfer_v2(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v, delta_v_in):
    sum_delta_v = 0
    min_t = 5000
    delta_v = 0
    R_max = 35680e3+EARTH_RADIUS

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

    check = True

    for i in range(1, t.shape[0]):
        print (t[i])
        v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_earth_transfer(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                         x[:, i - 1], y[:, i - 1], m, h, der_r, der_v, delta_v)
        if (check):
            if (t[i]>min_t and sqrt(x[0, i]**2+y[0, i]**2)<R_max):
                print ("Accelerating")
                delta_v = delta_v_in
                sum_delta_v += delta_v*h
            elif (sqrt(x[0, i]**2+y[0, i]**2)>=R_max):
                new_v = sqrt(G*m[-1]/sqrt(x[0, i]**2+y[0, i]**2))
                angle = arctan2(y[0, i], x[0, i])
                new_vx = -new_v*sin(angle)
                new_vy = new_v*cos(angle)
                sum_delta_v += sqrt((v_x[0, i]-new_vx)**2+(v_y[0, i]-new_vy)**2)
                v_x[0, i], v_y[0, i] = new_vx, new_vy
                check=False
                delta_v = 0
        if (sqrt(x[0, i]**2+y[0, i]**2)<EARTH_RADIUS):
            print ("CRASH")
            break

    return t, x, y, v_x, v_y, sum_delta_v

def RK4_exp_step_earth_transfer(n, v_x, v_y, x, y, m, h, der_r, der_v, acc):
    s1_x, s1_y, s1_vx, s1_vy = empty(n), empty(n), zeros(n), zeros(n)
    s2_x, s2_y, s2_vx, s2_vy = empty(n), empty(n), zeros(n), zeros(n)
    s3_x, s3_y, s3_vx, s3_vy = empty(n), empty(n), zeros(n), zeros(n)
    s4_x, s4_y, s4_vx, s4_vy = empty(n), empty(n), zeros(n), zeros(n)

    angle = arctan2(v_y, v_x)
    s1_vx += der_v(-x, -y) * m[-1] + acc*cos(angle)
    s1_vy += der_v(-y, -x) * m[-1] + acc*sin(angle)
    for i in range(n):
        s1_x[i] = der_r(v_x[i])
        s1_y[i] = der_r(v_y[i])
        for j in range(i+1, n):
            a = der_v(x[j] - x[i], y[j] - y[i])
            s1_vx[i] += a * m[j]
            s1_vx[j] -= a * m[i]

            a = der_v(y[j] - y[i], x[j] - x[i])
            s1_vy[i] += a * m[j]
            s1_vy[j] -= a * m[i]

    angle = arctan2(v_y+(h/2)*s1_vy, v_x+(h/2)*s1_vx)
    s2_vx += der_v(-(x + (h / 2) * s1_x), -(y + (h / 2) * s1_y)) * m[-1] + acc*cos(angle)
    s2_vy += der_v(-(y + (h / 2) * s1_y), -(x + (h / 2) * s1_x)) * m[-1] + acc*sin(angle)
    for i in range(n):
        s2_x[i] = der_r(v_x[i] + (h / 2) * s1_vx[i])
        s2_y[i] = der_r(v_y[i] + (h / 2) * s1_vy[i])
        for j in range(i+1, n):
            a = der_v(x[j] + (h / 2) * s1_x[j] - (x[i] + (h / 2) * s1_x[i]),
                              y[j] + (h / 2) * s1_y[j] - (y[i] + (h / 2) * s1_y[i]))
            s2_vx[i] += a * m[j]
            s2_vx[j] -= a * m[i]

            a = der_v(y[j] + (h / 2) * s1_y[j] - (y[i] + (h / 2) * s1_y[i]),
                              x[j] + (h / 2) * s1_x[j] - (x[i] + (h / 2) * s1_x[i]))
            s2_vy[i] += a * m[j]
            s2_vy[j] -= a * m[i]

    #------------------------------
    angle = arctan2(v_y+(h/2)*s2_vy, v_x+(h/2)*s2_vx)
    s3_vx += der_v(-(x + (h / 2) * s2_x), -(y + (h / 2) * s2_y)) * m[-1] + acc*cos(angle)
    s3_vy += der_v(-(y + (h / 2) * s2_y), -(x + (h / 2) * s2_x)) * m[-1] + acc*sin(angle)
    for i in range(n):
        s3_x[i] = der_r(v_x[i] + (h / 2) * s2_vx[i])
        s3_y[i] = der_r(v_y[i] + (h / 2) * s2_vy[i])

        for j in range(i+1, n):
            a = der_v(x[j] + (h / 2) * s2_x[j] - (x[i] + (h / 2) * s2_x[i]),
                              y[j] + (h / 2) * s2_y[j] - (y[i] + (h / 2) * s2_y[i]))
            s3_vx[i] += a * m[j]
            s3_vx[j] -= a * m[i]

            a = der_v(y[j] + (h / 2) * s2_y[j] - (y[i] + (h / 2) * s2_y[i]),
                              x[j] + (h / 2) * s2_x[j] - (x[i] + (h / 2) * s2_x[i]))
            s3_vy[i] += a * m[j]
            s3_vy[j] -= a * m[i]

    angle = arctan2(v_y+h*s3_vy, v_x+h*s3_vx)
    s4_vx += der_v(-(x + h * s3_x), -(y + h * s3_y)) * m[-1] + acc*cos(angle)
    s4_vy += der_v(-(y + h * s3_y), -(x + h * s3_x)) * m[-1] + acc*sin(angle)
    for i in range(n):
        s4_x[i] = der_r(v_x[i] + h * s3_vx[i])
        s4_y[i] = der_r(v_y[i] + h * s3_vy[i])
        for j in range(i+1, n):
            a = der_v(x[j] + h * s3_x[j] - (x[i] + h * s3_x[i]),
                              y[j] + h * s3_y[j] - (y[i] + h * s3_y[i]))
            s4_vx[i] += a * m[j]
            s4_vx[j] -= a * m[i]

            a = der_v(y[j] + h * s3_y[j] - (y[i] + h * s3_y[i]),
                              x[j] + h * s3_x[j] - (x[i] + h * s3_x[i]))
            s4_vx[i] += a*m[j]
            s4_vy[j] -= a*m[i]

    v_x_new, v_y_new, x_new, y_new = empty(n), empty(n), empty(n), empty(n)

    for i in range(n):
        v_x_new[i] = v_x[i] + (h / 6) * (s1_vx[i] + 2 * s2_vx[i] + 2 * s3_vx[i] + s4_vx[i])
        v_y_new[i] = v_y[i] + (h / 6) * (s1_vy[i] + 2 * s2_vy[i] + 2 * s3_vy[i] + s4_vy[i])
        x_new[i]   = x[i] + (h / 6) * (s1_x[i] + 2 * s2_x[i] + 2 * s3_x[i] + s4_x[i])
        y_new[i]   = y[i] + (h / 6) * (s1_y[i] + 2 * s2_y[i] + 2 * s3_y[i] + s4_y[i])
    return v_x_new, v_y_new, x_new, y_new

if __name__=="__main__":
    t_0 = 0
    h = 100

    m = array([720, 5.972*(10**24)])
    r = array([[EARTH_RADIUS+322e3, 0]])
    v = array([[0, sqrt(G*m[-1]/r[0, 0])]])

    delta_vs = geomspace(1e-4, 1, 20)
    t_ends = flip(array([1e5, 1.5e5, 1.5e5, 1.5e5, 2e5, 2.5e5, 2.5e5, 3.5e5, 5e5, 7.5e5,
    1e6, 1.65e6, 3e6, 4.5e6, 7e6, 1.2e7, 1.9e7,  3e7, 5e7,  5e7]), 0)

    sum_delta_vs = empty(delta_vs.shape[0])
    for i in range(delta_vs.shape[0]-1, -1, -1):
        if (delta_vs[i] < 2e-3):
            h = 100
        else:
            h = 10
        t, x, y, v_x, v_y, sum_delta_v = RK4_exp_static_earth_transfer_v2(1, t_0, t_ends[i], r, v, m, h, der_r, der_v_SI, delta_vs[i])

        sum_delta_vs[i] = sum_delta_v

        circle = plt.Circle((0, 0), EARTH_RADIUS, color='b')

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

    file = abs_path + r'/Data_files/Oppgave_3_2_cont.npz'
    savez(file, delta_vs=delta_vs, t_ends=t_ends, sum_delta_vs=sum_delta_vs)