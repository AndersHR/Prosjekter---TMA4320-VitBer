from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_step_static_sun
from Gravitational_force import der_r, der_v_SI
path.__delitem__(0)

from numpy import arange, empty, array, sqrt, arctan, cos, sin
from scipy.constants import G, pi
from matplotlib import pyplot as plt

def get_angle(x, y):
    if x<0:
        deg = arctan(y/x)+pi
    elif y<0:
        deg = arctan(y/x)+2*pi
    else:
        deg = arctan(y/x)
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
    while (t[i] < t_first_transfer):
        v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                         x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
        i += 1
    for i_r in range(r_transfer.shape[0]):
        r_initial = sqrt(x[0, i-1]**2+y[0, i-1]**2)

        delta_v = sqrt(G*m[-1]/r_initial)*(sqrt(2*r_transfer[i_r]/(r_initial + r_transfer[i_r]))-1)
        sum_delta_v += delta_v

        angle = get_angle(x[0, i-1], y[0, i-1])
        v_x [0, i-1] += delta_v*sin(angle)
        v_y [0, i-1] += delta_v*cos(angle)
        angle += pi
        print (angle)
        angle_round = 0

        new_angle = get_angle(x[0, i-1], y[0, i-1])
        while (new_angle<angle):
            print (new_angle)
            v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                                        x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
            new_angle = get_angle(x[0, i], y[0, i])+angle_round*2*pi
            if (new_angle<(angle-pi)):
                print ("TEST")
                angle_round += 1
                new_angle += 2*pi
            print (new_angle)
            i += 1
            if (i==t.shape[0]):
                break
        print (i)
        delta_v = sqrt(G*m[-1]/r_transfer[i_r]) * (1-sqrt(2*r_initial/(r_initial+r_transfer[i_r])))
        sum_delta_v += delta_v
        print (delta_v)
        angle = get_angle(v_x[0, i - 1], v_y[0, i - 1])
        v_x[0, i - 1] = v_x[0, i - 1] + delta_v * cos(angle)
        v_y[0, i - 1] = v_y[0, i - 1] + delta_v * sin(angle)
    while (i < t.shape[0]):
        v_x[:, i], v_y[:, i], x[:, i], y[:, i] = RK4_exp_step_static_sun(n_objects, v_x[:, i - 1], v_y[:, i - 1],
                                                                         x[:, i - 1], y[:, i - 1], m, h, der_r, der_v)
        i += 1


    return t, x, y, v_x, v_y, sum_delta_v

if __name__=="__main__":
    t_0, t_end = 0, 150000
    h = 1


    m =array([720, 5.972*(10**24)])
    r = array([[322e3, 0]])
    v = array([[0, sqrt(G*m[-1]/r[0, 0])]])
    #r_transfer = array([10000e3, 20000e3, 30000e3, 35680e3])
    r_transfer = array([35680e3])

    t, x, y, v_x, v_y, sum_delta_v = RK4_exp_static_earth_transfer(1, t_0, t_end, r, v, m, h, der_r, der_v_SI, r_transfer, 0.1)
    plt.figure('Test')
    plt.plot(x[0], y[0])
    plt.show()

    print (sum_delta_v)

