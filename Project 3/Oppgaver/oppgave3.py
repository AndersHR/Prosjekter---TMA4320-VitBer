from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, der_v_SI
path.__delitem__(0)

from numpy import array,sqrt,savez
import matplotlib.pyplot as plt

M_s = 1.989*(10**30)            # kg / M_s
AU = 1.496*(10**8)              # m / AU
G = 6.67408*(10**-11)           #
YEAR = 3.154*(10**7)            # s / year

M_j = 5.972*(10**24)

if __name__=='__main__':
    m_satellite = 720           # kg
    r_start = 328371            # m
    r_geostat = 35686371        # m

    t_0 = 0
    t_end1 = 59.2#4
    r = array([[r_start,0]])
    v = array([[0,sqrt((G*M_j)/r[0,0])]])
    h = 0.01


    #Èn periode ved r = 322 km: T = 57.5

    t_1, x_1, y_1, v_x1, v_y1 = RK4_exp_static_sun(1, t_0, t_end1, r, v, array([m_satellite,M_j]), h, der_r, der_v_SI)

    t_end2 = 12028.5#8
    delta_v1 = 14201.51

    print(v_x1[0,-1])
    print(v_y1[0,-1])

    t_2, x_2, y_2, v_x2, v_y2 = RK4_exp_static_sun(1, t_1[-1], t_1[-1]+t_end2,  array([[x_1[-1,0],y_1[-1,0]]]),
                                                                    array([[v_x1[-1,0],v_y1[-1,0]+delta_v1]]),
                                                                    array([m_satellite, M_j]), h, der_r, der_v_SI)

    t_end3 = 90000
    delta_v2 = 2889.80

    print(v_x2[0, -1])
    print(v_y2[0, -1])

    t_3, x_3, y_3, v_x3, v_y3 = RK4_exp_static_sun(1, t_2[-1], t_2[-1]+t_end3,  array([[x_2[0, -1], y_2[0, -1]]]),
                                                                    array([[v_x2[0, -1], v_y2[0, -1] - delta_v2]]),
                                                                    array([m_satellite, M_j]), h, der_r, der_v_SI)

    K_1 = (1/2)*m_satellite*(v_x1[0]**2 + v_y1[0]**2)
    V_1 = - G*m_satellite*M_j/sqrt(x_1[0]**2 + y_1[0]**2)
    E_1 = K_1 + V_1

    K_2 = (1 / 2) * m_satellite * (v_x2[0] ** 2 + v_y2[0] ** 2)
    V_2 = - G * m_satellite * M_j / sqrt(x_2[0] ** 2 + y_2[0] ** 2)
    E_2 = K_2 + V_2

    K_3 = (1 / 2) * m_satellite * (v_x3[0] ** 2 + v_y3[0] ** 2)
    V_3 = - G * m_satellite * M_j / sqrt(x_3[0] ** 2 + y_3[0] ** 2)
    E_3 = K_3 + V_3


    ax = plt.subplot(111)

    ax.plot(0,0,'ob',Markersize=2)
    ax.plot(x_1[0], y_1[0], '-k', Linewidth=1)
    ax.plot(x_2[0], y_2[0], '-r', Linewidth=1)
    ax.plot(x_3[0], y_3[0], '-b', Linewidth=1)
    plt.show()

    file = abs_path + r'/Data_files/Oppgave_3.npz'
    savez(file, t_1=t_1, x_1=x_1, y_1=y_1, v_x1=v_x1, v_y1=v_y1, K_1=K_1, V_1=V_1, E_1=E_1,
                t_2=t_2, x_2=x_2, y_2=y_2, v_x2=v_x2, v_y2=v_y2, K_2=K_2, V_2=V_2, E_2=E_2,
                t_3=t_3, x_3=x_3, y_3=y_3, v_x3=v_x3, v_y3=v_y3, K_3=K_3, V_3=V_3, E_3=E_3,
                delta_v1=delta_v1,delta_v2=delta_v2)
