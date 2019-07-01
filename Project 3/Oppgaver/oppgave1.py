from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary
path.__delitem__(0)

from numpy import array, pi, empty, abs, sqrt, savez

if __name__=="__main__":
    m = array([1/333480, 1]) #mass of sun must be the last parameter for the static_sun-solvers

    h_list = array([10**-1,10**-2,10**-3,10**-4,10**-5,10**-6])

    t_min = 0
    t_max = 1

    delta_r_RK4 = empty(len(h_list))
    delta_r_EC = empty(len(h_list))
    delta_r_LG = empty(len(h_list))

    delta_E_RK4 = empty(len(h_list))
    delta_E_EC = empty(len(h_list))
    delta_E_LG = empty(len(h_list))

    for i in range(len(h_list)):
        RK_t,RK_x,RK_y,RK_vx,RK_vy = RK4_exp_static_sun(1,0,1,array([[1,0]]),array([[0,-2*pi]]),[1/333490,1],h_list[i],der_r,der_v)
        delta_r_RK4[i] = sqrt((1-RK_x[0,-1])**2 + (RK_y[0,-1])**2)
        RK_E = (1/2)*m[0]*(RK_vx[0]**2 + RK_vy[0]**2) - m[0]*(4*(pi**2))/sqrt(RK_x[0]**2 + RK_y[0]**2)
        delta_E_RK4[i] = abs(RK_E[0] - RK_E[-1])

        EC_t,EC_x,EC_y,EC_vx,EC_vy = EC1_exp_static_sun(1, 0, 1, [1,0], [0,-2*pi], [1/333490,1], h_list[i], der_r, der_v)
        delta_r_EC[i] = sqrt((1-EC_x[0,-1]) ** 2 + (EC_y[0,-1] ** 2))
        EC_E = (1 / 2) * m[0] * (EC_vx[0] ** 2 + EC_vy[0] ** 2) - m[0] * (4 * (pi ** 2)) / sqrt(EC_x[0] ** 2 + EC_y[0] ** 2)
        delta_E_EC[i] = abs(EC_E[0] - EC_E[-1])

        LG_t, LG_x, LG_y, LG_vx, LG_vy = LG4_imp_static_sun(1, 0, 1, [1, 0], [0, -2 * pi], [1 / 333490, 1], h_list[i],
                                                            der_r, der_v, der_v_der_primary, der_v_der_secondary)
        delta_r_LG[i] = sqrt((1 - LG_x[0, -1]) ** 2 + (LG_y[0, -1] ** 2))
        LG_E = (1 / 2) * m[0] * (LG_vx[0] ** 2 + LG_vy[0] ** 2) - m[0] * (4 * (pi ** 2)) / sqrt(
            LG_x[0] ** 2 + LG_y[0] ** 2)
        delta_E_LG[i] = abs(LG_E[0] - LG_E[-1])

    file = abs_path + r'/Data_files/Oppgave_1.npz'
    savez(file, RK_t = RK_t, RK_x = RK_x, RK_y = RK_y, RK_vx = RK_vx, RK_vy = RK_vy, delta_r_RK4 = delta_r_RK4, RK_E=RK_E, delta_E_RK4=delta_E_RK4,
                EC_t = EC_t, EC_x = EC_x, EC_y = EC_y, EC_vx = EC_vx, EC_vy = EC_vy, delta_r_EC  = delta_r_EC,  EC_E=EC_E, delta_E_EC=delta_E_EC,
                LG_t = LG_t, LG_x = LG_x, LG_y = LG_y, LG_vx = LG_vx, LG_vy = LG_vy, delta_r_LG  = delta_r_LG,  LG_E=LG_E, delta_E_LG=delta_E_LG,
                h_list=h_list)

    print(delta_r_RK4)
    print()
    print(delta_r_EC)

