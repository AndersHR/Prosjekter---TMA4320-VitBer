from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary,get_start_vals,get_vals
path.__delitem__(0)

from numpy import array, empty, savez

planets = array(['Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'])

if __name__=='__main__':
    tau_list = empty(len(planets))
    a_list = empty(len(planets))
    for p in range(len(planets)):
        params = get_vals(planets[p])
        a_list[p] = params[0]
        r,v,m = get_start_vals(planets[p],1,0)
        t_start = 0
        t_end = 1.2*params[0]**(3/2)
        h = 10**-4

        t,x,y,vx,vy = RK4_exp_static_sun(1,t_start,t_end,array([r]),array([v]),[m,1],h,der_r,der_v)

        halfway = False
        for i in range(len(t)):
            if ((y[0,i] > 0) and (halfway)):
                tau = t[i]
                break
            elif (y[0,i] < 0) and not (halfway):
                halfway = True
        tau_list[p] = tau


    file = abs_path + r'/Data_files/Oppgave_1_keplers3rd.npz'
    savez(file,tau_list=tau_list,a_list=a_list)

    for j in range(len(tau_list)):
        c = (tau_list[j]**2)/(a_list[j]**3)
        print(planets[j],":",c)