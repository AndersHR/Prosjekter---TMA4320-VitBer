from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))
path.insert(0, abs_path+r'/Computing_files')
from ODE_solvers import LG4_imp, RK4_exp, EC1_exp
from ODE_solvers_static_sun import RK4_exp_static_sun, EC1_exp_static_sun, LG4_imp_static_sun
from Gravitational_force import der_r, der_v, der_v_der_primary, der_v_der_secondary, get_start_vals, der_v_einstein
path.__delitem__(0)

from numpy import array, pi, empty, abs, sqrt, savez, append, arctan, geomspace, log, exp, polyfit
from matplotlib import pyplot as plt
from scipy.stats import linregress


RAD_TO_ARCSEC = 648000/pi
def get_der_v_einstein(alfa):
    return lambda x, y: der_v_einstein(x, y, alfa)

def squared(x, y):
    return x**2+y**2

def get_angle(x, y):
    deg = empty(x.shape[0])
    round = 0
    for i in range(x.shape[0]):
        if x[i]<0:
            deg[i] = arctan(y[i]/x[i]) + pi
        elif y[i]<0:
            deg[i] = arctan(y[i]/x[i]) + 2*pi
        else:
            deg[i] = arctan(y[i]/x[i])
    return deg

def get_angle_test (x, y):
    if x<0:
        deg = arctan(y/x)+pi
    elif y<0:
        deg = arctan(y/x)+2*pi
    else:
        deg = arctan(y/x)
    return deg

def plot(x, y, alpha):
    plt.figure('Alpha = {}'.format(alpha))
    plt.plot(x,  y, 'b-')
    plt.plot(0, 0, 'yo', ms=20)
    plt.show()

if __name__== "__main__":
    t_0 = 0
    t_end = 100
    h = 1e-4

    alpha = geomspace(1e-2, 3e-6, num=50) #5e-5
    #alpha = array([1e-4])
    print (alpha)
    arcsec_per_year =empty(alpha.shape[0])

    r, v, m = empty((1, 2)), empty((1, 2)), empty(2)
    r[0], v[0], m[0] = get_start_vals('Mercury', -1, 0)
    m[-1] = 1
    r[0, 0] = r[0, 0]*-1
    v[0, 1] = v[0, 1]*-1

    for iter in range(alpha.shape[0]):
        t, x, y, v_x, v_y = RK4_exp_static_sun(1, t_0, t_end, r, v, m, h, der_r, get_der_v_einstein(alpha[iter]))

        r0 = squared(x[0], y[0])
        deg = get_angle(x[0], y[0])

        #periphelion_degs = array([0])
        #periphelion_time = array([0])

        round = 0
        last_periphelion = 0

        if (alpha[iter]>=1e-4):
            check = False
            for i in range(1, t.shape[0]):
                if (r0[i] > r0[i - 1]):
                    if (check):
                        if (deg[i-1]<deg[last_periphelion]):
                            round += 1
                        #periphelion_degs =append(periphelion_degs, deg[i-1])
                        last_periphelion = i - 1
                        check = False
                elif (r0[i] < r0[i - 1]):
                    check = True
        else:
            check = False
            for i in range(1, t.shape[0]):
                if (r0[i] > r0[i - 1]):
                    if (check):
                        last_periphelion = i - 1
                        check = False
                elif (r0[i] < r0[i - 1]):
                    check = True

        #print (periphelion_degs)

        last_deg = deg[last_periphelion]
        last_periphelion -= 1

        """
        check = False
        for i in range(1, t.shape[0]):
            if (r0[i] > r0[i-1]):
                if(check):
                    last_periphelion = i-2
                    #periphelion_degs = append(periphelion_degs, deg[i-1])
                    #periphelion_time = append(periphelion_degs, t[i-1])
                    check = False
            elif (r0[i] < r0[i-1]):
                check = True
        """

        print (iter, 0)

        t_new, x_new, y_new, v_x_new, v_y_new = RK4_exp_static_sun(1, t[last_periphelion], t[last_periphelion+3],
                        array([[x[0, last_periphelion], y[0, last_periphelion]]]), array([[v_x[0, last_periphelion], v_y[0, last_periphelion]]]), m, h*1e-5, der_r, get_der_v_einstein(alpha[iter]))


        r0 = squared(x_new[0], y_new[0])
        deg = get_angle(x_new[0], y_new[0])

        check = False
        for i in range(1, t_new.shape[0]):
            if (r0[i] > r0[i - 1]):
                if (check):
                    if (alpha[iter]>=1e-4):
                        if (deg[i-1] > pi and last_deg < pi):
                            print ("Close call")
                            print (last_deg, deg[i-1])
                            round -=1
                        elif (deg[i-1] < pi and last_deg > pi):
                            print ("close Call")
                            print (last_deg, deg[i-1])
                            round += 1
                    arcsec_per_year[iter] = (deg[i-1]+round*2*pi)/t_new[i-1]
                    print (i)
                    check = False
            elif (r0[i] < r0[i - 1]):
                check = True

        print (round)
        #print(periphelion_degs)

        """
        for i in range(1, periphelion_degs.shape[0]):
            if (periphelion_degs[i] < periphelion_degs[i-1]):
                for j in range(i, periphelion_degs.shape[0]):
                    periphelion_degs[j] += 2*pi

        #print(periphelion_degs)
        """

        #arcsec_per_year[iter] = periphelion_degs[-1]/periphelion_time[-1]

        print (iter)

    arcsec_per_year = arcsec_per_year*RAD_TO_ARCSEC

    vals = array([t_0, t_end, h])

    #file = abs_path + r'/Data_files/Oppgave_2.2_v2.npz'
    #savez(file, vals=vals, alpha=alpha, arcsec_per_year=arcsec_per_year)

    print (arcsec_per_year)

    slope, intercept, r_value, p_value, std_err = linregress(log(alpha), log(arcsec_per_year))

    log_arcpy = intercept+slope*log(alpha)
    arc_py_1 =exp(log_arcpy)

    plt.figure('loglog')
    plt.loglog(alpha, arc_py_1)
    plt.loglog(alpha, arcsec_per_year, 'r-')
    plt.show()

    print (exp(intercept+slope*log(1.1*10**-8))*100)


