from os.path import dirname

abs_path = dirname(dirname(__file__))
from numpy import load, array, log10, pi, sqrt
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt

if __name__=="__main__":
    file = abs_path + r'/Data_files/Oppgave_1.npz'
    npz_file = load(file)

    delta_r_RK4 = npz_file['delta_r_RK4']
    delta_r_EC = npz_file['delta_r_EC']
    h_list = npz_file['h_list']

    ax = plt.subplot(111)

    ax.set_xlim(10**-0.5,10**-6.5)

    ax.set_xlabel(r'$h$ / years')
    ax.set_ylabel(r'$error$ / AU')

    ax.loglog(h_list[0:6],delta_r_RK4[0:6],'--or')
    ax.loglog(h_list[0:6],delta_r_EC[0:6],'--ob')

    plt.show()


    m = 1/333480

    t = npz_file['RK_t']

    x = (npz_file['RK_x'])[0]
    y = (npz_file['RK_y'])[0]
    vx = (npz_file['RK_vx'])[0]
    vy = (npz_file['RK_vy'])[0]

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    #plot av banen
    ax = plt.subplot(111)

    ax.plot(0, 0, 'oy', Markersize=12)
    ax.plot(x, y, '-r')

    ax.set_xlabel('$x$ / AU',fontsize=17)
    ax.set_ylabel('$y$ / AU',fontsize=17)

    plt.rc('xtick',labelsize=15)
    plt.rc('ytick',labelsize=15)

    plt.show()

    M_s = 1.989*(10**30)    # kg / M_s
    AU = 1.496*(10**8)      # m / AU
    G = 6.67408*(10**-11)   #
    YEAR = 3.154*(10**7)    # s / year


    V = - 4*(pi**2)*m*1/sqrt(x**2 + y**2)
    K = (1/2)*m*(vx**2 + vy**2)

    E = V + K

    #plot av potensiell, kinetisk og total energi
    ax = plt.subplot(111)

    ax.set_xlabel(r'$t$ / years',fontsize=15)
    ax.set_ylabel(r'$E$ / M$_{s}$ AU year$^{-1}$',fontsize=15)

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    ax.plot(t,E,'-b',label='Total energy, $E$')
    ax.plot(t,V,'-r',label='Potential energy, $V$')
    ax.plot(t,K,'-g',label='Kinetic energy, $K$')

    ax.legend(loc=1,bbox_to_anchor=(0., 0.85, 1., .102),fontsize=13)

    plt.show()

    #Plot av feil i E

    delta_E_RK4 = npz_file['delta_E_RK4']
    delta_E_EC = npz_file['delta_E_EC']

    ax = plt.subplot(111)

    ax.set_xlim(10 ** -0.5, 10 ** -6.5)

    ax.set_xlabel(r'$h$ / years',fontsize=20)
    ax.set_ylabel(r'$error$ $\frac{\Delta E}{E}$',fontsize=20)

    ax.loglog(h_list,delta_E_RK4 / abs(E[0]),'--or',label='Runge-Kutta 4th order')
    ax.loglog(h_list,delta_r_EC / abs(E[0]),'--ob',label='Euler-Cromer')

    ax.legend(loc=1,fontsize=13)

    plt.show()