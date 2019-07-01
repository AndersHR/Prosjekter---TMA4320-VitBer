from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))


from numpy import load, log, copy, exp, geomspace, empty, log10
from matplotlib import pyplot as plt
from scipy.stats import linregress
import matplotlib.gridspec as grd


TRUE_ALPHA = 1.1*1e-8

if __name__=="__main__":

    ALPHAS_TO_OMIT = 3

    file = abs_path + r'/Data_files/Oppgave_2_2.npz'
    npzfile = load(file)
    vals, alpha, arcsec = copy(npzfile['vals']), copy(npzfile['alpha']), copy(npzfile['arcsec_per_year'])
    npzfile.close()

    #Create logarithimatcilly distributed x-values from the enndpoints of the alphas
    x_vals = geomspace(alpha[-1], alpha[0], num=100)

    lin_reg = empty((alpha.shape[0]-ALPHAS_TO_OMIT, 2)) #[[intercept, slope].......]

    #Linear regression of the logarithms, gradually decreasing upper alpha-value used in the regression
    for i in range(alpha.shape[0]-ALPHAS_TO_OMIT):
        lin_reg[i, 1], lin_reg[i, 0], r_value, p_value, std_err = linregress(log10(alpha[i:]), log10(arcsec[i:]))


    fig = plt.figure('LogLog_plot')
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)
    gs = grd.GridSpec(1, 1)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel(r'$\alpha$   \   $\mathrm{Au}^2$', fontsize=28)
    ax.set_ylabel(r'$\mathrm{Precession}$   \   $\frac{Arcseconds}{\mathrm{Century}}$', fontsize=28)
    ax.grid(color="lightgrey", linestyle='dashed')
    ax.loglog(alpha, arcsec*100, 'b-', linewidth=4, label='Num. results')
    ax.loglog(x_vals, 10**(lin_reg[0, 0]+lin_reg[0, 1]*log10(x_vals))*100, linestyle = '--', color='r', linewidth=2,
              label= r'Lin.reg: $\alpha\in[10^{-2}, 3*10^{-6}]$')
    ax.loglog(x_vals, 10**(lin_reg[-1, 0]+lin_reg[-1, 1]*log10(x_vals))*100, linestyle = '--', color='c', linewidth=2,
              label=r'Lin.reg: $\alpha\in[4.2*10^{-6}, 3*10^{-6}]$')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
               borderaxespad=0., fontsize=22,
               bbox_transform=ax.transAxes)
    plt.show()


    fig1 = plt.figure('Figur E')
    gs = grd.GridSpec(1, 1)
    plt.rc('xtick', labelsize=22)
    plt.rc('ytick', labelsize=22)

    ax1 = fig1.add_subplot(gs[0, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.set_xlabel(r'$\alpha$   \   $\mathrm{Au}^2$', fontsize=28)
    ax1.set_ylabel(r'$\mathrm{Precession}$   \   $\frac{Arcseconds}{\mathrm{Century}}$', fontsize=28)
    ax1.axhline(43, color='grey')
    ax1.semilogx(alpha[:alpha.shape[0]-ALPHAS_TO_OMIT], 10**(lin_reg[:, 0]+lin_reg[:, 1]*log10(TRUE_ALPHA))*100, 'r-')
    ax1.invert_xaxis()
    plt.show()

    print (10**((log10(0.43)-lin_reg[-1, 0])/(lin_reg[-1, 1])))
    print (10**(lin_reg[-1, 0]+log10(3.66106496936e-09)*lin_reg[-1, 1]))


