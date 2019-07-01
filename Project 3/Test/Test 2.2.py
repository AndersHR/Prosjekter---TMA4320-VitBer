from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))


from numpy import load, log, copy, exp, geomspace, empty
from matplotlib import pyplot as plt
from scipy.stats import linregress

TRUE_ALPHA = 1.1*1e-8

if __name__=="__main__":
    ALPHAS_TO_OMIT = 3

    file = abs_path + r'/Data_files/Oppgave_2_2.npz'
    npzfile = load(file)
    vals, alpha, arcsec = copy(npzfile['vals']), copy(npzfile['alpha']), copy(npzfile['arcsec_per_year'])
    npzfile.close()


    #Create logarithimatcilly distributed x-values from the enndpoints of the alphas
    x_vals = geomspace(alpha[0], alpha[-1], num=100)

    lin_reg = empty((alpha.shape[0]-ALPHAS_TO_OMIT, 2)) #[[intercept, slope].......]
    #Linear regression of the logarithms
    for i in range(alpha.shape[0]-ALPHAS_TO_OMIT):
        lin_reg[i, 1], lin_reg[i, 0], r_value, p_value, std_err = linregress(log(alpha[i:]), log(arcsec[i:]))


    plt.figure('LogLog_plot')
    ax = plt.subplot(111)
    ax.plot(alpha, arcsec, 'b-', linewidth=2)
    #for i in range(alpha.shape[0]-ALPHAS_TO_OMIT):
        #plt.loglog(x_vals, exp(lin_reg[i, 0]+lin_reg[i, 1]*log(x_vals)), 'g--')
    plt.show()

    for i in range(alpha.shape[0]-ALPHAS_TO_OMIT):
        print (exp(lin_reg[i, 0]+lin_reg[i, 1]*log(TRUE_ALPHA))*100)

    plt.figure('test')
    for i in range(alpha.shape[0]-ALPHAS_TO_OMIT):
        plt.plot(i, exp(lin_reg[i, 0]+lin_reg[i, 1]*log(TRUE_ALPHA))*100, 'ro')
    plt.show()


