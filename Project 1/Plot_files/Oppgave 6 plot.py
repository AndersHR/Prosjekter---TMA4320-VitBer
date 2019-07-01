from numpy import load, pi, sin, exp, linspace
from sys import path
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from Plotting import init_plot, plot_normal, show_plot, save_fig, plot_logy_test, save_fig_legend_on_top, set_axis_specs, plot_small
path.__delitem__(0)

def rho(x):
    w, y = 3*pi, -2
    return sin(w*x)*exp(y*x)


if __name__ == "__main__":
    npzfile_d0025 = load(rel_path+r'/Data_files/Oppgave6_d0025.npz')
    npzfile_d025 = load(rel_path+r'/Data_files/Oppgave6_d025.npz')
    npzfile_d25 = load(rel_path+r'/Data_files/Oppgave6_d25.npz')

    fig1 = plt.figure(r'b and \tilde{b}')
    plt.rc('xtick', labelsize=17)
    plt.rc('ytick', labelsize=17)
    gs = grd.GridSpec(90, 1)

    ax1 = plt.subplot(gs[:30, 0])
    ax1.set_ylabel(r'F(x)', fontsize=24)
    ax1.set_xlim(0, 1)
    ax1.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        labelbottom='off')  # labels along the bottom edge are off
    ax1.yaxis.set_label_coords(-0.05, -0.5)
    ax1.grid(color="lightgrey")
    ax1.plot(npzfile_d0025['x_s'], npzfile_d0025['b_tilde'], 'b-', linewidth=3, label=r'$\tilde{b}$')
    ax1.plot(npzfile_d0025['x_s'], npzfile_d0025['b'], 'y--', linewidth=3, label=r'$b$')
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol = 2, mode="expand", borderaxespad=0., fontsize=26,
              bbox_transform=ax1.transAxes)
    ax1.text(0.99, 0.95, r'd=0.025', fontsize=24, horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))

    ax2 = plt.subplot(gs[31:60, 0])
    ax2.grid(color="lightgrey")
    ax2.set_xlim(0, 1)
    ax2.plot(npzfile_d025['x_s'], npzfile_d025['b_tilde'], 'b-', linewidth=3)
    ax2.plot(npzfile_d025['x_s'], npzfile_d025['b'], 'y--', linewidth=3)
    ax2.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        labelbottom='off')  # labels along the bottom edge are off
    ax2.text(0.99, 0.95, r'd=0.25', fontsize=24, horizontalalignment='right', verticalalignment='top',
             transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))

    ax3 = plt.subplot((gs[61:, 0]))
    ax3.grid(color="lightgrey")
    ax3.set_xlim(0, 1)
    ax3.plot(npzfile_d25['x_s'], npzfile_d25['b_tilde'], 'b-', linewidth=3)
    ax3.plot(npzfile_d25['x_s'], npzfile_d25['b'], 'y--', linewidth=3)
    ax3.set_xlabel(r'x', fontsize=24)
    ax3.text(0.99, 0.95, r'd=2.5', fontsize=24, horizontalalignment='right', verticalalignment='top',
             transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))


    plt.show()

    #------------------------------------------------------------------------------------------------------
    fig2 = plt.figure('rho - d=0.025')
    plt.rc('xtick', labelsize=17)
    plt.rc('ytick', labelsize=17)
    gs = grd.GridSpec(1, 60)

    ax1 = plt.subplot(gs[0, :20])
    ax1.set_ylabel(r'$\rho$(x)', fontsize=24)
    ax1.set_xlabel(r'x', fontsize=24)
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.7, 0.7)
    ax1.set_yticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax1.grid(color="lightgrey")
    ax1.plot(npzfile_d0025['x_s'], rho(npzfile_d0025['x_s']), 'r-', linewidth=5, label=r'Analytic $\rho$')
    ax1.plot(npzfile_d0025['x_s'], npzfile_d0025['rho_hat'], 'y--', linewidth=3, label=r'$\hat{\rho}$')
    ax1.plot(npzfile_d0025['x_s'], npzfile_d0025['rho_hat_tilde'], 'ko--', linewidth=3, label=r'$\tilde{\rho}$')
    ax1.legend(bbox_to_anchor=(0.12, 0.91, .8, .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=26,
               bbox_transform=fig2.transFigure)
    ax1.text(0.5, 0.97, r'd=0.025', fontsize=24, horizontalalignment='center', verticalalignment='top',
             transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))


    ax2 = plt.subplot(gs[0, 21:40])
    ax2.set_xlabel(r'x', fontsize=24)
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-0.7*1e12, 0.7*1e12)
    ax2.set_yticks([-0.6*1e12, -0.3*1e12, 0, 0.3*1e12, 0.6*1e12])
    ax2.grid(color="lightgrey")
    ax2.tick_params(
        axis='y',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        labelleft='off')  # labels along the bottom edge are off
    ax2.plot(npzfile_d025['x_s'], rho(npzfile_d025['x_s']), 'r-', linewidth=5, label=r'Analytic $\rho$')
    ax2.plot(npzfile_d025['x_s'], npzfile_d025['rho_hat'], 'y--', linewidth=3, label=r'$\hat{\rho}$')
    ax2.plot(npzfile_d025['x_s'], npzfile_d025['rho_hat_tilde'], 'ko--', linewidth=3, label=r'$\tilde{\rho}$')
    ax2.text(0.5, 0.97, r'd=0.25', fontsize=24, horizontalalignment='center', verticalalignment='top',
             transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))

    ax3 = plt.subplot(gs[0, 41:])
    ax3.set_xlabel(r'x', fontsize=24)
    ax3.set_xlim(-0.05, 1.05)
    ax3.set_ylim(-0.7*1e14, 0.7*1e14)
    ax3.set_yticks([-0.6*1e14, -0.3*1e14, 0, 0.3*1e14, 0.6*1e14])
    ax3.grid(color="lightgrey")
    ax3.tick_params(
        axis='y',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        labelleft='off')  # labels along the bottom edge are off
    ax3.plot(npzfile_d25['x_s'], rho(npzfile_d25['x_s']), 'r-', linewidth=5, label=r'Analytic $\rho$')
    ax3.plot(npzfile_d25['x_s'], npzfile_d25['rho_hat'], 'y--', linewidth=3, label=r'$\hat{\rho}$')
    ax3.plot(npzfile_d25['x_s'], npzfile_d25['rho_hat_tilde'], 'ko--', linewidth=3, label=r'$\tilde{\rho}$')
    ax3.text(0.5, 0.97, r'd=2.5', fontsize=24, horizontalalignment='center', verticalalignment='top',
             transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5))
    plt.show()





