from numpy import load
from matplotlib import pyplot as plt
from matplotlib import gridspec as grd

from sys import path
from os.path import dirname

abs_path = dirname(dirname(__file__))

if __name__=="__main__":
    file = abs_path + r'/Data_files/Oppgave_4_1.npz'
    npzfile = load(file)

    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)

    fig1 = plt.figure('Figur E')
    gs = grd.GridSpec(1, 1)

    ax1 = fig1.add_subplot(gs[0, 0])
    ax1.grid(color="lightgrey", linestyle='dashed')
    ax1.plot(npzfile['bend_list'], npzfile['E_list'], 'r-')
    plt.show()

    fig2 = plt.figure('Figur L')
    gs = grd.GridSpec(1, 1)

    ax2 = fig2.add_subplot(gs[0, 0])
    ax2.grid(color="lightgrey", linestyle='dashed')
    ax2.plot(npzfile['bend_list'], npzfile['L_list'], 'r-')
    plt.show(fig2)