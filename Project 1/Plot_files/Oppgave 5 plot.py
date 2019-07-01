from numpy import load
from sys import path
from os.path import dirname
rel_path = dirname(dirname(__file__))

path.insert(0, rel_path+r'/Help_files')
from Plotting import init_plot, plot_normal, show_plot, save_fig, plot_logy_test, save_fig_legend_on_top, set_axis_specs, plot_small
path.__delitem__(0)


file = rel_path + r'/Data_files/Oppgave5.npz'
if __name__ == "__main__":

    npz_file = load(file)
    fig, ax = init_plot('Question 5 - Error')
    ax = set_axis_specs(ax, r'$N_c$', r'Max $|$Error$_j|$', 36, 36)

    N_c_list, d_0025, d_025, d_25 = npz_file['N_c_list'], npz_file['max_diff_1'], npz_file['max_diff_2'], npz_file['max_diff_3']

    plot_logy_test(ax, N_c_list, d_0025, 'ro--', r'Err(N$_{c}$) - d = ' + '{}'.format(0.025), linewidth=4, ms=10)
    plot_logy_test(ax, N_c_list, d_025, 'bo--', r'Err($N_{c}$) - d = ' + '{}'.format(0.25),  linewidth=4, ms=10)
    plot_logy_test(ax, N_c_list, d_25, 'go--', r'Err($N_c$) - d ={}'.format(2.5),  linewidth=4, ms=10)

    ax.set_xlim(4.8, 30.2)

    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=28,
              bbox_transform=ax.transAxes)



