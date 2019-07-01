from numpy import empty, zeros, int8, amin, array, exp, sqrt, float64, copy, amax
from numpy.random import uniform, randint
from scipy.constants import k
from matplotlib import pyplot as plt
import matplotlib.gridspec as grd

class Protein:
    def __init__(self, length, T):
        self.grid = zeros((length + 2, length + 2), dtype=int8)
        for i in range(1, length + 1):
            self.grid[int((length + 1) / 2)][i] = i
        self.aa_pos = array([[int((length + 1) / 2), i + 1] for i in range(length)])  # [y, x]
        self.backup_grid = copy(self.grid)
        self.aa_pos_backup = copy(self.aa_pos)

        self.potentials = zeros((length, length), dtype=float64)
        for i in range(length - 2):
            self.potentials[i][i + 2:] = uniform(-10.4e-21, -3.47e-21, length - 2 - i)
        for i in range(2, length):
            for j in range(i - 1):
                self.potentials[i][j] = self.potentials[j][i]

        self.fixed_point = int((length + 1) / 2)
        self.length = length
        self.T = T

        self.u = 0.0
        self.L = length
        self.valid_bends = 0

        self.axis = 0
        self.direction = 1
        self.rotation = 0

    def choose_axis(self):
        self.axis = randint(2, self.length)
        if self.axis < (self.fixed_point):
            self.direction = -1
        else:
            self.direction = 1
        self.rotation = randint(0, 2, 1)[0]
        if self.rotation == 0:
            self.rotation = -1

    def check_for_collisions(self):
        if (self.direction == 1):
            end = self.length + 1
        else:
            end = 0
        for i in range(self.axis + self.direction, end, self.direction):
            new_pos = array([(self.aa_pos[i - 1, 1] - self.aa_pos[self.axis - 1, 1]) * self.rotation + self.aa_pos[
                self.axis - 1, 0],
                             (self.aa_pos[self.axis - 1, 0] - self.aa_pos[i - 1, 0]) * self.rotation + self.aa_pos[
                                 self.axis - 1, 1]])
            if (self.grid[new_pos[0], new_pos[1]] == 0 or self.grid[ new_pos[0], new_pos[1]] * self.direction > self.axis * self.direction):
                self.aa_pos[i - 1] = new_pos
            else:
                self.aa_pos = copy(self.aa_pos_backup)
                return False
        return True

    def bend(self):
        if (self.direction == 1):
            end = self.length + 1
        else:
            end = 0
        for i in range(self.axis + self.direction, end, self.direction):
            self.grid[self.aa_pos_backup[i - 1][0], self.aa_pos_backup[i - 1][1]] = 0
        for i in range(self.axis + self.direction, end, self.direction):
            self.grid[self.aa_pos[i - 1][0], self.aa_pos[i - 1][1]] = i

    def compute_E(self):
        new_E = 0
        for i in range(self.length):
            for j in range(-1, 2, 2):
                posy, posx = self.aa_pos[i][0], self.aa_pos[i][1]
                if (posy + j < self.length+2 and posy + j > -1):
                    if (self.grid[posy + j, posx] != 0):
                        new_E += self.potentials[i, self.grid[posy + j, posx] - 1]
                if (posx + j < self.length+2 and posx + j > -1):
                    if (self.grid[posy, posx + j] != 0):
                        new_E += self.potentials[i, self.grid[posy, posx + j] - 1]
        return new_E / 2

    def compute_L(self):
        L1 = 0
        for i in range(self.length - 1):
            for j in range(i + 1, self.length):
                cur_L = (self.aa_pos[i][0] - self.aa_pos[j][0]) ** 2 + (self.aa_pos[i][1] - self.aa_pos[j][1]) ** 2
                if cur_L > L1:
                    L1 = cur_L
        self.L = sqrt(L1)

    def check_for_fluctuations(self):
        new_E = self.compute_E()
        delta_u = new_E - self.u
        if (delta_u > 0):
            p_ratio = exp(-delta_u / (k * self.T))
            fluctuations = uniform()
            if (p_ratio < fluctuations):
                self.aa_pos = copy(self.aa_pos_backup)
                self.grid = copy(self.backup_grid)
                return False

        self.backup_grid = copy(self.grid)
        self.aa_pos_backup = copy(self.aa_pos)
        self.u = new_E
        return True

    def bend_protein(self):
        self.choose_axis()
        if (self.check_for_collisions()):
            self.bend()
            self.valid_bends += 1
            if (self.check_for_fluctuations()):
                self.compute_L()
                return True
        return False

    def get_pos_as_lists(self):
        x, y = empty(self.length), empty(self.length)

        for i in range(self.length):
            x[i] = (self.aa_pos[i][1])
            y[i] = (self.length + 1 - self.aa_pos[i][0])
        return x, y


    def plot_protein(self):

        Label = "length = {}, E = {:.2E}, L = {:.2E}".format(self.length, self.u, self.L)

        figure = plt.figure('Figur av protein')
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        gs = grd.GridSpec(1, 1)


        ax = plt.subplot(gs[0, 0])

        ax.tick_params(
            axis='both',
            which='both',
            bottom='on',
            top='on',
            right = 'on',
            labelbottom='on',
            labelleft='on')

        ax.set_aspect('equal')

        x, y = empty(self.length), empty(self.length)

        for i in range(self.length):
            x[i] = (self.aa_pos[i][1])
            y[i] = (self.length + 1 - self.aa_pos[i][0])

        ax.set_xlim(amin(x) - 1, amax(x) + 1)
        ax.set_ylim(amin(y) - 1, amax(y) + 1)
        #ax.set_ylim(5, 10)
        #ax.set_xlim(0, 11)
        #ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        #ax.set_yticks([5, 6, 7, 8, 9, 10])

        ax.set_axisbelow(True)

        ax.grid(color='lightgray', linestyle='dashed')

        label, = ax.plot(x, y, '-or', label=Label)

        legend = ax.legend(bbox_to_anchor=(0., 1.04, 1., .102), loc=3, ncol=1, mode="expand", borderaxespad=0., fontsize=17.3,
                  bbox_transform=ax.transAxes)

        plt.show()