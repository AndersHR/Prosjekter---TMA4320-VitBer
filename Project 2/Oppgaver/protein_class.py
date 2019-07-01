from numpy import empty, zeros, int8, array, exp, sqrt, float64, copy, append
from numpy.random import uniform, randint, seed
from scipy.constants import k
from matplotlib.pyplot import plot, subplot, show, legend, figure
import matplotlib.gridspec as grd
from matplotlib import pyplot as plt

class Protein:

    def __init__(self,length):
        self.grid = zeros((length+2, length+2), dtype=int8)
        for i in range(1, length+1):
            self.grid[int((length+1)/2)][i] = i
        self.aa_pos = array([[int((length +1)/ 2), i + 1] for i in range(length)]) #[y, x]
        self.backup_grid = copy(self.grid)
        self.aa_pos_backup = copy(self.aa_pos)


        self.potentials = zeros((length, length), dtype = float64)
        for i in range(length-2):
            self.potentials[i][i+2:] = uniform(-10.4e-21, -3.47e-21, length-2-i)
        for i in range(2, length):
            for j in range(i-1):
                self.potentials[i][j] = self.potentials[j][i]

        self.fixed_point = int((length +1)/ 2)
        self.u = 0.0
        self.d = length
        self.axis = 0
        self.direction = 1
        self.length = length
        self.rotation = 0
        self.T = 100

    def choose_axis(self):
        self.axis = randint(2,self.length)
        if self.axis < (self.fixed_point):
            self.direction = -1
        else:
            self.direction = 1
        self.rotation = randint(0, 2, 1)[0]
        if self.rotation == 0:
            self.rotation = -1

    def bend(self):
        if (self.direction==1):
            end = self.length+1
        else:
            end = 0
        for i in range(self.axis+self.direction, end, self.direction):
            self.grid[self.aa_pos_backup[i-1][0],self.aa_pos_backup[i-1][1]] = 0
        for i in range(self.axis+self.direction, end, self.direction):
            self.grid[self.aa_pos[i-1][0], self.aa_pos[i-1][1]] = i

    def check_for_collisions(self):
        if (self.direction==1):
            end = self.length+1
        else:
            end = 0
        for i in range(self.axis+self.direction, end, self.direction):
            new_pos = array([(self.aa_pos[i-1, 1]-self.aa_pos[self.axis-1, 1])*self.rotation + self.aa_pos[self.axis-1, 0],
                       (self.aa_pos[self.axis-1, 0] - self.aa_pos[i-1, 0])*self.rotation + self.aa_pos[self.axis-1, 1]])
            if (self.grid[new_pos[0],new_pos[1]] == 0 or self.grid[new_pos[0],new_pos[1]]*self.direction > self.axis*self.direction):
                self.aa_pos[i-1] = new_pos
            else:
                self.aa_pos = copy(self.aa_pos_backup)
                return False
        return True

    def compute_u(self):
        ny_u = 0
        for i in range(self.length):
            for j in range(-1, 2, 2):
                posy, posx = self.aa_pos[i][0], self.aa_pos[i][1]
                if (posy+j<self.length and posy+j>-1):
                    if (self.grid[posy+j, posx] != 0):
                        ny_u += self.potentials[i, self.grid[posy+j, posx]-1]
                if (posx+j<self.length and posx+j>-1):
                    if (self.grid[posy, posx+j] != 0):
                        ny_u += self.potentials[i, self.grid[posy, posx+j]-1]
        return ny_u/2

    def check_for_fluctuations(self):
        ny_u = self.compute_u()
        delta_u = ny_u - self.u
        if (delta_u > 0):
            p_ratio = exp(-delta_u/(k*self.T))
            fluctuations = uniform()
            if (p_ratio < fluctuations):
                self.aa_pos = copy(self.aa_pos_backup)
                self.grid = copy(self.backup_grid)
                return False
            
        self.backup_grid = copy(self.grid)
        self.aa_pos_backup = copy(self.aa_pos)
        self.u = ny_u
        return True

    def bend_protein(self):
        self.choose_axis()
        if (self.check_for_collisions()):
            self.bend()
            if (self.check_for_fluctuations()):
                return True
        return False

    def find_d(self):
        d1=0
        for i in range(self.length-1):
            for j in range(i+1, self.length):
                cur_d = (self.aa_pos[i][0]-self.aa_pos[j][0])**2 + (self.aa_pos[i][1]-self.aa_pos[j][1])**2
                if cur_d > d1:
                    d1 = cur_d
        self.d = sqrt(d1)

    def plot_protein(self):

        Label = "length = {}, E = {:.2E}, L = {:.2E}".format(self.length,self.u,0)

        ax = subplot(111)
        x, y = empty(self.length),empty(self.length)
        for i in range(self.length):
            x[i] = (self.aa_pos[i][1])
            y[i] = (self.length+1-self.aa_pos[i][0])

        ax.set_xlim(0,self.length+1)
        ax.set_ylim(0,self.length+2)

        ax.set_axisbelow(True)

        #ax.set_xticks([i for i in range(1,self.length + 1)])
        #ax.set_yticks([i for i in range(1,self.length + 2)])

        ax.yaxis.grid(color='gray', linestyle='dashed')
        ax.xaxis.grid(color='gray', linestyle='dashed')

        plot(x, y, '-or', label=Label)

        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=14,
                  bbox_transform=ax.transAxes)

        show()

        print ("RrrrRrRRRrrrrRRRRRRRRRRRRrrrrr")


    def __str__(self):
        print (self.grid)
        print (self.potentials)
        return ""

if __name__ == "__main__":
    seed(12425)                  #Interessante seeds: 1344(100), 41250(100)
    protein = Protein(300)
    e_vals = zeros(1)
    bends = 0
    turns = 0
    while (bends < 500):
        print (bends, turns)
        if protein.bend_protein():
            bends += 1
        turns += 1
        e_vals = append(e_vals, protein.u)
    print (protein.grid)
    print(" ")
    protein.plot_protein()

    x_vals = array([i for i in range(e_vals.shape[0])])

    fig = plt.figure('Figur E')
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    gs = grd.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.grid(color="lightgrey", linestyle = 'dashed')
    ax.invert_yaxis()
    ax.plot(x_vals, e_vals, 'r-')
    plt.show()








