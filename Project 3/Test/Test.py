from os.path import dirname

abs_path = dirname(dirname(__file__))

from numpy import load, copy

if __name__=="__main__":
    file = abs_path + r'/Data_files/Oppgave_1_keplers3rd.npz'
    file_npz = load(file)
    a = copy(file_npz['a_list'])
    tau = copy(file_npz['tau_list'])
    for i in range(a.shape[0]):
        print (tau[i]**2/a[i]**3)


