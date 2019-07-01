from numpy import empty

#Lager et newton-polynom v.h.a. en liste X med x-komponenter
#og den tilhørende listen Y med y-komponenter
def get_newton_poly(X, Y):
    n = len(X)
    f_y = empty(n)
    f_y[0] = Y[0]
    f1 = Y
    for i in range(1, n):
        f2 = empty(n - i)
        for j in range(n - i):
            f2[j] = (f1[j + 1] - f1[j]) / (X[j + i] - X[j])
        f_y[i] = f2[0]
        f1 = f2
    return f_y

def evaluate_newton(X, f_y, x_to_eval): #Evaluerer newtonpolynom med noder X(array) med vekter f_y i punktet x_to_eval
    n = len(X)
    Ysum = 0
    for i in range(n):
        Ysum = Ysum * (x_to_eval - X[n - i-1]) + f_y[n - i-1]
    return Ysum

def get_lagrange_poly_weights(X, Y):
    n = len(X)
    r_y = empty(n)
    for i in range(n):
        inv_weight = 1
        for j in range(0, i):
            inv_weight = inv_weight*(X[i]-X[j])
        for j in range(i+1, n):
            inv_weight = inv_weight*(X[i]-X[j])
        r_y[i] = Y[i]/inv_weight
    return r_y

def evaluate_lagrange(x_l, x, n_l): #Tar inn nodene x_l til lagrangepolynomet som skal evalueres i x, og returner
    N_l = len(x_l)                    #og returnerer PI(i=/=l)(x-x_i)
    sum_l = 1
    for j in range(0, n_l):
        sum_l = sum_l*(x - x_l[j])
    for j in range(n_l + 1, N_l):
        sum_l = sum_l * (x - x_l[j])
    return sum_l


