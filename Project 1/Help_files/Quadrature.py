from numpy import arange, array, empty
from numpy.polynomial.legendre import leggauss

def midpoint(a, b, n):
    x_q = arange(a+(b-a)/(2*n), b, (b-a)/n)
    w_q = array([(b-a)/n for i in range(n)])
    return x_q, w_q

def simpson(a, b, n):
    x_q = arange(a, b+(b-a)/(3*n), (b-a)/(2*n))
    w_q = empty(2*n+1)
    w_q[0], w_q[2*n] = 1, 1
    for i in range(1, 2*n, 2):
        w_q[i] = 4
    for i in range(2, 2*n-1, 2):
        w_q[i] = 2
    w_q = w_q*(b-a)/(6*n)
    return x_q, w_q

def legendre_gauss(a, b, N_q):
    x_lg, w_lg = leggauss(N_q)
    x_lg = ((b + a) / 2) + ((b - a) / 2) * x_lg       #Fikset livet ditt Vemund,takk meg senere
    w_lg = w_lg*(b-a)/2
    return x_lg ,w_lg

def compute_legendre_gauss_quadratures(a,b,N_q,f):
    x,w = legendre_gauss(a,b,N_q)
    quadrature_sum = 0
    for i in range(N_q):
        quadrature_sum += w[i]*f(x[i])
    return quadrature_sum


if __name__ == "__main__":
    a, b = legendre_gauss(0, 1, 5)
    print (a, b)