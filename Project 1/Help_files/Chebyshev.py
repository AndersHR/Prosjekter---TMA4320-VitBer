from numpy import empty, cos, pi

def chebyshev_n(n):
    x=empty(n)
    for i in range (n):
        x[i]=cos(pi*(1+2*(n-1-i))/(2*n))
    return x

def shift_chebyshev(a, b, x):
    for i in range(len(x)):
        x[i] = ((b-a)/2)*x[i]+(b+a)/2
    return x