from numpy import array, sqrt, dot, sin, cos
from scipy.constants import pi

g = 4*pi**2
G = 6.67*(10**(-11))

def der_r(v):
    return v

def der_v(primary_coord_diff, secondary_coord_diff):
    return (g * (primary_coord_diff) / ((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2) ** (3 / 2))

def der_v_SI(primary_coord_diff, secondary_coord_diff):
    return (G * (primary_coord_diff) / ((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2) ** (3 / 2))

def der_v_einstein(primary_coord_diff, secondary_coord_diff, alfa):
    return (g * (primary_coord_diff) / ((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2) ** (3 / 2))*(1 + alfa/((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2))

#Dem to under her trengs bare for den implisette løseren, btw
def der_v_der_primary(primary_coord_diff, secondary_coord_diff):
    return g*((1 - 3 * (primary_coord_diff) ** 2) / ((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2)) / \
           (((primary_coord_diff) ** 2 + (secondary_coord_diff) ** 2) ** (3 / 2))

def der_v_der_secondary(primary_coord_diference, secondary_cood_diff):
    return -g*(3*(primary_coord_diference*secondary_cood_diff))/(((primary_coord_diference)**2+(secondary_cood_diff)**2)**(5/2))

def get_vals(planet):
    #Dict contains [Semimajor axis of orbit, Period, Eccentricity, Mass]
    interesting_shit = {'Mercury':array([0.3871, 0.2408, 0.2056, 0.0553/333480]),
                        'Venus':  array([0.7233, 0.6152, 0.0068, 0.8250/333480]),
                        'Earth':  array([1.0000, 1.0000, 0.0167, 1.0000/333480]),
                        'Mars':   array([1.5237, 1.8809, 0.0934, 0.1074/333480]),
                        'Jupiter':array([5.2028, 11.862, 0.0483, 317.89/333480]),
                        'Saturn': array([9.5388, 29.456, 0.0560, 95.159/333480]),
                        'Uranus': array([19.191, 84.07,  0.0461, 14.56/333480]),
                        'Neptune':array([30.061, 164.81, 0.0100, 17.15/333480])}
    return interesting_shit[planet]

def get_start_vals(planet, x, theta):
    vals = get_vals(planet)
    m = vals[3]
    a = vals[0]
    e = vals[2]
    if (x == 1):
        r = array([a * (1 + vals[2]), 0])
        v = array([0.0, sqrt(g) * sqrt((1 - vals[2]) * (1 + m) / (vals[0] * (1 + vals[2])))])
    elif (x == -1):
        r = array([-vals[0] * (1 - vals[2]), 0])
        v = array([0.0, -sqrt(g) * sqrt((1 + vals[2]) * (1 + m) / (vals[0] * (1 - vals[2])))])
    else:
        a = vals[0]
        e = vals[2]
        c = a * e
        b = sqrt(a ** 2 - c ** 2)
        if (x > 0):
            x_center = ((a+c)*x-c)
        else:
            print ("TEST")
            x_center = ((a-c)*x-c) #-((a-c)(-x)-c
        y = b*sqrt(1-x_center**2/a**2)
        r = array([x_center+c, y])
        dy_dx = (x_center/y)*(b/a)**2
        v = sqrt(g)*sqrt(2/(sqrt(y**2+(x_center+c)**2))-1/a)*sqrt(1+m)\
            * sqrt(1/(1+dy_dx**2)) * array([-1, dy_dx])
    if (theta!=0):
        r = dot(array([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]]), r)
        v = dot(array([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]]), v)
    return r, v, m
    #r=[x, y], v=[v_x, v_y], m=m



