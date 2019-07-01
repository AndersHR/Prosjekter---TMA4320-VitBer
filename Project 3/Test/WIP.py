from numpy import array, sqrt
from scipy.constants import pi

g=4*pi**2

def get_vals(planet):
    # Dict contains [Semimajor axis of orbit, Period, Eccentricity, Mass]
    interesting_shit = {'Mercury': array([0.3871, 0.2408, 0.2056, 0.0553 / 333480]),
                        'Venus': array([0.7233, 0.6152, 0.0068, 0.8250 / 333480]),
                        'Earth': array([1.0000, 1.0000, 0.0167, 1.0000 / 333480]),
                        'Mars': array([1.5237, 1.8809, 0.0934, 0.1074 / 333480]),
                        'Jupiter': array([5.2028, 11.862, 0.0483, 317.89 / 333480]),
                        'Saturn': array([9.5388, 29.456, 0.0560, 95.159 / 333480]),
                        'Uranus': array([19.191, 84.07, 0.0461, 14.56 / 333480]),
                        'Neptune': array([30.061, 164.81, 0.0100, 17.15 / 333480])}
    return interesting_shit[planet]


def getStartVals(planet, x, theta):
    vals = get_vals(planet)
    if (x == 1):
        r = array([vals[0] * (1 + vals[3]), 0])
        v = array([0, sqrt(g) * sqrt((1 + vals[3]) * (1 + vals[4]) / (vals[0] * (1 - vals[3])))])
        m = vals[4]
    elif (x == -1):
        r = array([vals[0] * (1 - vals[3]), 0])
        v = array([0, sqrt(g) * sqrt((1 - vals[3]) * (1 + vals[4]) / (vals[0] * (1 + vals[3])))])
        m = vals[4]
    else:
        a = vals[0]
        e = vals[3]
        c = a * e
        b = sqrt(a ** 2 - c ** 2)
        if (x > 1):
            x_new = ((a+c)*x-c)
        else:
            x_new = ((a-c)*x-c) #-((a-c)(-x)-c
        y = b*sqrt(1-x**2/a**2)


