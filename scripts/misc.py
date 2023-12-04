"""
Misc functions

"""
import math
import numpy as np


def shell_volume(h, rho=17.5):
    """
    Volume function for an orbital shell. Takes the height as an 
    input and returns the volume of a spherical shell with thickness 35km. 
    The "rho" here is shell thickness, what is now \Delta in the SI.
    Units: km^3

    """
    volume = 4/3 * math.pi * (6371 + h + rho)**3 - 4/3 * math.pi * (6371 + h - rho)**3
    
    return round(volume)


def velocity(h):
    """
    # Velocity function for an orbital shell. Takes the height as an input and 
    # returns the velocity of an object at that altitude. From Wiki on orbital 
    # speed: "For orbits with small eccentricity, the length of the orbit is close 
    # to that of a circular one, and the mean orbital speed can be approximated 
    # either from observations of the orbital period and the semimajor axis of its 
    # orbit, or from knowledge of the masses of the two bodies and the semimajor 
    # axis." The formula given is v \approx sqrt(\mu/a), where \mu is the 
    # gravitational constant and a is the semimajor axis. The semimajor axis 
    # is the average distance from the center of the orbit to the center of the 
    # Earth, which is 6371 + h. The gravitational constant is 
    # 3.986004415e+5 km^3/s^2. (Source: https://en.wikipedia.org/wiki/Orbital_speed)
    # Units: km/s
    velocity <- function(h) {
        GM <- 3.986004415e+5 # gravitational constant in units of km^3/s^2
        sqrt(GM/(6371+h))
    }

    """
    gm = 3.986004415e+5 # gravitational constant in units of km^3/s^2

    velocity = np.sqrt(gm / (6371 + h))

    return velocity