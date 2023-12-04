"""
Misc functions

"""



def shell_volume(h, rho=17.5):
    """
    Volume function for an orbital shell. Takes the height as an 
    input and returns the volume of a spherical shell with thickness 35km. 
    The "rho" here is shell thickness, what is now \Delta in the SI.
    Units: km^3

    """
    return    {
        4/3 * np.pi * (6371+h+rho)^3 - 4/3 * pi * (6371+h-rho)^3
    }