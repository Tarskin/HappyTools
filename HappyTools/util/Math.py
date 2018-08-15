from numpy import exp

class Math(object):
    def __init__(self, master):
        self.master = master

    def gauss_function(self, master, x, *p):
        """Define and return a Gaussian function.

        This function returns the value of a Gaussian function, using the
        A, mu and sigma value that is provided as *p.

        Keyword arguments:
        x -- number
        p -- A, mu and sigma numbers
        """
        A, mu, sigma = p
        return A*exp(-(x-mu)**2/(2.*sigma**2))

    def powerLaw(self, master, x,a,b,c):
        """ TODO
        """
        penalty = 0
        if b > 2.:
            penalty = abs(b-1.)*10000
        if b < 0.:
            penalty = abs(2.-b)*10000
        return a*x**b + c + penalty
