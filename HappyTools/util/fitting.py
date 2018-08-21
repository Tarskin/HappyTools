from numpy import exp


class Fitting(object):

    def gauss_function(self, master):
        """Define and return a Gaussian function.

        This function returns the value of a Gaussian function, using the
        A, mu and sigma value that is provided as *p.

        Keyword arguments:
        x -- number
        p -- A, mu and sigma numbers
        """
        A, mu, sigma = master.p
        x = master.x
        return A*exp(-(x-mu)**2/(2.*sigma**2))
