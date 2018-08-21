from numpy import exp


class Fitting(object):

    @classmethod
    def gauss_function(self, x, *p):
        """Define and return a Gaussian function.

        This function returns the value of a Gaussian function, using the
        A, mu and sigma value that is provided as *p.

        Keyword arguments:
        x -- number
        p -- A, mu and sigma numbers
        """
        A, mu, sigma = p
        return A*exp(-(x-mu)**2/(2.*sigma**2))
