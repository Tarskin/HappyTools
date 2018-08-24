from HappyTools.util.fitting import Fitting
from bisect import bisect_left, bisect_right
from scipy.optimize import curve_fit
from numpy import argmax, array, linspace, mean, std
from numpy import max as numpy_max
from numpy import exp # Temporary
from sys import maxint
from math import sqrt, log


class Peak(object):
    def __init__(self, master):
        self.settings = master.settings

        self.peak = master.peak
        self.time = master.time
        self.window = master.window
        self.peak_area = 0.
        self.gaussian_area = 0.
        self.signal_noise = 0.
        self.background_area = 0.
        self.peak_noise = 0.
        self.residual = 0.
        self.background = 0.
        self.noise = 0.
        self.fwhm = 0.
        self.width = 0.
        self.height = 0.
        self.center = 0.
        self.actual_time = 0.
        self.total_area = 0.
        self.coeff = 0.

        time, _ = zip(*master.data.data)
        self.low = bisect_left(time, self.time-self.window)
        self.high = bisect_right(time, self.time+self.window)

        self.low_background = bisect_left(time, max(
            self.time-self.settings.background_window, self.settings.start))
        self.high_background = bisect_right(time, min(
            self.time+self.settings.background_window, self.settings.end))

    def background_correct(self, master):
        time, intensity = zip(*master.data.data)
        intensity = [x+abs(self.background) for x in intensity]
        master.data.data = zip(time, intensity)

    def determine_background_and_noise(self, master):
        _, intensity = zip(*master.data.data[self.low_background:
            self.high_background])

        if self.settings.background_noise_method == "NOBAN":
            raise NotImplementedError("This feature is not implemented in "+
                "the refactor yet.")

        elif self.settings.background_noise_method == "MT":
            background = maxint
            noise = 0
            for index, i in enumerate(intensity[:-self.settings.slicepoints]):
                buffer = intensity[index:index+self.settings.slicepoints]
                if mean(buffer) < background:
                    background = mean(buffer)
                    if self.settings.noise == "MM":
                        noise = max(buffer)-min(buffer)
                    elif self.settings.noise == "RMS":
                        noise = std(buffer)
            if noise == 0:
                noise = 1

        self.background = background
        self.noise = noise

    def determine_background_area(self, master):
        background_area = 0
        time, intensity = zip(*master.data.data)
        for index, j in enumerate(intensity[self.low:self.high]):
            try:
                background_area += max(self.background, 0) * (
                    time[self.low+index]-time[self.low+index-1])
            except IndexError:
                continue

        self.background_area = background_area

    def determine_gaussian_area(self, master):
        time, intensity = zip(*master.data.data)
        gaussian_area = 0.

        for index, j in enumerate(intensity[self.low:self.high]):
            gaussian_area += max(Fitting().gauss_function(time[self.low+index],
                *self.coeff), 0) * (time[self.low+index]-
                time[self.low+index-1])

        self.gaussian_area = gaussian_area

    def determine_gaussian_coefficients(self, master):
        coeff = []

        x_data, y_data = zip(*master.data_subset)
        peak = array(x_data)[y_data > exp(-0.5)*max(y_data)]
        guess_sigma = 0.5*(max(peak) - min(peak))
        new_gauss_x = linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-
            x_data[0]))
        p0 = [numpy_max(y_data), x_data[argmax(y_data)], guess_sigma]
        try:
            coeff, var_matrix = curve_fit(Fitting().gauss_function, x_data, y_data,
                p0)

        except TypeError:
            master.log("Not enough data points to fit a Gaussian to peak: "+
                str(self.peak))

        except RuntimeError:
            master.log("Unable to determine residuals for peak: "+
                str(self.peak))

        self.coeff = coeff

    def determine_gaussian_parameters(self, master):
        """Calculate the FWHM.

        This function will calculate the FWHM based on the following formula
        FWHM = 2*sigma*sqrt(2*ln(2)). The function will return a dictionary
        with the fwhm ('fwhm'), the Gaussian peak center ('center') and the
        +/- width, from the peak center ('width').

        Keyword arguments:
        coeff -- coefficients as calculated by SciPy curve_fit
        """
        fwhm = abs(2*self.coeff[2]*sqrt(2*log(2)))
        width = 0.5*fwhm
        center = self.coeff[1]

        self.fwhm = fwhm
        self.width = width
        self.center = center

    def determine_height(self, master):
        edge = self.center+self.width
        height = Fitting().gauss_function(edge, *self.coeff)

        self.height = height

    def determine_peak_area(self, master):
        peak_area = 0.
        time, intensity = zip(*master.data.data)

        for index, j in enumerate(intensity[self.low:self.high]):
            try:
                peak_area += max(j, 0) * (time[self.low+index]-time[
                    self.low+index-1])
            except IndexError:
                continue

        self.peak_area = peak_area

    def determine_peak_noise(self, master):
        time, intensity = zip(*master.data.data)
        peak_noise = std(intensity[self.low:self.high])

        self.peak_nose = peak_noise

    def determine_residual(self, master):
        residual = 0.

        try:
            if self.gaussian_area != 0:
                residual = min(self.gaussian_area / self.total_area, 1.0)
        except ZeroDivisionError:
            pass

        self.residual = residual

    def determine_signal_noise(self, master):
        time, intensity = zip(*master.data.data)
        maximum_point = max(intensity[self.low:self.high])
        signal_noise = (maximum_point - self.background) / self.noise

        self.signal_noise = signal_noise

    def determine_total_area(self, master):
        total_area = 0.
        time, intensity = zip(*master.data.data)

        for index, j in enumerate(intensity[self.low:self.high]):
            total_area += max(j-self.background, 0) * (time[self.low+index]-
                time[self.low+index-1])

        self.total_area = total_area
