from HappyTools.util.fitting import gauss_function
from bisect import bisect_left, bisect_right, bisect
from math import sqrt, log
from numpy import amax, argmax, array, exp, greater, less, linspace, mean, std
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from sys import maxsize


class Peak(object):
    def __init__(self, master):
        self.master = master
        self.settings = master.settings
        self.logger = master.logger

        self.peak_name = master.peak_name
        self.peak_time = master.peak_time
        self.peak_window = master.peak_window
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
        self.coeff = array([])
        self.peak_data = None
        self.first_derivative_data = None
        self.univariate_spline_data = None
        self.breaks = []
        self.peak_maximum_data = None

        # Inherit full data
        time, _ = zip(*master.chrom_data)

        low_background = bisect_left(time, max(
            self.peak_time-max(self.settings.background_window,
            self.peak_window), self.settings.start))
        high_background = bisect_right(time, min(
            self.peak_time+max(self.settings.background_window,
            self.peak_window), self.settings.end))

        # Inherit only the required data (based on background window)
        time, intensity = zip(*master.chrom_data[low_background:high_background])

        self.low = bisect_left(time, self.peak_time-self.peak_window)
        self.high = bisect_right(time, self.peak_time+self.peak_window)

        self.peak_data = list(zip(time, intensity))

    def background_correct(self):
        time, intensity = zip(*self.peak_data)
        intensity = [x-self.background for x in intensity]
        self.peak_data = list(zip(time, intensity))

    def determine_actual_time(self):
        time, intensity = zip(*self.peak_data)
        if self.coeff.any():
            intensity = gauss_function(time, *self.coeff)
            intensity = intensity.tolist()

        max_intensity_index = intensity.index(max(intensity))
        self.actual_time = time[max_intensity_index]

    def determine_background_and_noise(self):
        _, intensity = zip(*self.peak_data)

        if self.settings.background_noise_method == 'NOBAN':
            raise NotImplementedError('This feature is not implemented in '+
                'the refactor yet.')

        elif self.settings.background_noise_method == 'MT':
            background = maxsize
            noise = 0
            for index, _ in enumerate(intensity[:-self.settings.slicepoints]):
                buffer = intensity[index:index+self.settings.slicepoints]
                if mean(buffer) < background:
                    background = mean(buffer)
                    if self.settings.noise == 'MM':
                        noise = max(buffer)-min(buffer)
                    elif self.settings.noise == 'RMS':
                        noise = std(buffer)
            if noise == 0:
                noise = 1

        self.background = background
        self.noise = noise

    def determine_breakpoints(self):
        time, intensity = zip(*self.first_derivative_data)

        maxm = argrelextrema(array(intensity), greater)
        minm = argrelextrema(array(intensity), less)

        breaks = maxm[0].tolist() + minm[0].tolist()
        breaks = sorted(breaks)
        self.breaks = breaks

    def determine_background_area(self):
        background_area = 0
        time, intensity = zip(*self.peak_data)
        for index, _ in enumerate(intensity[self.low:self.high]):
            try:
                background_area += max(self.background, 0) * (
                    time[self.low+index]-time[self.low+index-1])
            except IndexError:
                continue

        self.background_area = background_area

    def determine_gaussian_area(self):
        time, intensity = zip(*self.peak_data)
        gaussian_area = 0.

        for index, _ in enumerate(intensity[self.low:self.high]):
            gaussian_area += max(gauss_function(time[self.low+index],
                *self.coeff), 0) * (time[self.low+index]-
                time[self.low+index-1])

        self.gaussian_area = gaussian_area

    def determine_gaussian_coefficients(self):
        x_data, y_data = zip(*self.peak_maximum_data)
        peak = array(x_data)[y_data > exp(-0.5)*max(y_data)]
        guess_sigma = 0.5*(max(peak) - min(peak))
        p0 = [amax(y_data), x_data[argmax(y_data)], guess_sigma]
        try:
            coeff, _ = curve_fit(gauss_function, x_data, y_data, p0)
            self.coeff = coeff

        except TypeError:
            self.logger.warn('Not enough data points to fit a Gaussian to '+
                             'peak: '+str(self.peak))

        except RuntimeError:
            self.logger.error('Unable to determine residuals for peak: '+
                              str(self.peak))

    def determine_gaussian_parameters(self):
        '''Calculate the FWHM.

        This function will calculate the FWHM based on the following formula
        FWHM = 2*sigma*sqrt(2*ln(2)). The function will return a dictionary
        with the fwhm ('fwhm'), the Gaussian peak center ('center') and the
        +/- width, from the peak center ('width').

        Keyword arguments:
        coeff -- coefficients as calculated by SciPy curve_fit
        '''
        fwhm = abs(2*self.coeff[2]*sqrt(2*log(2)))
        width = 0.5*fwhm
        center = self.coeff[1]

        self.fwhm = fwhm
        self.width = width
        self.center = center

    def determine_height(self):
        edge = self.center+self.width
        height = gauss_function(edge, *self.coeff)

        self.height = height

    def determine_peak_area(self):
        peak_area = 0.
        time, intensity = zip(*self.peak_data)

        for index, j in enumerate(intensity[self.low:self.high]):
            try:
                peak_area += max(j, 0) * (time[self.low+index]-
                    time[self.low+index-1])
            except IndexError:
                continue

        self.peak_area = peak_area

    def determine_peak_noise(self):
        _, intensity = zip(*self.peak_data)
        peak_noise = std(intensity[self.low:self.high])

        self.peak_noise = peak_noise

    def determine_residual(self):
        residual = 0.

        try:
            if self.gaussian_area != 0:
                residual = min(self.gaussian_area / self.total_area, 1.0)
        except ZeroDivisionError:
            pass

        self.residual = residual

    def determine_signal_noise(self):
        _, intensity = zip(*self.peak_data)
        maximum_point = max(intensity[self.low:self.high])
        signal_noise = (maximum_point - self.background) / self.noise

        self.signal_noise = signal_noise

    def determine_spline_and_derivative(self):
        time, intensity = zip(*self.peak_data)
        low = bisect_left(time, self.peak_time-self.peak_window)
        high = bisect_right(time, self.peak_time+self.peak_window)

        # Failsafe
        if high == len(time):
            high =- 1

        new_x = linspace(time[low], time[high],
                         int(2500*(time[high]-time[low])))

        f = InterpolatedUnivariateSpline(time[low:high],
                                         intensity[low:high])
        f_prime = f.derivative()

        self.univariate_spline_data = list (zip(new_x, f(new_x)))
        self.first_derivative_data = list(zip(new_x, f_prime(new_x)))

    def determine_total_area(self):
        total_area = 0.
        time, intensity = zip(*self.peak_data)

        for index, j in enumerate(intensity[self.low:self.high]):
            total_area += max(j-self.background, 0) * (time[self.low+index]-
                time[self.low+index-1])

        self.total_area = total_area

    def subset_data(self):
        # Todo, use bisect_left and bisect_right here for consistency
        time, intensity = zip(*self.univariate_spline_data)

        max_intensity = max(intensity)
        max_value = intensity.index(max_intensity)

        insert = bisect(self.breaks, max_value)

        if insert == 0:
            start = 0
        else:
            start = self.breaks[insert-1]

        try:
            end = self.breaks[insert]
        except IndexError:
            end = None

        time = time[start:end]
        intensity = intensity[start:end]

        self.peak_maximum_data = list(zip(time, intensity))
