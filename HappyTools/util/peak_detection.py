from HappyTools.util.fitting import Fitting
from HappyTools.bin.peak import Peak
import logging
import tkinter.filedialog as filedialog
from numpy import linspace
from pathlib import Path
from bisect import bisect_left, bisect_right


class PeakDetection(object):
    def __init__(self, master):
        self.settings = master.settings
        self.logger = logging.getLogger(__name__)
        self.functions = master.functions
        self.chrom = master.chrom
        self.detected_peaks = []

    def detect_peaks(self, master):

        orig_time, orig_intensity = zip(*self.chrom.trace.chrom_data)
        curr_intensity = orig_intensity

        time_start = bisect_left(orig_time, self.settings.start)
        time_end = bisect_right(orig_time, self.settings.end)

        # Loop till intensity falls below specified threshold
        while (max(curr_intensity[time_start:time_end]) >
               self.settings.peak_detection_min * max(
               orig_intensity[time_start:time_end])):

            self.window = (self.settings.end - self.settings.start) / 2
            self.time = self.window + self.settings.start

            # Determine breaks and get subset of data
            self.breaks = self.functions.determine_breakpoints(self)
            self.data_subset = self.functions.subset_data(self)

            # Get time and intensity lists
            x_data, _ = zip(*self.data_subset)

            # Create Peak() object
            self.peak = None
            self.window = (x_data[-1] - x_data[0]) / 2
            self.time = self.window + x_data[0]
            self.peak = Peak(self)

            # Gaussian fit
            self.peak.determine_gaussian_coefficients(self)
            if self.peak.coeff.any():
                new_intensity = []
                for index, i in enumerate(curr_intensity):
                    new_intensity.append(i - Fitting().gauss_function(
                                         orig_time[index], *self.peak.coeff))
                curr_intensity = new_intensity

            # Subtract Gaussian intensity from the current intensity
            self.chrom.trace.chrom_data = list(zip(orig_time, curr_intensity))

            # Create Gaussian data at 3 sigma width
            gauss_start = self.time-3*self.peak.coeff[2]
            gauss_end = self.time+3*self.peak.coeff[2]
            gauss_time = linspace(gauss_start, gauss_end, (gauss_end-
                                  gauss_start)*1000)
            gauss_intensity = Fitting().gauss_function(gauss_time,
                                                       *self.peak.coeff)

            # Store detected peak
            self.detected_peaks.append(zip(gauss_time, gauss_intensity))

        # Restore original data
        self.chrom.trace.chrom_data = list(zip(orig_time, orig_intensity))

    def plot_peaks(self, master):
        for index, peak in enumerate(self.detected_peaks):
            x_array, y_array = zip(*peak)
            master.axes.fill_between(x_array, 0, y_array, alpha=0.5,
                              label="Peak "+str(index+1))

    def write_peaks(self, master):
        save_file = filedialog.asksaveasfilename(title="Save Annotation File")
        with Path(save_file).open('w') as fw:
            fw.write("foo")
        
