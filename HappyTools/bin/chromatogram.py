from HappyTools.bin.peak import Peak
from HappyTools.util.pdf import Pdf
from HappyTools.util.file_parser import file_parser
from bisect import bisect_left, bisect_right
from pathlib import PurePath
from numpy import polyfit, poly1d
from scipy.signal import savgol_filter
import operator


class Chromatogram(object):
    def __init__(self, master):
        self.filename = master.filename
        self.settings = master.settings
        self.master = master
        self.logger = master.logger
        self.axes = master.axes
        self.chrom_data = None
        self.peaks = []

    def baseline_correction(self):
        # Background determination
        background = []
        chunks = [self.chrom_data[x:x+self.settings.points] for x in range(
                  0, len(self.chrom_data), self.settings.points)]
        for j in chunks:
            buff1, buff2 = zip(*j)
            min_index, _ = min(enumerate(buff2), key=operator.itemgetter(1))
            if buff1[0] > self.settings.start and buff1[-1] < self.settings.end:
                background.append((buff1[min_index], buff2[min_index]))

        # Baseline function
        time, intensity = zip(*background)
        func = polyfit(time, intensity, self.settings.baseline_order)
        p = poly1d(func)

        # Transform
        time = [a for a, b in self.chrom_data]
        new_chrom_intensity = [b-p(a) for a, b in self.chrom_data]

        # Uplift
        low = bisect_left(time, self.settings.start)
        high = bisect_right(time, self.settings.end)
        offset = abs(min(min(new_chrom_intensity[low:high]), 0))

        self.chrom_data = list(zip(time,
                               [x+offset for x in new_chrom_intensity]))

    def calibrate_chromatogram(self):
        time, intensity = zip(*self.chrom_data)
        time = self.calibration_function(time)
        self.chrom_data = list(zip(time, intensity))

    def determine_calibration_timepairs(self):
        self.calibration_time_pairs = []

        for i in self.master.reference:

            self.peak_name = i[0]
            self.peak_time = i[1]
            self.peak_window = i[2]

            self.peak = Peak(self)
            self.peak.determine_background_and_noise()
            self.peak.determine_signal_noise()

            time, intensity = zip(*self.peak.peak_data[self.peak.low:self.peak.high])
            max_value = max(intensity)
            max_index = intensity.index(max_value)

            if self.peak.signal_noise >= self.settings.min_peak_SN:
                self.calibration_time_pairs.append((i[1], time[max_index]))

    def determine_calibration_function(self):
        expected, observed = zip(*self.calibration_time_pairs)
        z = polyfit(observed, expected, 2)
        self.calibration_function = poly1d(z)

    def generate_pdf_report(self):
        pdf = Pdf(self)
        pdf.plot_overview()
        for self.peak in self.peaks:
            pdf.plot_individual()
        pdf.close()

    def normalize_chromatogram(self):
        time, intensity = zip(*self.chrom_data)

        # Normalize to maximum intensity
        maximum = max(intensity[bisect_left(
            time, self.settings.start):bisect_right(
            time, self.settings.end)])
        normalized_intensity = [b/maximum for a, b, in self.chrom_data]

        self.chrom_data = list(zip(time, normalized_intensity))

    def open_chromatogram(self):
        file_parser(self)

    def plot_chromatogram(self):
        label = PurePath(self.filename).stem
        time, intensity = zip(*self.chrom_data)
        self.axes.plot(time, intensity, label=str(label))

    def quantify_chromatogram(self):
        time, _ = zip(*self.chrom_data)
        self.peaks = []

        # Iterate over peaks
        for i in self.master.reference:

            self.peak_name = i[0]
            self.peak_time = i[1]
            self.peak_window = i[2]

            self.peak = Peak(self)

            # Ignore peaks outside the Trace RT window
            if self.peak_time < self.settings.start+self.settings.background_window \
                    or self.peak_time > self.settings.end-self.settings.background_window:
                continue

            # Background correct
            self.peak.determine_background_and_noise()
            if self.peak.background < 0.:
                self.peak.background_correct()
                self.peak.determine_background_and_noise()

            # Determine Areas, background and noise
            self.peak.determine_peak_area()
            self.peak.determine_background_area()
            self.peak.determine_peak_noise()
            self.peak.determine_signal_noise()
            self.peak.determine_total_area()

            # Data Subset (based on second derivative)
            self.peak.determine_spline_and_derivative()
            self.peak.determine_breakpoints()
            self.peak.subset_data()

            # Gaussian fit
            self.peak.determine_gaussian_coefficients()
            if self.peak.coeff.size > 0:
                self.peak.determine_gaussian_area()
                self.peak.determine_gaussian_parameters()
                self.peak.determine_height()
            self.peak.determine_actual_time()
            self.peak.determine_residual()

            # Results
            self.peaks.append(self.peak)

    def smooth_chromatogram(self):
        # Apply Savitzky-Golay filter
        # TODO: Allow user to specify Smoothing method
        time, intensity = zip(*self.chrom_data)
        new = savgol_filter(intensity, 21, 3)
        self.chrom_data = list(zip(time, new))

    def save_chromatogram(self):
        with open(self.filename, 'w') as fw:
            for data_point in self.trace.chrom_data:
                fw.write(
                    str(format(data_point[0], '0.' +
                        str(self.settings.decimal_numbers)+'f'))+'\t' +
                    str(format(data_point[1], '0.' +
                        str(self.settings.decimal_numbers)+'f'))+'\n')

def finalize_plot(master):
    master.axes.set_xlabel('Time [m]')
    master.axes.set_ylabel('Intensity [au]')
    handles, labels = master.axes.get_legend_handles_labels()
    master.axes.legend(handles, labels)
    master.axes.get_xaxis().get_major_formatter().set_useOffset(False)
    master.canvas.draw()
