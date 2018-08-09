from HappyTools.bin.Chromatogram import Chromatogram
import HappyTools.gui.Settings as settings
import HappyTools.gui.ProgressBar as progressbar
from numpy import std, mean, polyfit, poly1d
from bisect import bisect_left, bisect_right
import glob
import os
import sys

class Functions(object):

    def batch_process(self, cal_file, anal_file, batch_folder):
        if batch_folder.get():
            bar = progressbar.ProgressBar(self)

            if cal_file.get():
                calibrants = self.read_peak_list(cal_file.get())
                files = self.get_calibration_files(batch_folder)
                for index, file in enumerate(files):
                    data = Chromatogram(file)
                    self.calibrate_chrom(calibrants, data, batch_folder)
                    bar.update_progress_bar(bar.progressbar, bar.calibration_percentage, index, len(files))
            bar.update_progress_bar(bar.progressbar, bar.calibration_percentage, 1, 1)

            if anal_file.get():
                analytes = self.read_peak_list(anal.file.get())
                files = self.get_quantitation_files(batch_folder)
                for index, file in enumerate(files):
                    data = Chromatogram(file)
                    bar.update_progress_bar(bar.progressbar2, bar.quantitation_percentage, index, len(files))
                bar.update_progress_bar(bar.progressbar2, bar.quantitation_percentage, 1, 1)

    def get_calibration_files(self, batch_folder):
        """ TODO
        """
        calibration_files = []
        for files in settings.calibration_filetypes:
            for file in glob.glob(os.path.join(batch_folder.get(),files)):
                if file not in settings.exclusion_files:
                    calibration_files.append(os.path.join(batch_folder.get(),file))
        return calibration_files

    def get_quantitation_files(self, batch_folder):
        """ TODO
        """
        quantitation_files = []
        for files in settings.quantiation_filetypes:
            for file in glob.glob(os.path.join(batch_folder.get(),files)):
                if file not in settings.exclusion_files:
                    quantitation_files.append(os.path.join(batch_folder.get(),file))
        return quantitation_files

    def calibrate_chrom(self, calibrants, data, batch_folder):
        # TODO: Decide if baseline correction should be built in, or if
        # this is something that we expect the user to do manually and
        # check in the GUI.
        #data = Trace().baseline_correction(data)
        time_pairs = self.find_peak(calibrants, data)
        if len(time_pairs) >= settings.min_peaks:
            f = self.determine_calibration_function(time_pairs)
            data = self.apply_calibration_function(f, data)
            data.filename = os.path.join(batch_folder.get(), "calibrated_"+os.path.basename(data.filename))
        else:
            data.filename = os.path.join(batch_folder.get(), "uncalibrated_"+os.path.basename(data.filename))
        self.write_data(data)

    def quantify_chrom(self, master):
        pass

    def read_peak_list(self, file_name):
        """Read and parse the peak file and return a list of peaks.
        
        This function opens the file that is specified in 'fileName', and 
        reads it on a line per line basis. The function will split each 
        line on '\t' prior to trying to append the parts to the 'peaks'
        list. The function will write a message to the logFile if logging
        is enabled if the previous mentioned try goes to the except clause.
        
        Keyword argments:
        fileName -- string
        """
        peaks = []
        try:
            with open(file_name,'r') as fr:
                for line in fr:
                    line = line.rstrip("\n").split("\t")
                    try:
                        peaks.append((str(line[0]), float(line[1]), float(line[2])))
                    except ValueError:
                        self.log("Ignoring line: "+str(line)+" from file: "+str(file_name))

                        pass
        except IOError:
            tkMessageBox.showinfo("File Error","The selected reference file could not be opened.")
        return peaks

    def find_peak(self, reference, data):
        """ TODO
        """
        # TODO: Migrate rest of code to use commented version

        time_pairs = []
        #signal_noise = []

        time, intensity = zip(*data.data)
    
        for i in reference:

            low = bisect_left(time, i[1]-i[2])
            high = bisect_right(time, i[1]+i[2])

            low_background = bisect_left(time, max(i[1]-settings.background_window, settings.start))
            high_background = bisect_right(time, min(i[1]+settings.background_window, settings.end))

            if settings.background_noise_method == "NOBAN":
                NOBAN = self.noban(intensity[low_background:high_background])
            elif settings.background_noise_method == "MT":
                NOBAN = self.background_noise(intensity[low_background:high_background])

            max_value = max(intensity[low:high])
            max_index = intensity[low:high].index(max_value)

            if ((max_value - NOBAN['Background'])/NOBAN['Noise']) >= settings.min_peak_SN:
                time_pairs.append((i[1],time[low+max_index]))

            #time_pairs.append((i[1],time[low+max_index]))
            #signal_noise.append((max_value - NOBAN['Background'])/NOBAN['Noise'])
        
        return time_pairs
        
        #return {'time_pars': time_pairs, 'signal_noise': signal_noise}

    def background_noise(self, data):
        """Return the background and noise.

        This function determines the average and the standard deviation or
        the maximum difference of all segments of data, where each segment
        has the length specified in the slicepoints parameter.

        Keyword arguments:
        data -- list of intensities
        """
        background = sys.maxint
        currNoise = 0
        for index,i in enumerate(data[:-settings.slicepoints]):
            buffer = data[index:index+settings.slicepoints]
            if mean(buffer) < background:
                background = mean(buffer)
                if settings.noise == "MM":
                    currNoise = max(buffer)-min(buffer)
                elif settings.noise == "RMS":
                    currNoise = std(buffer)
        if currNoise == 0:
            currNoise = 1
        return {'Background': background, 'Noise': currNoise}

    def determine_calibration_function(self, time_pairs):
        """ TODO
        """
        expected, observed = zip(*time_pairs)
        if settings.use_UPC == "False":
            z = polyfit(observed, expected, 2)
            function = poly1d(z)
        elif settings.use_UPC == "True":
            raise NotImplementedError("Ultra performance calibration has not been implemented in the refactoring yet.")
        return function

    def apply_calibration_function(self, function, data):
        """ TODO
        """
        time, intensity = zip(*data.data)
        data.data = zip(function(time), intensity)
        return data

    def write_data(self, data):
        """ TODO
        """
        try:
            with open(data.filename,'w') as fw:
                for data_point in data.data:
                    fw.write(str(format(data_point[0],'0.'+str(settings.decimal_numbers)+'f'))+"\t"+str(format(data_point[1],'0.'+str(settings.decimal_numbers)+'f'))+"\n")
        except IOError:
            self.log("File: "+str(os.path.basename(data.filename))+" could not be opened.")

    def log(self, message):
        """ TODO
        """
        if debug.logging == True and debug.logLevel >= 1:
            with open(debug.logFile,'a') as fw:
                fw.write(str(datetime.now().replace(microsecond=0))+str(message)+"\n")
