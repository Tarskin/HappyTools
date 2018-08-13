import HappyTools.gui.ProgressBar as progressbar
import HappyTools.gui.Version as version
import HappyTools.gui.Debug as debug
from HappyTools.gui.Settings import Settings
from HappyTools.bin.Chromatogram import Chromatogram
from datetime import datetime
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from numpy import std, mean, polyfit, poly1d, linspace, argmax, greater, less, array
from numpy import max as numpy_max
from numpy import exp # only till gauss function gets moved to math
from bisect import bisect_left, bisect_right
from os import path
from glob import glob
from sys import maxint
import math

class Functions(object):

    def apply_calibration_function(self, master, data):
        """ TODO
        """
        time, intensity = zip(*data.data)
        data.data = zip(master.function(time), intensity)
        return data

    def background_noise(self, master, data):
        """Return the background and noise.

        This function determines the average and the standard deviation or
        the maximum difference of all segments of data, where each segment
        has the length specified in the slicepoints parameter.

        Keyword arguments:
        data -- list of intensities
        """
        background = maxint
        currNoise = 0
        for index,i in enumerate(data[:-master.settings.slicepoints]):
            buffer = data[index:index+master.settings.slicepoints]
            if mean(buffer) < background:
                background = mean(buffer)
                if master.settings.noise == "MM":
                    currNoise = max(buffer)-min(buffer)
                elif master.settings.noise == "RMS":
                    currNoise = std(buffer)
        if currNoise == 0:
            currNoise = 1
        return {'Background': background, 'Noise': currNoise}

    def batch_process(self, master):
        if master.batch_folder.get():

            self.settings = master.settings
            bar = progressbar.ProgressBar(self)

            if master.cal_file.get():

                self.cal_file = master.cal_file
                self.batch_folder = master.batch_folder

                self.calibrants = self.read_peak_list(self.cal_file.get())
                files = self.get_calibration_files(master)
                for index, file in enumerate(files):
                    data = Chromatogram(file)
                    self.calibrate_chrom(self, data)
                    bar.update_progress_bar(bar.progressbar, bar.calibration_percentage, index, len(files))

            bar.update_progress_bar(bar.progressbar, bar.calibration_percentage, 1, 1)

            if master.anal_file.get():

                analytes = self.read_peak_list(master.anal_file.get())
                files = self.get_quantitation_files(master.batch_folder)

                for index, file in enumerate(files):
                    data = Chromatogram(file)
                    bar.update_progress_bar(bar.progressbar2, bar.quantitation_percentage, index, len(files))

                bar.update_progress_bar(bar.progressbar2, bar.quantitation_percentage, 1, 1)

    def determine_breakpoints(self, master, data):

        time, intensity = zip(*data)
        f = InterpolatedUnivariateSpline(time[master.low:master.high], intensity[master.low:master.high])
        f_prime = f.derivative()

        new_x = linspace(time[master.low], time[master.high], 2500*(time[master.high]-time[master.low]))
        new_prime_y = f_prime(new_x)

        maxm = argrelextrema(new_prime_y, greater)
        minm = argrelextrema(new_prime_y, less)

        breaks = maxm[0].tolist() + minm[0].tolist()
        breaks = sorted(breaks)

        return breaks

    def calibrate_chrom(self, master, data):
        # TODO: Decide if baseline correction should be built in, or if
        # this is something that we expect the user to do manually and
        # check in the GUI.
        
        # I broke it!
        
        #data = Trace().baseline_correction(data)
        self.time_pairs = self.find_peak(master, data)
        if len(self.time_pairs) >= master.settings.min_peaks:
            self.function = self.determine_calibration_function(self)
            data = self.apply_calibration_function(self, data)
            data.filename = path.join(master.batch_folder.get(), "calibrated_"+path.basename(data.filename))
        else:
            data.filename = path.join(master.batch_folder.get(), "uncalibrated_"+path.basename(data.filename))
        self.write_data(master, data)

    def combine_results(self, master):
        """
        """
        # Read the raw files and construct a data structure
        Results = []
        for file in glob(path.join(master.batch_folder.get(),"*.raw")):
            Buffer = []
            with open(file,'r') as fr:
                fr.readline()
                for line in fr:
                    chunks = line.rstrip('\n').split('\t')
                    Buffer.append({'Peak':str(chunks[0]),'Time':float(chunks[1]),'Area':float(chunks[2]),'S/N':float(chunks[3]),
                                'Background':float(chunks[4]),'Noise':float(chunks[5]),'Residual':float(chunks[6]),'PeakNoise':float(chunks[7]),
                                'BackgroundArea':float(chunks[8]),'ActualTime':float(chunks[9]),'fwhm':float(chunks[10])})
            #with open(os.path.join(batchFolder.get(),os.path.splitext(os.path.basename(file))[0]+".cal")) as fr:
            #    formula = fr.readline()
            #Results.append({'File':str(os.path.splitext(os.path.basename(file))[0]), 'Calibration':str(formula), 'Data':Buffer})
            Results.append({'File':str(path.splitext(path.basename(file))[0]), 'Data':Buffer})

        # Construct the filename for the output
        utc_datetime = datetime.utcnow()
        s = utc_datetime.strftime("%Y-%m-%d-%H%MZ")
        filename = s + "_" + master.settings.output

        # Construct header
        header = ""
        for i in Results:
            for j in i['Data']:
                header = header + "\t"+str(j['Peak'])
            header = header + "\n"
            for j in i['Data']:
                header = header + "\t"+str(j['Time'])
            header = header + "\n"
            break

        # Write results, settings and version information
        with open(path.join(master.batch_folder.get(),filename),'w') as fw:
            # Metadata
            fw.write("HappyTools Settings\n")
            fw.write("Version:\t"+str(version.version)+"\n")
            fw.write("Build:\t"+str(version.build)+"\n")
            fw.write("Start Time:\t"+str(master.settings.start)+"\n")
            fw.write("End Time:\t"+str(master.settings.end)+"\n")
            fw.write("Baseline Order:\t"+str(master.settings.baseline_order)+"\n")
            fw.write("Background Window:\t"+str(master.settings.background_window)+"\n")
            fw.write("Background and noise method:\t"+str(master.settings.background_noise_method)+"\n")
            if master.settings.background_noise_method == "MT":
                fw.write("MT Slice Points:\t"+str(master.settings.slicepoints)+"\n")
            elif master.settings.background_noise_method == "NOBAN":
                fw.write("NOBAN Initial Estimate:\t"+str(master.settings.noban_start)+"\n")
            fw.write("Noise:\t"+str(master.settings.noise)+"\n")
            fw.write("\n")

            # Area (non background subtracted)
            if master.abs_int.get() == 1 and master.bck_sub.get() == 0:
                fw.write("Peak Area")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['Area']))
                    fw.write("\n")
                fw.write("\n")

            # Area (Background subtracted)
            if  master.abs_int.get() == 1 and master.bck_sub.get() == 1:
                fw.write("Peak Area (Background Subtracted)")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(max(j['Area']-j['BackgroundArea'],0)))
                    fw.write("\n")
                fw.write("\n")

            # Relative Area
            if master.rel_int.get() == 1 and master.bck_sub.get() == 0:
                fw.write("Relative Peak Area (TAN)")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    total = 0.
                    for j in i['Data']:
                        total += j['Area']
                    for j in i['Data']:
                        try:
                            fw.write("\t"+str(j['Area']/total))
                        except ZeroDivisionError:
                            fw.write("\t"+str(0.0))
                    fw.write("\n")
                fw.write("\n")

            # Relative Area (Background subtracted)
            if master.rel_int.get() == 1 and master.bck_sub.get() == 1:
                fw.write("Relative Peak Area (TAN, Background Subtracted)")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    total = 0.
                    for j in i['Data']:
                        total += max(j['Area']-j['BackgroundArea'],0)
                    for j in i['Data']:
                        try:
                            fw.write("\t"+str(max(j['Area']-j['BackgroundArea'],0)/total))
                        except ZeroDivisionError:
                            fw.write("\t"+str(0.0))
                    fw.write("\n")
                fw.write("\n")

            # Peak Noise (standard deviation of the integration window)
            if master.bck_noise.get() == 1:
                fw.write("Peak Noise (standard deviation of integration window)")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    total = 0.
                    for j in i['Data']:
                        fw.write("\t"+str(j['PeakNoise']))
                    fw.write("\n")
                fw.write("\n")

            # Background
            if master.bck_noise.get() == 1:
                fw.write("Background")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['Background']))
                    fw.write("\n")
                fw.write("\n")

            # Noise
            if master.bck_noise.get() == 1:
                fw.write("Noise")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['Noise']))
                    fw.write("\n")
                fw.write("\n")

            # S/N
            if master.peak_qual.get() == 1:
                fw.write("Signal-to-Noise")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['S/N']))
                    fw.write("\n")
                fw.write("\n")

            # GPQ
            if master.peak_qual.get() == 1:
                fw.write("GPQ (Gaussian Peak Quality)")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['Residual']))
                    fw.write("\n")
                fw.write("\n")

            # FWHM
            if master.peak_qual.get() == 1:
                fw.write("FWHM")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(j['fwhm']))
                    fw.write("\n")
                fw.write("\n")

            # Tr residual
            if master.peak_qual.get() == 1:
                fw.write("Retention Time Residual")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    #if i['Calibration']:
                    #    fw.write(" ["+str(i['Calibration'])+"]")
                    for j in i['Data']:
                        residualTime = abs(float(j['ActualTime']) - float(j['Time']))
                        fw.write("\t"+str(residualTime))
                    fw.write("\n")
                fw.write("\n")

            # Peak Tr
            if master.peak_qual.get() == 1:
                fw.write("Retention Time")
                fw.write(header)
                for i in Results:
                    fw.write(i['File'])
                    for j in i['Data']:
                        fw.write("\t"+str(float(j['ActualTime'])))
                    fw.write("\n")
                fw.write("\n")

    def determine_calibration_function(self, master):
        """ TODO
        """
        expected, observed = zip(*master.time_pairs)
        if master.settings.use_UPC == "False":
            z = polyfit(observed, expected, 2)
            function = poly1d(z)
        elif master.settings.use_UPC == "True":
            raise NotImplementedError("Ultra performance calibration has not been implemented in the refactoring yet.")
        return function

    def determine_gaussian_parameters(self, coeff):
        """Calculate the FWHM.
        
        This function will calculate the FWHM based on the following formula
        FWHM = 2*sigma*sqrt(2*ln(2)). The function will return a dictionary
        with the fwhm ('fwhm'), the Gaussian peak center ('center') and the
        +/- width, from the peak center ('width').
        
        Keyword arguments:
        coeff -- coefficients as calculated by SciPy curve_fit
        """
        fwhm = abs(2*coeff[2]*math.sqrt(2*math.log(2)))
        width = 0.5*fwhm
        center = coeff[1]
        return {'fwhm':fwhm, 'width':width, 'center':center}

    def find_peak(self, master, data):
        """ TODO
        """
        # TODO: Migrate rest of code to use commented version

        time_pairs = []
        #signal_noise = []

        time, intensity = zip(*data.data)
    
        for i in master.calibrants:

            low = bisect_left(time, i[1]-i[2])
            high = bisect_right(time, i[1]+i[2])

            low_background = bisect_left(time, max(i[1]-master.settings.background_window, master.settings.start))
            high_background = bisect_right(time, min(i[1]+master.settings.background_window, master.settings.end))

            if master.settings.background_noise_method == "NOBAN":
                NOBAN = self.noban(master, intensity[low_background:high_background])
            elif master.settings.background_noise_method == "MT":
                NOBAN = self.background_noise(master, intensity[low_background:high_background])

            max_value = max(intensity[low:high])
            max_index = intensity[low:high].index(max_value)

            if ((max_value - NOBAN['Background'])/NOBAN['Noise']) >= master.settings.min_peak_SN:
                time_pairs.append((i[1],time[low+max_index]))

            #time_pairs.append((i[1],time[low+max_index]))
            #signal_noise.append((max_value - NOBAN['Background'])/NOBAN['Noise'])
        
        return time_pairs
        
        #return {'time_pars': time_pairs, 'signal_noise': signal_noise}

    def quantify_peak(self, master, data):
        """ TODO
        """

        peaks = []
        times = []
        signal_noises = []
        peak_areas = []
        background_areas = []
        peak_noises = []
        residuals = []
        backgrounds = []
        noises = []
        fwhms = []
        actual_times = []
        results = []

        time, intensity = zip(*data.data)

        for i in master.reference:
            peak_area = 0
            background_area = 0
            gauss_area = 0
            total_area = 0

            # Find Insertion Points
            self.low = bisect_left(time, i[1]-i[2])
            self.high = bisect_right(time, i[1]+i[2])

            self.low_background = bisect_left(time, max(i[1]-master.settings.background_window, master.settings.start))
            self.high_background = bisect_right(time, min(i[1]+master.settings.background_window, master.settings.end))

            # Determine Background and Noise
            if master.settings.background_noise_method == "NOBAN":
                NOBAN = self.noban(master, intensity[self.low_background:self.high_background])
            elif master.settings.background_noise_method == "MT":
                NOBAN = self.background_noise(master, intensity[self.low_background:self.high_background])

            max_value = max(intensity[self.low:self.high])

            # Peak and Background areas
            for index,j in enumerate(intensity[self.low:self.high]):
                try:
                    peak_area += max(j,0) * (time[self.low+index]-time[self.low+index-1])
                    background_area += max(NOBAN['Background'],0) * (time[self.low+index]-time[self.low+index-1])
                except IndexError:
                    continue
            peak_noise = std(intensity[self.low:self.high])

            # Horrible Fitting Section

            #breaks = self.determine_breakpoints(self, data)

            f = InterpolatedUnivariateSpline(time[self.low:self.high], intensity[self.low:self.high])
            f_prime = f.derivative()

            new_x = linspace(time[self.low], time[self.high], 2500*(time[self.high]-time[self.low]))
            new_y = f(new_x)
            new_prime_y = f_prime(new_x)

            maxm = argrelextrema(new_prime_y, greater)
            minm = argrelextrema(new_prime_y, less)

            breaks = maxm[0].tolist() + minm[0].tolist()
            breaks = sorted(breaks)
            
            # Initialize maxPoint, xData and yData
            #new_x = linspace(time[low], time[high], 2500*(time[high]-time[low])
            newX = new_x
            newY = new_y
            
            maxPoint = 0
            xData = newX
            yData = [x - NOBAN['Background'] for x in newY]

            # Subset the data
            # Region from newY[0] to breaks[0]
            try:
                if max(newY[0:breaks[0]]) > maxPoint:
                    maxPoint = max(newY[0:breaks[0]])
                    xData = newX[0:breaks[0]]
                    yData = [x - NOBAN['Background'] for x in newY[0:breaks[0]]]
            except IndexError:
                pass
            # Regions between breaks[x] and breaks[x+1]
            try:
                for index,j in enumerate(breaks):
                    if max(newY[breaks[index]:breaks[index+1]]) > maxPoint:
                        maxPoint = max(newY[breaks[index]:breaks[index+1]])
                        xData = newX[breaks[index]:breaks[index+1]]
                        yData = [x - max(NOBAN['Background'],0) for x in newY[breaks[index]:breaks[index+1]]]
            except IndexError:
                pass
            # Region from break[-1] to newY[-1]
            try:
                if max(newY[breaks[-1]:-1]) > maxPoint:
                    maxPoint = max(newY[breaks[-1]:-1])
                    xData = newX[breaks[-1]:-1]
                    yData = [x - NOBAN['Background'] for x in newY[breaks[-1]:-1]]
            except IndexError:
                pass

            # Gaussian fit on main points
            peak = xData[yData > exp(-0.5)*max(yData)]
            guess_sigma = 0.5*(max(peak) - min(peak))
            newGaussX = linspace(xData[0], xData[-1], 2500*(xData[-1]-xData[0]))
            p0 = [numpy_max(yData), xData[argmax(yData)],guess_sigma]
            try:
                coeff, var_matrix = curve_fit(self.gauss_function, xData, yData, p0)
                newGaussY = self.gauss_function(newGaussX, *coeff)
                newGaussY = [x + NOBAN['Background'] for x in newGaussY]
                for index,j in enumerate(intensity[self.low:self.high]):
                    gauss_area += max(self.gauss_function(time[self.low+index],*coeff),0) * (time[self.low+index]-time[self.low+index-1])

                fwhm = self.determine_gaussian_parameters(coeff)
                height = self.gauss_function(fwhm['center']+fwhm['width'], *coeff)+NOBAN['Background']

            except TypeError:
                self.log("Not enough data points to fit a Gaussian to peak: "+str(i[0]))

            except RuntimeError:
                self.log("Unable to determine residuals for peak: "+str(i[1]))

            # Determine Area
            for index,j in enumerate(intensity[self.low:self.high]):
                total_area += max(j-NOBAN['Background'],0) * (time[self.low+index]-time[self.low+index-1])

            # Determine Residual
            try:
                if gauss_area != 0:
                    residual = min(gauss_area / total_area, 1.0)
            except ZeroDivisionError:
                pass

            # Append
            peaks.append(i[0])
            times.append(i[1])
            peak_areas.append(peak_area)
            signal_noises.append((max_value - NOBAN['Background'])/NOBAN['Noise'])
            background_areas.append(background_area)
            peak_noises.append(peak_noise)
            residuals.append(residual)
            backgrounds.append(NOBAN['Background'])
            noises.append(NOBAN['Noise'])
            fwhms.append(fwhm['fwhm'])
            actual_times.append(fwhm['center'])

            results.append({'peak': i[0], 'time': i[1], 'peak_area': peak_area, 'signal_noise': (max_value - NOBAN['Background'])/NOBAN['Noise'], 
                'background_area': background_area, 'peak_noise': peak_noise,
                'residual': residual, 'background': NOBAN['Background'], 'noise': NOBAN['Noise'], 
                'fwhm': fwhm['fwhm'], 'actual_time': fwhm['center']})
        
        return results
        #return {'peak': peaks, 'time': times, 'peak_areas': peak_areas, 'signal_noises': signal_noises, 
                #'background_areas': background_areas, 'peak_noises': peak_noises,
                #'residuals': residuals, 'backgrounds': backgrounds, 'noises': noises, 
                #'actual_times': actual_times, 'fwhms': fwhms, 'actual_times': actual_times}

    def get_calibration_files(self, master):
        """ TODO
        """
        calibration_files = []
        for files in master.settings.calibration_filetypes:
            for file in glob(path.join(master.batch_folder.get(),files)):
                if file not in master.settings.exclusion_files:
                    calibration_files.append(path.join(master.batch_folder.get(),file))
        return calibration_files

    def get_quantitation_files(self, master):
        """ TODO
        """
        quantitation_files = []
        for files in master.settings.quantiation_filetypes:
            for file in glob(path.join(master.batch_folder.get(),files)):
                if file not in master.settings.exclusion_files:
                    quantitation_files.append(path.join(master.batch_folder.get(),file))
        return quantitation_files

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
        
    def log(self, message):
        """ TODO
        """
        if debug.logging == True and debug.logLevel >= 1:
            with open(debug.logFile,'a') as fw:
                fw.write(str(datetime.now().replace(microsecond=0))+str(message)+"\n")

    def noban(self, master):
        raise NotImplementedError("This feature is not implemented in the refactor yet.")

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

    def write_data(self, master, data):
        """ TODO
        """
        try:
            with open(data.filename,'w') as fw:
                for data_point in data.data:
                    fw.write(str(format(data_point[0],'0.'+str(master.settings.decimal_numbers)+'f'))+"\t"+str(format(data_point[1],'0.'+str(master.settings.decimal_numbers)+'f'))+"\n")
        except IOError:
            self.log("File: "+str(path.basename(data.filename))+" could not be opened.")

    def write_results(self, master, data, results):
        """ TODO
        """
        result_file = path.join(master.batch_folder.get(), path.splitext(path.basename(data.filename))[0]+".raw")
        with open(result_file,'w') as fw:
            fw.write("Name\tTime\tPeak Area\tS/N\tBackground\tNoise\tGaussian Residual RMS\tPeak Noise\tBackground Area\tPeak Time\tFWHM\n")
            for i in results:
                fw.write(str(i['peak'])+"\t"+str(i['time'])+"\t"+str(i['peak_area'])+"\t"+str(i['signal_noise'])+"\t"+str(i['background'])+"\t"+
                    str(i['noise'])+"\t"+str(i['residual'])+"\t"+str(i['peak_noise'])+"\t"+str(i['background_area'])+"\t"+
                    str(i['actual_time'])+"\t"+str(i['fwhm'])+"\n")
