from os import path, getcwd
from HappyTools.util.functions import create_tooltip
import tkinter as tk


class Settings(object):
    def __init__(self, master):
        self.master = master

        self.points = 100
        self.start = 10
        self.end = 60
        self.baseline_order = 1
        self.background_window = 1
        self.noban_start = 0.25
        self.slicepoints = 5
        self.peak_detection_min = 0.05
        self.peak_detection_edge = 'Sigma'
        self.peak_detection_edge_value = 2.0
        self.create_figure = 'True'
        self.min_peaks = 4
        self.min_peak_SN = 27
        self.use_UPC = 'False'

        self.decimal_numbers = 6
        self.min_improvement = 0.05
        self.use_interpolation = False
        self.noise = 'RMS'
        self.background_noise_method = 'MT'

        self.output = 'summary.results'
        self.date_format = '%Y-%m-%d-%H%MZ'
        self.settings = 'HappyTools.ini'

        self.exclusion_files = ['LICENSE.txt', 'CHANGELOG.txt']
        self.calibration_filetypes = ['*.txt', '*.arw', '*.csv']
        self.quantitation_filetypes = ['calibrated*.txt']

    def read_settings(self, master):
        '''Read the settings file.

        This function opens the settings file (default is HappyTools.ini),
        parses the lines of the settings file and takes the value from the
        settings file as a value for the changeable settings (e.g. the start
        variable can be read from the settings file).

        Keyword arguments:
        none
        '''
        with open(path.join(getcwd(), master.settings), 'r') as fr:
            for line in fr:
                chunks = line.strip('\n').split('\t')
                if chunks[0] == 'points:':
                    master.points = int(chunks[1])
                elif chunks[0] == 'start:':
                    master.start = float(chunks[1])
                elif chunks[0] == 'end:':
                    master.end = float(chunks[1])
                elif chunks[0] == 'baselineOrder:':
                    master.baseline_order = int(chunks[1])
                elif chunks[0] == 'backgroundWindow:':
                    master.background_window = float(chunks[1])
                elif chunks[0] == 'noise:':
                    master.noise = str(chunks[1])
                elif chunks[0] == 'slicepoints:':
                    master.slicepoints = int(chunks[1])
                elif chunks[0] == 'createFigure:':
                    master.create_figure = str(chunks[1])
                elif chunks[0] == 'minPeaks:':
                    master.min_peaks = int(chunks[1])
                elif chunks[0] == 'minPeakSN:':
                    master.min_peak_SN = int(chunks[1])
                elif chunks[0] == 'peakDetectionMin:':
                    master.peak_detection_min = float(chunks[1])
                elif chunks[0] == 'peakDetectionEdge:':
                    master.peak_detection_edge = str(chunks[1])
                elif chunks[0] == 'peakDetectionEdgeValue:':
                    master.peak_detection_edge_value = float(chunks[1])

    def settings_popup(self, master):
        ''' TODO: Redesign settings window (it's f'ing ugly)
        '''

        figure_variable = tk.StringVar()
        figure_variable.set(master.create_figure)
        peak_detection_edge_var = tk.StringVar()
        peak_detection_edge_var.set(master.peak_detection_edge)
        ultra_performance_calibration_var = tk.StringVar()
        ultra_performance_calibration_var.set(master.use_UPC)
        options = ['True', 'False']

        def close(self, master):
            ''' TODO
            '''
            master.points = int(points_window.get())
            master.start = float(start_window.get())
            master.end = float(end_window.get())
            master.baseline_order = int(baseline_order_window.get())
            master.background_window = float(background_window_window.get())
            master.slicepoints = int(slicepoints_window.get())
            master.create_figure = str(figure_variable.get())
            master.min_peaks = int(min_peak_window.get())
            master.min_peak_SN = int(min_peak_SN_window.get())
            master.peak_detection_min = float(peak_detection.get())
            master.peak_detection_edge = str(peak_detection_edge_var.get())
            master.peak_detection_edge_value = float(
                peak_detection_edge_value_window.get())
            master.use_UPC = str(ultra_performance_calibration_var.get())
            top.destroy()

        def save(self, master):
            ''' TODO
            '''
            with open(path.join(getcwd(), master.settings), 'w') as fw:
                fw.write('points:\t'+str(int(points_window.get()))+'\n')
                fw.write('start:\t'+str(float(start_window.get()))+'\n')
                fw.write('end:\t'+str(float(end_window.get()))+'\n')
                fw.write('baselineOrder:\t'+str(int(
                    baseline_order_window.get()))+'\n')
                fw.write('backgroundWindow:\t'+str(float(
                    background_window_window.get()))+'\n')
                fw.write('noise:\t'+str(master.noise)+'\n')
                fw.write('slicepoints:\t'+str(slicepoints_window.get())+'\n')
                fw.write('createFigure:\t'+str(figure_variable.get())+'\n')
                fw.write('minPeaks:\t'+str(min_peak_window.get())+'\n')
                fw.write('minPeakSN:\t'+str(min_peak_SN_window.get())+'\n')
                fw.write('peakDetectionMin:\t'+str(peak_detection.get())+'\n')
                fw.write('peakDetectionEdge:\t'+str(
                    peak_detection_edge_var.get())+'\n')
                fw.write('peakDetectionEdgeValue:\t'+str(
                    peak_detection_edge_value_window.get())+'\n')
                fw.write('useUPC:\t'+str(
                    ultra_performance_calibration_var.get())+'\n')

        top = tk.top = tk.Toplevel()
        top.title('Settings')
        top.protocol('WM_DELETE_WINDOW', lambda: close(self, master))

        # General Settings
        general_label = tk.Label(top, text='General Settings', font=(
                                 'Helvetica', 16))
        general_label.grid(row=0, columnspan=2, sticky=tk.W)

        start_label = tk.Label(top, text='Start Time', font=('Helvetica', 12))
        start_label.grid(row=1, column=0, sticky=tk.W)
        start_window = tk.Entry(top)
        start_window.insert(0, master.start)
        start_window.grid(row=1, column=1, sticky=tk.W)

        end_label = tk.Label(top, text='End Time', font=('Helvetica', 12))
        end_label.grid(row=2, column=0, sticky=tk.W)
        end_window = tk.Entry(top)
        end_window.insert(0, master.end)
        end_window.grid(row=2, column=1, sticky=tk.W)

        # Peak Detection Settings
        peak_detection_label = tk.Label(top, text='Peak Detection Settings',
                                        font=('Helvetica', 16))
        peak_detection_label.grid(row=3, columnspan=2, sticky=tk.W)

        peak_detection_label = tk.Label(top, text='Minimum Intensity',
                                        font=('Helvetica', 12))
        peak_detection_label.grid(row=4, column=0, sticky=tk.W)
        peak_detection = tk.Entry(top)
        peak_detection.insert(0, master.peak_detection_min)
        peak_detection.grid(row=4, column=1, sticky=tk.W)

        peak_detection_edge_label = tk.Label(top, text='Edge Method',
                                             font=('Helvetica', 12))
        peak_detection_edge_label.grid(row=5, column=0, sticky=tk.W)
        peak_detection_edge_window = tk.OptionMenu(
            top, peak_detection_edge_var, 'Sigma', 'FWHM')
        peak_detection_edge_window.grid(row=5, column=1, sticky=tk.W)

        peak_detection_edge_value_label = tk.Label(
            top, text='Sigma Value', font=('Helvetica', 12))
        peak_detection_edge_value_label.grid(row=6, column=0, sticky=tk.W)
        peak_detection_edge_value_window = tk.Entry(top)
        peak_detection_edge_value_window.insert(
            0, master.peak_detection_edge_value)
        peak_detection_edge_value_window.grid(row=6, column=1, sticky=tk.W)

        # Calibration Settings
        calibration_label = tk.Label(top, text='Calibration Settings',
                                     font=('Helvetica', 16))
        calibration_label.grid(row=7, columnspan=2, sticky=tk.W)

        min_peak_label = tk.Label(top, text='Minimum Peaks',
                                  font=('Helvetica', 12))
        min_peak_label.grid(row=8, column=0, sticky=tk.W)
        min_peak_window = tk.Entry(top)
        min_peak_window.insert(0, master.min_peaks)
        min_peak_window.grid(row=8, column=1, sticky=tk.W)

        min_peak_SN_label = tk.Label(top, text='Minimum S/N',
                                     font=('Helvetica', 12))
        min_peak_SN_label.grid(row=9, column=0, sticky=tk.W)
        min_peak_SN_window = tk.Entry(top)
        min_peak_SN_window.insert(0, master.min_peak_SN)
        min_peak_SN_window.grid(row=9, column=1, sticky=tk.W)

        ultra_performance_calibration_label = tk.Label(
            top, text='UltraPerformanceCalibration', font=('Helvetica', 12))
        ultra_performance_calibration_label.grid(row=10, column=0, sticky=tk.W)
        ultra_performance_calibration_window = tk.OptionMenu(
            top, ultra_performance_calibration_var, *options)
        ultra_performance_calibration_window.grid(row=10, column=1,
                                                  sticky=tk.W)

        # Quantitation Settings
        quantitation_label = tk.Label(top, text='Quantitation Settings',
                                      font=('Helvetica', 16))
        quantitation_label.grid(row=11, columnspan=2, sticky=tk.W)

        points_label = tk.Label(top, text='Datapoints', font=('Helvetica', 12))
        points_label.grid(row=12, column=0, sticky=tk.W)
        points_window = tk.Entry(top)
        points_window.insert(0, master.points)
        points_window.grid(row=12, column=1, sticky=tk.W)

        baseline_order_label = tk.Label(top, text='Baseline Order',
                                        font=('Helvetica', 12))
        baseline_order_label.grid(row=13, column=0, sticky=tk.W)
        baseline_order_window = tk.Entry(top)
        baseline_order_window.insert(0, master.baseline_order)
        baseline_order_window.grid(row=13, column=1, sticky=tk.W)

        background_window_label = tk.Label(top, text='Background Window',
                                           font=('Helvetica', 12))
        background_window_label.grid(row=14, column=0, sticky=tk.W)
        background_window_window = tk.Entry(top)
        background_window_window.insert(0, master.background_window)
        background_window_window.grid(row=14, column=1, sticky=tk.W)

        slicepoints_label = tk.Label(top, text='MT Slice points',
                                     font=('Helvetica', 12))
        slicepoints_label.grid(row=15, column=0, sticky=tk.W)
        slicepoints_window = tk.Entry(top)
        slicepoints_window.insert(0, master.slicepoints)
        slicepoints_window.grid(row=15, column=1, sticky=tk.W)

        figure_label = tk.Label(top, text='Create figure for each analyte',
                                font=('Helvetica', 12))
        figure_label.grid(row=16, column=0, sticky=tk.W)
        figure_window = tk.OptionMenu(top, figure_variable, *options)
        figure_window.grid(row=16, column=1, sticky=tk.W)

        # Close/Save Buttons
        save_button = tk.Button(top, text='Save',
                                command=lambda: save(self, master))
        save_button.grid(row=17, column=0, sticky=tk.W)
        close_button = tk.Button(top, text='Close',
                                 command=lambda: close(self, master))
        close_button.grid(row=17, column=1, sticky=tk.E)

        # Tooltips
        create_tooltip(
            points_label, 'The number of data points that is used to ' +
            'determine the baseline. Specifically, this setting specifies ' +
            'how large each segment of the whole chromatogram will be to ' +
            'identify the lowest data point per window, i.e. a setting of ' +
            '100 means that the chromatogram is split into segments of 100 ' +
            'data points per segment.')

        create_tooltip(
            start_label, 'This setting tells the program from which ' +
            'time point it is supposed to begin processing, this setting ' +
            'should be set in such a way that it is before the analytes of ' +
            'interest but after any potential big increases or decrease in ' +
            'intensity.')

        create_tooltip(
            end_label, 'This setting tells the program until which ' +
            'time point it is supposed to begin processing, this setting ' +
            'should be set in such a way that it is after the analytes of ' +
            'interest but before any potential big increases or decrease ' +
            'in intensity.')

        create_tooltip(
            baseline_order_label, 'This setting tells the program ' +
            'what sort of function should be used to correct the baseline. ' +
            'A value of 1 refers to a linear function, while a value of 2 ' +
            'refers to a quadratic function. We advise to use a linear ' +
            'function as the need for any higher order function indicates ' +
            'an unexpected event in the chromatography.')

        create_tooltip(
            background_window_label, 'This setting tells the program ' +
            'the size of the region that will be examined to determine the ' +
            'background. A value of 1 means that the program will look ' +
            'from 20.0 to 21.0 minutes and 21.4 to 22.4 for an analyte that ' +
            'elutes from 21.0 to 21.4 minutes.')

        create_tooltip(
            slicepoints_label, 'The number of conscutive data points ' +
            'that will be used to determine the background and noise using ' +
            'the MT method. The MT method will scan all datapoints that ' +
            'fall within the background window (specified above) to find ' +
            'the here specified number of consecutive data points that ' +
            'yield the lowest average intensity, the average of these data ' +
            'points is then used as background while the standard deviation ' +
            'of these data points is used as the noise.')

        create_tooltip(
            figure_label, 'This setting specifies if HappyTools ' +
            'should create a figure for each integrated peak, showing the ' +
            'raw datapoints, background, noise, S/N and GPQ values. This is ' +
            'a very performance intensive option and it is recommended to ' +
            'only use this on a subset of your samples (e.g. less than 25 ' +
            'samples).')

        create_tooltip(
            min_peak_label, 'This setting specifies the minimum ' +
            'number of calibrant peaks that have to pass the specified S/N ' +
            'value that must be present in a chromatogram. A chromatogram ' +
            'for which there are not enough calibrant peaks passing the ' +
            'specified criteria will not be calibrated and excluded from ' +
            'further quantitation.')

        create_tooltip(
            min_peak_SN_label, 'This setting specifies the minimum ' +
            'S/N value a calibrant peak must surpass to be included in the ' +
            'calibration. The actual S/N value that is determined by ' +
            'HappyTools depends heavily on which method to determine signal ' +
            'and noise is used, the default method being rather conservative.')

        create_tooltip(
            peak_detection_label, 'This setting specifies the minimum ' +
            'intensity, relative to the main peak in a chromatogram, that ' +
            'the peak detection algorithm will try to annotate. For ' +
            'example, a value of 0.01 means that the program will attempt ' +
            'to annotate peaks until the next highest peak is below 1% of ' +
            'the intensity of the main peak in the chromatogram.')

        create_tooltip(
            peak_detection_edge_label, 'This setting specifies if ' +
            'HappyTools will determine the integration window using either ' +
            'the full width at half maximum (FWHM) or a specified sigma ' +
            'value. The Sigma value has to be specified in the field below ' +
            'if Sigma is the selected method.')

        create_tooltip(
            peak_detection_edge_value_label, 'This setting specifies ' +
            'the Sigma value that will be used for determining the border ' +
            'of the integration window. A value of 1.0 means that ' +
            'HappyTools will set the integration window so that 68.3% of ' +
            'the Gaussian peak will be quantified (2 Sigma = 95.5% and 3 ' +
            'sigma = 99.7%). Please note that this value should depend on ' +
            'how complex the chromatogram is, for instance a low sigma will ' +
            'yield better results in a complex chromatogram.')

        create_tooltip(
            ultra_performance_calibration_label, 'This setting ' +
            'specifies if HappyTools will use the standard calibration ' +
            'method using a second degree polynomial (False) or if ' +
            'HappyTools is allowed to use the experimental ' +
            'ultraPerformanceCalibration (True) that attempts to identify ' +
            'the optimal function for calibration based on the residual ' +
            'time error of the calibrants.\n\nNOTE: The ' +
            'ultraPerformanceCalibration is not ready to be used in a ' +
            'production environment yet!')
