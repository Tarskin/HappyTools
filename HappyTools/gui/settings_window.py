from HappyTools.gui.tooltip import create_tooltip
import tkinter as tk


class SettingsWindow(object):
    def __init__(self, master):
        self.settings = master.settings
        self.create_window()

    def create_window(self):
        root = tk.Toplevel()
        root.protocol('WM_DELETE_WINDOW', self.close_settings_window)

        # Buffer Variables
        peak_edge_type_buffer = tk.StringVar()
        peak_edge_type_buffer.set(self.settings.peak_detection_edge)

        # General Settings
        general_label = tk.Label(root, text='General Settings',
                                 font='bold')
        general_label.grid(row=0, columnspan=2, sticky=tk.W)

        start_time_label = tk.Label(root, text='Start Time')
        start_time_label.grid(row=1, column=0, sticky=tk.W)
        start_time = tk.Entry(root)
        start_time.insert(0, self.settings.start)
        start_time.grid(row=1, column=1, sticky=tk.W)

        end_time_label = tk.Label(root, text='End Time')
        end_time_label.grid(row=2, column=0, sticky=tk.W)
        end_time = tk.Entry(root)
        end_time.insert(0, self.settings.end)
        end_time.grid(row=2, column=1, sticky=tk.W)

        baseline_datapoints_label = tk.Label(root, text='Number of '+
                                             'Datapoints')
        baseline_datapoints_label.grid(row=3, column=0, sticky=tk.W)
        baseline_datapoints = tk.Entry(root)
        baseline_datapoints.insert(0, self.settings.points)
        baseline_datapoints.grid(row=3, column=1, sticky=tk.W)

        baseline_order_label = tk.Label(root, text='Baseline Function '+
                                        'order')
        baseline_order_label.grid(row=4, column=0, sticky=tk.W)
        baseline_order = tk.Entry(root)
        baseline_order.insert(0, self.settings.baseline_order)
        baseline_order.grid(row=4, column=1, sticky=tk.W)

        background_window_label = tk.Label(root, text='Background '+
                                           'Window')
        background_window_label.grid(row=5, column=0, sticky=tk.W)
        background_window = tk.Entry(root)
        background_window.insert(0, self.settings.background_window)
        background_window.grid(row=5, column=1, sticky=tk.W)

        # Peak Detection Settings
        peak_detection_label = tk.Label(root, text='Peak Detection '+
                                        'Settings', font='bold')
        peak_detection_label.grid(row=6, columnspan=2, sticky=tk.W)

        minimum_peak_intensity_label = tk.Label(root, text='Minimum '+
                                                'Peak Intensity')
        minimum_peak_intensity_label.grid(row=7, column=0, sticky=tk.W)
        minimum_peak_intensity = tk.Entry(root)
        minimum_peak_intensity.insert(
                0, self.settings.peak_detection_min)
        minimum_peak_intensity.grid(row=7, column=1, sticky=tk.W)

        peak_edge_type_label = tk.Label(root, text='Peak Edge Type')
        peak_edge_type_label.grid(row=8, column=0, sticky=tk.W)
        peak_edge_type = tk.OptionMenu(
            root, peak_edge_type_buffer, 'Sigma', 'FWHM')
        peak_edge_type.grid(row=8, column=1, sticky=tk.W)

        peak_edge_value_label = tk.Label(root, text='Peak Edge Value')
        peak_edge_value_label.grid(row=9, column=0, sticky=tk.W)
        peak_edge_value = tk.Entry(root)
        peak_edge_value.insert(
                0, self.settings.peak_detection_edge_value)
        peak_edge_value.grid(row=9, column=1, sticky=tk.W)

        # Calibration Settings
        calibration_label = tk.Label(root, text='Calibration Settings',
                                     font='bold')
        calibration_label.grid(row=10, columnspan=2, sticky=tk.W)

        minimum_number_calibrants_label = tk.Label(
                root, text='Minimum Number of Calibrants')
        minimum_number_calibrants_label.grid(row=11, column=0,
                                             sticky=tk.W)
        minimum_number_calibrants = tk.Entry(root)
        minimum_number_calibrants.insert(0, self.settings.min_peaks)
        minimum_number_calibrants.grid(row=11, column=1, sticky=tk.W)

        minimum_calibrant_sn_label = tk.Label(root, text='Minimum '+
                                              'Calibrant S/N')
        minimum_calibrant_sn_label.grid(row=12, column=0, sticky=tk.W)
        minimum_calibrant_sn = tk.Entry(root)
        minimum_calibrant_sn.insert(0, self.settings.min_peak_SN)
        minimum_calibrant_sn.grid(row=12, column=1, sticky=tk.W)

        # Close/Save Buttons
        save_button = tk.Button(root, text='Save',
                                command=self.save_settings)
        save_button.grid(row=17, column=0, sticky=tk.W)
        close_button = tk.Button(root, text='Close',
                                 command=self.close_settings_window)
        close_button.grid(row=17, column=1, sticky=tk.E)

        # Inheritance
        self.root = root
        self.start_time = start_time
        self.end_time = end_time
        self.baseline_datapoints = baseline_datapoints
        self.baseline_order = baseline_order
        self.background_window = background_window
        self.minimum_peak_intensity = minimum_peak_intensity
        self.peak_edge_type = peak_edge_type_buffer
        self.peak_edge_value = peak_edge_value
        self.minimum_number_calibrants = minimum_number_calibrants
        self.minimum_calibrant_sn = minimum_calibrant_sn

        # General Settings Tooltips
        create_tooltip(
            start_time_label, 'This setting tells the program from '+
            'which time point it is supposed to begin processing, '+
            'this setting should be set in such a way that it is '+
            'before the analytes of interest but after any potential '+
            'big increases or decrease in intensity.')

        create_tooltip(
            end_time_label, 'This setting tells the program until '+
            'which time point it is supposed to begin processing, '+
            'this setting should be set in such a way that it is '+
            'after the analytes of interest but before any potential '+
            'big increases or decrease in intensity.')

        create_tooltip(
            baseline_datapoints_label, 'The number of data points '+
            'that is used to determine the baseline. Specifically, '+
            'this setting specifies how large each segment of the '+
            'whole chromatogram will be to identify the lowest data '+
            'point per window, i.e. a setting of 100 means that the '+
            'chromatogram is split into segments of 100 data points '+
            'per segment.')

        create_tooltip(
            baseline_order_label, 'This setting tells the program '+
            'what sort of function should be used to correct the '+
            'baseline. A value of 1 refers to a linear function, '+
            'while a value of 2 refers to a quadratic function. We '+
            'advise to use a linear function as the need for any '+
            'higher order function indicates an unexpected event in '+
            'the chromatography.')

        create_tooltip(
            background_window_label, 'This setting tells the program '+
            'the size of the region that will be examined to '+
            'determine the background. A value of 1 means that the '+
            'program will look from 20.0 to 21.0 minutes and 21.4 to '+
            '22.4 for an analyte that elutes from 21.0 to 21.4 '+
            'minutes.')

        # Peak Detection Settings Tooltips
        create_tooltip(
            minimum_peak_intensity_label, 'This setting specifies the '+
            'minimum intensity, relative to the main peak in a '+
            'chromatogram, that the peak detection algorithm will try '+
            'to annotate. For example, a value of 0.01 means that the '+
            'program will attempt to annotate peaks until the next '+
            'highest peak is below 1% of the intensity of the main '+
            'peak in the chromatogram.')

        create_tooltip(
            peak_edge_type_label, 'This setting specifies if ' +
            'HappyTools will determine the integration window using '+
            'either the full width at half maximum (FWHM) or a '+
            'specified sigma value. The Sigma value has to be '+
            'specified in the field below if Sigma is the selected '+
            'method.')

        create_tooltip(
            peak_edge_value_label, 'This setting specifies the Sigma '+
            'value that will be used for determining the border of '+
            'the integration window. A value of 1.0 means that '+
            'HappyTools will set the integration window so that 68.3% '+
            'of the Gaussian peak will be quantified (2 Sigma = 95.5% '+
            'and 3 sigma = 99.7%). Please note that this value should '+
            'depend on how complex the chromatogram is, for instance '+
            'a low sigma will yield better results in a complex '+
            'chromatogram.')

        # Calibration Settings Tooltips
        create_tooltip(
            minimum_number_calibrants_label, 'This setting specifies '+
            'the minimum number of calibrant peaks that have to pass '+
            'the specified S/N value that must be present in a '+
            'chromatogram. A chromatogram for which there are not '+
            'enough calibrant peaks passing the specified criteria '+
            'will not be calibrated and excluded from further '+
            'quantitation.')

        create_tooltip(
            minimum_calibrant_sn_label, 'This setting specifies the '+
            'minimum S/N value a calibrant peak must surpass to be '+
            'included in the calibration. The actual S/N value that '+
            'is determined by HappyTools depends heavily on which '+
            'method to determine signal and noise is used, the '+
            'default method being rather conservative.')

    # Rewrite
    def save_settings(self):
        self.store_settings()
        self.settings.save_to_disk()

    def store_settings(self):
        # General Settings
        self.settings.points = int(
            self.baseline_datapoints.get()
        )

        self.settings.baseline_order = int(
            self.baseline_order.get()
        )

        self.settings.start = float(
            self.start_time.get()
        )

        self.settings.end = float(
            self.end_time.get()
        )

        self.settings.background_window = float(
            self.background_window.get()
        )

        # Peak Detection Settings
        self.settings.peak_detection_min = float(
            self.minimum_peak_intensity.get()
        )

        self.settings.peak_detection_edge = str(
            self.peak_edge_type.get()
        )

        self.settings.peak_detection_edge_value = float(
            self.peak_edge_value.get()
        )

        # Calibration Settings
        self.settings.min_peaks = int(
            self.minimum_number_calibrants.get()
        )

        self.settings.min_peak_SN = float(
            self.minimum_calibrant_sn.get()
        )

    def close_settings_window(self):
        self.store_settings()
        self.root.destroy()
