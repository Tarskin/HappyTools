from os import path, getcwd
try:
    # Python 2
    import Tkinter as tk
except ImportError:
    # Python 3
    import tkinter as tk

class Settings(object):
    def __init__(self, master):

        self.points = 100
        self.start = 10
        self.end = 60
        self.baseline_order = 1
        self.background_window = 1
        self.noban_start = 0.25
        self.slicepoints = 5
        self.peak_detection_min = 0.05
        self.peak_detection_edge = "Sigma"
        self.peak_detection_edge_value = 2.0
        self.create_figure = "True"
        self.min_peaks = 4
        self.min_peak_SN = 27
        self.use_UPC = "False"

        self.decimal_numbers = 6
        self.min_improvement = 0.05
        self.use_interpolation = False
        self.noise = "RMS"
        self.background_noise_method = "MT"

        self.output = "summary.results"
        self.settings = "HappyTools.ini"

        self.exclusion_files = ["LICENSE.txt","CHANGELOG.txt"]
        self.calibration_filetypes = ["*.txt","*.arw"]
        self.integration_filetypes = ["calibrated*.txt"]

    def read_settings(self, master):
        """Read the settings file.

        This function opens the settings file (default is HappyTools.ini),
        parses the lines of the settings file and takes the value from the
        settings file as a value for the changeable settings (e.g. the start
        variable can be read from the settings file).

        Keyword arguments:
        none
        """
        with open(path.join(getcwd(), master.settings),'r') as fr:
            for line in fr:
                chunks = line.strip('\n').split('\t')
                if chunks[0] == "points:":
                    master.points = int(chunks[1])
                elif chunks[0] == "start:":
                    master.start = float(chunks[1])
                elif chunks[0] == "end:":
                    master.end = float(chunks[1])
                elif chunks[0] == "baselineOrder:":
                    master.baseline_order = int(chunks[1])
                elif chunks[0] == "backgroundWindow:":
                    master.background_window = float(chunks[1])
                elif chunks[0] == "noise:":
                    master.noise = str(chunks[1])
                elif chunks[0] == "slicepoints:":
                    master.slicepoints = int(chunks[1])
                elif chunks[0] == "createFigure:":
                    master.create_figure = str(chunks[1])
                elif chunks[0] == "minPeaks:":
                    master.min_peaks = int(chunks[1])
                elif chunks[0] == "minPeakSN:":
                    master.min_peak_SN = int(chunks[1])
                elif chunks[0] == "peakDetectionMin:":
                    master.peak_detection_min = float(chunks[1])
                elif chunks[0] == "peakDetectionEdge:":
                    master.peak_detection_edge = str(chunks[1])
                elif chunks[0] == "peakDetectionEdgeValue:":
                    master.peak_detection_edge_value = float(chunks[1])

    def settings_popup(self, master):
        """ TODO: Redesign settings window (it's f'ing ugly)
        """

        self.figure_variable = tk.StringVar()
        self.figure_variable.set(master.create_figure)
        self.peak_detection_edge_var = tk.StringVar()
        self.peak_detection_edge_var.set(master.peak_detection_edge)
        self.ultra_performance_calibration_var = tk.StringVar()
        self.ultra_performance_calibration_var.set(master.use_UPC)
        options = ["True", "False"]

        def close(self, master):
            """ TODO
            """
            master.points = int(self.points_window.get())
            master.start = float(self.start_window.get())
            master.end = float(self.end_window.get())
            master.baseline_order = int(self.baseline_order_window.get())
            master.background_window = float(self.background_window_window.get())
            master.slicepoints = int(self.slicepoints_window.get())
            master.create_figure = str(self.figure_variable.get())
            master.min_peaks = int(self.min_peak_window.get())
            master.min_peak_SN = int(self.min_peak_SN_window.get())
            master.peak_detection_min = float(self.peak_detection.get())
            master.peak_detection_edge = str(self.peak_detection_edge_var.get())
            master.peak_detection_edge_value = float(self.peak_detection_edge_value_window.get())
            master.use_UPC = str(self.ultra_performance_calibration_var.get())
            top.destroy()
        
        def save(self, master):
            """ TODO
            """
            with open(path.join(getcwd(), master.settings),'w') as fw:
                fw.write("points:\t"+str(int(self.points_window.get()))+"\n")
                fw.write("start:\t"+str(float(self.start_window.get()))+"\n")
                fw.write("end:\t"+str(float(self.end_window.get()))+"\n")
                fw.write("baselineOrder:\t"+str(int(self.baseline_order_window.get()))+"\n")
                fw.write("backgroundWindow:\t"+str(float(self.background_window_window.get()))+"\n")
                fw.write("noise:\t"+str(self.noise)+"\n")
                fw.write("slicepoints:\t"+str(self.slicepoints)+"\n")
                fw.write("createFigure:\t"+str(self.figure_variable.get())+"\n")
                fw.write("minPeaks:\t"+str(self.min_peak_window.get())+"\n")
                fw.write("minPeakSN:\t"+str(self.min_peak_SN_window.get())+"\n")
                fw.write("peakDetectionMin:\t"+str(self.peak_detection.get())+"\n")
                fw.write("peakDetectionEdge:\t"+str(self.peak_detection_edge_var.get())+"\n")
                fw.write("peakDetectionEdgeValue:\t"+str(self.peak_detection_edge_value_window.get())+"\n")
                fw.write("useUPC:\t"+str(self.ultra_performance_calibration_var.get())+"\n")

        top = tk.top = tk.Toplevel()
        top.title("Settings")
        top.protocol( "WM_DELETE_WINDOW", lambda: close(self, master))

        # General Settings
        self.general_label = tk.Label(top, text="General Settings", font=("Helvetica", 16))
        self.general_label.grid(row=0, columnspan=2, sticky=tk.W)
        
        self.start_label = tk.Label(top, text="Start Time", font=("Helvetica", 12))
        self.start_label.grid(row=1, column=0, sticky=tk.W)
        self.start_window = tk.Entry(top)
        self.start_window.insert(0, master.start)
        self.start_window.grid(row=1, column=1, sticky=tk.W)
        
        self.end_label = tk.Label(top, text="End Time", font=("Helvetica", 12))
        self.end_label.grid(row=2, column=0, sticky=tk.W)
        self.end_window = tk.Entry(top)
        self.end_window.insert(0, master.end)
        self.end_window.grid(row=2, column=1, sticky=tk.W)
        
        # Peak Detection Settings
        self.peak_detection_label = tk.Label(top, text="Peak Detection Settings", font=("Helvetica", 16))
        self.peak_detection_label.grid(row=3, columnspan=2, sticky=tk.W)

        self.peak_detection_label = tk.Label(top, text="Minimum Intensity", font=("Helvetica", 12))
        self.peak_detection_label.grid(row=4, column=0, sticky=tk.W)
        self.peak_detection = tk.Entry(top)
        self.peak_detection.insert(0, master.peak_detection_min)
        self.peak_detection.grid(row=4, column=1, sticky=tk.W)

        self.peak_detection_edge_label = tk.Label(top, text="Edge Method", font=("Helvetica", 12))
        self.peak_detection_edge_label.grid(row=5, column=0, sticky=tk.W)
        self.peak_detection_edge_window = tk.OptionMenu(top, master.peak_detection_edge_var, "Sigma", "FWHM")
        self.peak_detection_edge_window.grid(row=5, column=1, sticky=tk.W)

        self.peak_detection_edge_value_label = tk.Label(top, text="Sigma Value", font=("Helvetica", 12))
        self.peak_detection_edge_value_label.grid(row=6, column=0, sticky=tk.W)
        self.peak_detection_edge_value_window = tk.Entry(top)
        self.peak_detection_edge_value_window.insert(0, master.peak_detection_edge_value)
        self.peak_detection_edge_value_window.grid(row=6, column=1, sticky=tk.W)

        # Calibration Settings
        self.calibration_label = tk.Label(top, text="Calibration Settings", font=("Helvetica", 16))
        self.calibration_label.grid(row=7, columnspan=2, sticky=tk.W)

        self.min_peak_label = tk.Label(top, text="Minimum Peaks", font=("Helvetica", 12))
        self.min_peak_label.grid(row=8, column=0, sticky=tk.W)
        self.min_peak_window = tk.Entry(top)
        self.min_peak_window.insert(0, master.min_peaks)
        self.min_peak_window.grid(row=8, column=1, sticky=tk.W)

        self.min_peak_SN_label = tk.Label(top, text="Minimum S/N", font=("Helvetica", 12))
        self.min_peak_SN_label.grid(row=9, column=0, sticky=tk.W)
        self.min_peak_SN_window = tk.Entry(top)
        self.min_peak_SN_window.insert(0, master.min_peak_SN)
        self.min_peak_SN_window.grid(row=9, column=1, sticky=tk.W)
        
        self.ultra_performance_calibration_label = tk.Label(top, text="UltraPerformanceCalibration", font=("Helvetica", 12))
        self.ultra_performance_calibration_label.grid(row=10, column=0, sticky=tk.W)
        self.ultra_performance_calibration_window = tk.OptionMenu(top, master.ultra_performance_calibration_var, *options)
        self.ultra_performance_calibration_window.grid(row=10, column=1, sticky=tk.W)
        
        # Quantitation Settings
        self.quantitation_label = tk.Label(top, text="Quantitation Settings", font=("Helvetica", 16))
        self.quantitation_label.grid(row=11, columnspan=2, sticky=tk.W)
        
        self.points_label = tk.Label(top, text="Datapoints", font=("Helvetica", 12))
        self.points_label.grid(row=12, column=0, sticky=tk.W)
        self.points_window = tk.Entry(top)
        self.points_window.insert(0, master.points)
        self.points_window.grid(row=12, column=1, sticky=tk.W)

        self.baseline_order_label = tk.Label(top, text="Baseline Order", font=("Helvetica", 12))
        self.baseline_order_label.grid(row=13, column=0, sticky=tk.W)
        self.baseline_order_window = tk.Entry(top)
        self.baseline_order_window.insert(0, master.baseline_order)
        self.baseline_order_window.grid(row=13, column=1, sticky=tk.W)
        
        self.background_window_label = tk.Label(top, text="Background Window", font=("Helvetica", 12))
        self.background_window_label.grid(row=14, column=0, sticky=tk.W)
        self.background_window_window = tk.Entry(top)
        self.background_window_window.insert(0, master.background_window)
        self.background_window_window.grid(row=14, column=1, sticky=tk.W)

        self.slicepoints_label = tk.Label(top, text="MT Slice points", font=("Helvetica", 12))
        self.slicepoints_label.grid(row=15, column=0, sticky=tk.W)
        self.slicepoints_window = tk.Entry(top)
        self.slicepoints_window.insert(0, master.slicepoints)
        self.slicepoints_window.grid(row=15, column=1, sticky=tk.W)

        self.figure_label = tk.Label(top, text="Create figure for each analyte", font=("Helvetica", 12))
        self.figure_label.grid(row=16, column=0, sticky=tk.W)
        self.figure_window = tk.OptionMenu(top, master.figure_variable, *options)
        self.figure_window.grid(row=16, column=1, sticky=tk.W)

        # Close/Save Buttons
        self.save_button = tk.Button(top, text="Save", command=lambda:save(self, master))
        self.save_button.grid(row=17, column=0, sticky=tk.W)
        self.close_button = tk.Button(top, text="Close", command=lambda:close(self, master))
        self.close_button.grid(row=17, column=1, sticky=tk.E)
