from HappyTools.gui.output_window import OutputWindow
import HappyTools.gui.version as version
try:
    # Python 2
    import Tkinter as tk
    import tkFileDialog as filedialog
except ImportError:
    # Python 3
    import tkinter as tk
    import tk.filedialog as filedialog

class batchWindow(object):
    def __init__(self, master):
        """Create a batch processing pop-up.

        This function creates a new tkinter window that is used to control
        the batch processing. Specifically, it allows the user to select a
        calibration file, an analyte file, select the desired outputs (by
        calling the outputPopup function) and starting the batch process.

        Keyword arguments:
        none
        """
        self.master = master
        self.functions = master.functions
        self.settings = master.settings

        self.output_window = tk.IntVar()
        self.abs_int = tk.IntVar()
        self.rel_int = tk.IntVar()
        self.gauss_int = tk.IntVar()
        self.bck_sub = tk.IntVar()
        self.bck_noise = tk.IntVar()
        self.peak_qual = tk.IntVar()

        self.cal_file = tk.StringVar()
        self.anal_file = tk.StringVar()
        self.batch_folder = tk.StringVar()

        def close():
            """Close the batch processing pop-up.
            """
            top.destroy()

        def set_calibration_file():
            """Ask for the calibration file.
            """
            self.cal_file.set(filedialog.askopenfilename(title="Calibration File"))

        def set_analyte_file():
            """Ask for the analyte file.
            """
            self.anal_file.set(filedialog.askopenfilename(title="Analyte File"))

        def set_batch_folder():
            """Ask for the batch folder.
            """
            self.batch_folder.set(filedialog.askdirectory(title="Batch Folder"))

        def run():
            """Start the batch process.
            """
            self.functions.batch_process(self)

        top = tk.top = tk.Toplevel()
        top.title("Batch Process")
        top.protocol("WM_DELETE_WINDOW", lambda: close())

        calibrationButton = tk.Button(top, text="Calibration File", width=20, command=lambda: set_calibration_file())
        calibrationButton.grid(row=1, column=0, sticky=tk.W)
        calibrationLabel = tk.Label(tk.top, textvariable=self.cal_file, width=20)
        calibrationLabel.grid(row=1, column=1)

        analyteButton = tk.Button(top, text="Analyte File", width=20, command=lambda: set_analyte_file())
        analyteButton.grid(row=2, column=0, sticky=tk.W)
        analyteLabel = tk.Label(tk.top, textvariable=self.anal_file, width=20)
        analyteLabel.grid(row=2, column=1)

        batchButton = tk.Button(tk.top, text="Batch Directory", width=20, command=lambda: set_batch_folder())
        batchButton.grid(row=3, column=0, sticky=tk.W)
        batchLabel = tk.Label(tk.top, textvariable=self.batch_folder, width=20)
        batchLabel.grid(row=3, column=1, sticky=tk.W)
        
        outputButton = tk.Button(tk.top, text="Output Options", command=lambda: OutputWindow(self))
        outputButton.grid(row=4, column=0, columnspan=2, sticky=tk.E+tk.W)

        runButton = tk.Button(tk.top, text="Run", width=20, command=lambda: run())
        runButton.grid(row=5, column=0, sticky=tk.W)
        closeButton = tk.Button(tk.top, text="Close", width=20, command=lambda: close())
        closeButton.grid(row=5, column=1, sticky=tk.E)

        #self.master.lift()
        # Tooltips
        #createToolTip(calibrationButton,"This button will allow you to select your calibration file, the program expects a "+
                #"tab separated text file where each line consists of a peak ID, peak RT and a RT window.")
        #createToolTip(analyteButton,"This button will allow you to select your analyte file, the program expects a tab separated "+
                #"text file where each line consists of a peak ID, peak RT and a RT window.")
        #createToolTip(outputButton,"This button will open another window in which you can select which outputs you want "+
                #"HappyTools to show in the final summary.")
