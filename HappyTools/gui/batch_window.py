from HappyTools.gui.output_window import OutputWindow
from HappyTools.util.functions import create_tooltip
import tkinter as tk
import tkinter.filedialog as filedialog


class batchWindow(object):
    def __init__(self, master):
        """Create a batch processing pop-up.

        This function creates a new tkinter window that is used to
        control the batch processing. Specifically, it allows the user
        to select a calibration file, an analyte file, select the
        desired outputs (by calling the outputPopup function) and
        starting the batch process.

        Keyword arguments:
        none
        """
        self.master = master
        self.functions = master.functions
        self.settings = master.settings

        self.output_window_open = tk.IntVar(value=0)
        self.abs_int = tk.IntVar(value=0)
        self.rel_int = tk.IntVar(value=0)
        self.gauss_int = tk.IntVar(value=0)
        self.bck_sub = tk.IntVar(value=0)
        self.bck_noise = tk.IntVar(value=0)
        self.peak_qual = tk.IntVar(value=0)

        self.cal_file = tk.StringVar()
        self.anal_file = tk.StringVar()
        self.batch_folder = tk.StringVar()


        top = tk.top = tk.Toplevel()
        top.title('Batch Process')
        top.protocol('WM_DELETE_WINDOW', self.close)

        calibrationButton = tk.Button(top, text='Calibration File',
                                      width=20,
                                      command=self.set_calibration_file)
        calibrationButton.grid(row=1, column=0, sticky=tk.W)
        calibrationLabel = tk.Label(tk.top, textvariable=self.cal_file,
                                    width=20)
        calibrationLabel.grid(row=1, column=1)

        analyteButton = tk.Button(top, text='Analyte File',
                                  width=20,
                                  command=self.set_analyte_file)
        analyteButton.grid(row=2, column=0, sticky=tk.W)
        analyteLabel = tk.Label(tk.top, textvariable=self.anal_file,
                                width=20)
        analyteLabel.grid(row=2, column=1)

        batchButton = tk.Button(tk.top, text='Batch Directory',
                                width=20, command=self.set_batch_folder)
        batchButton.grid(row=3, column=0, sticky=tk.W)
        batchLabel = tk.Label(tk.top, textvariable=self.batch_folder,
                              width=20)
        batchLabel.grid(row=3, column=1, sticky=tk.W)

        outputButton = tk.Button(tk.top, text='Output Options',
                                 command=self.open_output_window)
        outputButton.grid(row=4, column=0, columnspan=2,
                          sticky=tk.E+tk.W)

        runButton = tk.Button(tk.top, text='Run', width=20,
                              command=self.run)
        runButton.grid(row=5, column=0, sticky=tk.W)
        closeButton = tk.Button(tk.top, text='Close', width=20,
                                command=self.close)
        closeButton.grid(row=5, column=1, sticky=tk.E)

        top.lift()
        self.top = top

        # Tooltips
        create_tooltip(
            calibrationButton, 'This button will allow you to select ' +
            'your calibration file, the program expects a tab separated ' +
            'text file where each line consists of a peak ID, peak RT and ' +
            'a RT window.')

        create_tooltip(
            analyteButton, 'This button will allow you to select your ' +
            'analyte file, the program expects a tab separated text file ' +
            'where each line consists of a peak ID, peak RT and a RT window.')

        create_tooltip(
            batchButton, 'This button will allow you to select the ' +
            'folder where the chromatograms are stored that HappyTools will ' +
            'process.')

        create_tooltip(
            outputButton, 'This button will open another window in ' +
            'which you can select which outputs you want HappyTools to show ' +
            'in the final summary.')

    def close(self):
        """Close the batch processing pop-up.
        """
        self.top.destroy()

    def open_output_window(self):
        OutputWindow(self)

    def run(self):
        """Start the batch process.
        """
        self.functions.batch_process(self)

    def set_batch_folder(self):
        """Ask for the batch folder.
        """
        self.batch_folder.set(filedialog.askdirectory(
            title='Batch Folder'))

    def set_calibration_file(self):
        """Ask for the calibration file.
        """
        self.cal_file.set(filedialog.askopenfilename(
            title='Calibration File'))

    def set_analyte_file(self):
        """Ask for the analyte file.
        """
        self.anal_file.set(filedialog.askopenfilename(
            title='Analyte File'))
