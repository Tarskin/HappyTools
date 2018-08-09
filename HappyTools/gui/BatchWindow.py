from __future__ import absolute_import
from HappyTools.util.Functions import Functions

import HappyTools.gui.Settings as settings
import HappyTools.gui.Version as version
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

        self.output_window = tk.IntVar()
        self.abs_int = tk.IntVar()
        self.rel_int = tk.IntVar()
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
            Functions().batch_process(self.cal_file, self.anal_file, self.batch_folder)

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
        
        outputButton = tk.Button(tk.top, text="Output Options", command=lambda: self.output_popup(self))
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

    def output_popup(self, master):
        """Create a pop-up enabling output selection.

        This function creates a pop up box that allows the user to specify 
        what output should be shown in the final summary. The default value 
        for all variables is off (0) and by ticking a box it is set to on
        (1).

        Keyword arguments:
        none
        """
        if master.output_window.get() == 1:
            return
        master.output_window.set(1)

        def select_all():
            """Set all variables to on (1).
            """
            master.abs_int.set(1)
            master.rel_int.set(1)
            master.bck_sub.set(1)
            master.bck_noise.set(1)
            master.peak_qual.set(1)

        def select_none():
            """Set all variables to off (0).
            """
            master.abs_int.set(0)
            master.rel_int.set(0)
            master.bck_sub.set(0)
            master.bck_noise.set(0)
            master.peak_qual.set(0)

        def close():
            """Close the output pop-up.
            """
            master.output_window.set(0)
            top.destroy()

        top = tk.Toplevel()
        top.protocol("WM_DELETE_WINDOW", lambda: close())
        top.title("Output Options")
        selAll = tk.Button(top, text="Select All", command=lambda: select_all())
        selAll.grid(row=0, column=0, sticky=tk.W)
        none = tk.Button(top, text="Select None", command=lambda: select_none())
        none.grid(row=0, column=1, sticky=tk.E)
        text1 = tk.Label(top, text="Base Outputs", font="bold")
        text1.grid(row=1, column=0, sticky=tk.W)
        text2 = tk.Label(top, text="Output Modifiers", font="bold")
        text2.grid(row=1, column=1, sticky=tk.W)
        ai = tk.Checkbutton(top, text=u"Analyte Intensity\u00B9", variable=master.abs_int, onvalue=1, offvalue=0)
        ai.grid(row=2, column=0, sticky=tk.W)
        ri = tk.Checkbutton(top, text=u"Relative Intensity\u00B9", variable=master.rel_int, onvalue=1, offvalue=0)
        ri.grid(row=3, column=0, sticky=tk.W)
        pq = tk.Checkbutton(top, text="Peak Quality Criteria", variable=master.peak_qual, onvalue=1, offvalue=0)
        pq.grid(row=4, column=0, sticky=tk.W)
        bn = tk.Checkbutton(top, text="Background and Noise", variable=master.bck_noise, onvalue=1, offvalue=0)
        bn.grid(row=5, column=0, sticky=tk.W)
        bck = tk.Checkbutton(top, text=u"\u00B9Background subtracted Intensities", variable=master.bck_sub, onvalue=1, offvalue=0)
        bck.grid(row=2, column=1, sticky=tk.W)
        button = tk.Button(top,text='Ok',command = lambda: close())
        button.grid(row = 6, column = 0, columnspan = 2)
        top.lift()
        return

    def close(self):
        self.master.destroy()
