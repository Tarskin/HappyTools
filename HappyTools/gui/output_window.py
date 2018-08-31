import tkinter as tk


class OutputWindow(object):
    def __init__(self, master):
        """Create a pop-up enabling output selection.

        This function creates a pop up box that allows the user to specify
        what output should be shown in the final summary. The default value
        for all variables is off (0) and by ticking a box it is set to on
        (1).

        Keyword arguments:
        none
        """
        if master.output_window_open.get() == 1:
            return
        master.output_window_open.set(1)

        top = tk.Toplevel()
        top.protocol("WM_DELETE_WINDOW", self.close)
        top.title("Output Options")
        selAll = tk.Button(top, text="Select All",
                           command=self.select_all)
        selAll.grid(row=0, column=0, sticky=tk.W)
        none = tk.Button(top, text="Select None",
                         command=self.select_none)
        none.grid(row=0, column=1, sticky=tk.E)
        text1 = tk.Label(top, text="Base Outputs", font="bold")
        text1.grid(row=1, column=0, sticky=tk.W)
        text2 = tk.Label(top, text="Output Modifiers", font="bold")
        text2.grid(row=1, column=1, sticky=tk.W)
        ai = tk.Checkbutton(top, text=u"Analyte Intensity\u00B9",
                            variable=master.abs_int, onvalue=1, offvalue=0)
        ai.grid(row=2, column=0, sticky=tk.W)
        ri = tk.Checkbutton(top, text=u"Relative Intensity\u00B9",
                            variable=master.rel_int, onvalue=1, offvalue=0)
        ri.grid(row=3, column=0, sticky=tk.W)
        pq = tk.Checkbutton(top, text="Peak Quality Criteria",
                            variable=master.peak_qual, onvalue=1, offvalue=0)
        pq.grid(row=4, column=0, sticky=tk.W)
        bn = tk.Checkbutton(top, text="Background and Noise",
                            variable=master.bck_noise, onvalue=1, offvalue=0)
        bn.grid(row=5, column=0, sticky=tk.W)
        bck = tk.Checkbutton(top, text=u"\u00B9Background subtracted " +
                             "Intensities", variable=master.bck_sub,
                             onvalue=1, offvalue=0)
        bck.grid(row=2, column=1, sticky=tk.W)
        exp = tk.Checkbutton(top, text=u"\u00B9Gaussian Intensities",
                             variable=master.gauss_int, onvalue=1, offvalue=0)
        exp.grid(row=3, column=1, sticky=tk.W)
        button = tk.Button(top, text='Ok', command=self.close)
        button.grid(row=6, column=0, columnspan=2)
        top.lift()

        self.master = master
        self.top = top

    def close(self):
        self.master.output_window_open.set(0)
        self.top.destroy()

    def select_all(self):
        """Set all variables to on (1).
        """
        self.master.abs_int.set(1)
        self.master.rel_int.set(1)
        self.master.gauss_int.set(1)
        self.master.bck_sub.set(1)
        self.master.bck_noise.set(1)
        self.master.peak_qual.set(1)

    def select_none(self):
        """Set all variables to off (0).
        """
        self.master.abs_int.set(0)
        self.master.rel_int.set(0)
        self.master.gauss_int.set(0)
        self.master.bck_sub.set(0)
        self.master.bck_noise.set(0)
        self.master.peak_qual.set(0)
