import tkinter as tk
import tkinter.ttk as ttk


class ProgressBar(object):
    def __init__(self, master):

        self.master = master

        self.calibration_percentage = tk.StringVar()
        self.quantitation_percentage = tk.StringVar()

        self.calibration_percentage.set("0%")
        self.quantitation_percentage.set("0%")

        self.bar_window = tk.Toplevel()
        self.bar_window.title("Progress Bar")

        self.cal = tk.Label(self.bar_window, text="Calibration", padx=25)
        self.cal.grid(row=0, column=0, sticky=tk.W)
        self.ft = ttk.Frame(self.bar_window)
        self.ft.grid(row=1, columnspan=2, sticky="")
        self.perc1 = tk.Label(self.bar_window,
                              textvariable=self.calibration_percentage)
        self.perc1.grid(row=0, column=1, padx=25)
        self.progressbar = ttk.Progressbar(self.ft, length=100,
                                           mode='determinate')
        self.progressbar.grid(row=1, columnspan=2, sticky="")

        self.ext = tk.Label(self.bar_window, text="Quantitation", padx=25)
        self.ext.grid(row=2, column=0, sticky=tk.W)
        self.ft2 = ttk.Frame(self.bar_window)
        self.ft2.grid(row=3, columnspan=2, sticky="")
        self.perc2 = tk.Label(self.bar_window,
                              textvariable=self.quantitation_percentage)
        self.perc2.grid(row=2, column=1, padx=25)
        self.progressbar2 = ttk.Progressbar(self.ft2, length=100,
                                            mode='determinate')
        self.progressbar2.grid(row=3, columnspan=2, sticky="")

    @staticmethod
    def update_progress_bar(bar, variable, index, length):
        variable.set(str(int((float(index)/float(length))*100))+"%")
        bar["value"] = int((float(index)/float(length))*100)
        bar.update()

    def close(self):
        self.master.destroy()


class SimpleProgressBar(object):
    def __init__(self, master):
        self.master = master

        style = ttk.Style(master.master)
        style.layout('text.Horizontal.TProgressbar', 
                     [('Horizontal.Progressbar.trough',
                      {'children': [('Horizontal.Progressbar.pbar',
                                    {'side': 'left', 'sticky': 'ns'})],
                       'sticky': 'nswe'}), 
                     ('Horizontal.Progressbar.label', {'sticky': ''})])
        style.configure('text.Horizontal.TProgressbar', text='0 %')

        self.bar = ttk.Progressbar(master.master, orient="horizontal",
                                   style='text.Horizontal.TProgressbar',
                                   length=1000, mode="determinate",
                                   variable=master.counter, maximum=100)

        self.style=style

    def update_progress_bar(self):
        self.update_progress_counter()
        self.bar.update()

    def update_progress_counter(self):
        self.style.configure('text.Horizontal.TProgressbar', 
                             text='{:g} %'.format(self.master.counter.get()))

    def reset_bar(self):
        self.master.counter.set(0)
        self.bar.update()

    def fill_bar(self):
        self.master.counter.set(100)
        self.update_progress_counter()
        self.bar.update()
