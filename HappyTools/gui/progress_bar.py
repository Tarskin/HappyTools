import tkinter as tk
import tkinter.ttk as ttk


class ProgressBar(object):
    def __init__(self, master):
        self.master = master
        self.counter = tk.IntVar()

        style = ttk.Style(self.master)
        style.layout('text.Horizontal.TProgressbar',
                     [('Horizontal.Progressbar.trough',
                      {'children': [('Horizontal.Progressbar.pbar',
                                    {'side': 'left', 'sticky': 'ns'})],
                       'sticky': 'nswe'}),
                     ('Horizontal.Progressbar.label', {'sticky': ''})])
        style.configure('text.Horizontal.TProgressbar', text='0 %')

        self.bar = ttk.Progressbar(self.master, orient="horizontal",
                                   style='text.Horizontal.TProgressbar',
                                   length=1000, mode="determinate",
                                   variable=self.counter, maximum=100)

        self.style=style

    def update_progress_bar(self):
        self.update_progress_counter()
        self.bar.update()

    def update_progress_counter(self):
        self.style.configure('text.Horizontal.TProgressbar',
                             text='{:g} %'.format(self.counter.get()))

    def reset_bar(self):
        self.counter.set(0)
        self.update_progress_counter()
        self.bar.update()

    def fill_bar(self):
        self.counter.set(100)
        self.update_progress_counter()
        self.bar.update()
