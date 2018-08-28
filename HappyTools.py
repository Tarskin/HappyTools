#! /usr/bin/env python
#
# Copyright 2017-2018 Bas C. Jansen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.
#
# See the License for the specific language governing permissions and
# limitations under the License.
#
# You should have received a copy of the Apache 2.0 license along
# with this program; if not, see
# http://www.apache.org/licenses/LICENSE-2.0

# General imports
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg
)
try:
    # Python 2
    import Tkinter as tk
    import tkFileDialog as filedialog
except ImportError:
    # Python 3
    import tkinter as tk
    import tk.filedialog as filedialog
import tkMessageBox
from matplotlib import image, figure
from os import path, getcwd

# Custom libraries
import HappyTools.plugins as plugins
from HappyTools.util.functions import Functions
from HappyTools.util.output import Output

# Gui elements
from HappyTools.gui.custom_toolbar import CustomToolbar
from HappyTools.gui.about_window import AboutWindow
from HappyTools.gui.batch_window import batchWindow
from HappyTools.gui.settings import Settings
from HappyTools.gui.output_window import OutputWindow
import HappyTools.gui.progress_bar as progressbar
import HappyTools.gui.version as version

# Class imports
from HappyTools.bin.chromatogram import Chromatogram
from HappyTools.bin.trace import Trace

# Directories
directories = [
    path.join(getcwd(), "HappyTools", "plugins"),
    path.join(getcwd(), "HappyTools", "gui"),
    path.join(getcwd(), "HappyTools", "bin"),
    path.join(getcwd(), "HappyTools", "util")
]


# Function overwrites
def dynamic_update(foo):
    pass
NavigationToolbar2TkAgg.dynamic_update = dynamic_update


# Applicatiom
class HappyToolsGui(object):
    @classmethod
    def run(cls):
        root = tk.Tk()
        HappyToolsGui(root)
        root.mainloop()

    def __init__(self, master):
        self.output_window = tk.IntVar(value=0)
        self.batch_folder = tk.StringVar(value=getcwd())
        self.abs_int = tk.IntVar(value=0)
        self.rel_int = tk.IntVar(value=0)
        self.gauss_int = tk.IntVar(value=0)
        self.bck_sub = tk.IntVar(value=0)
        self.bck_noise = tk.IntVar(value=0)
        self.peak_qual = tk.IntVar(value=0)
        self.create_figure = "True"

        self.master = master
        self.counter = tk.IntVar(value=0)
        self.functions = Functions(self)

        # ACCESS CHECK
        self.directories = directories
        if not self.functions.check_disk_access(self):
            tkMessageBox.showinfo(
                "Access Error", "HappyTools does " +
                "not have sufficient disk access rights. Please close " +
                "HappyTools and check if the current user has read/" +
                "write access to all folders in the Happytools folder.")

        # SETTINGS
        self.settings = Settings(self)
        if path.isfile(path.join(getcwd(), self.settings.settings)):
            self.settings.read_settings(self.settings)

        # CANVAS
        self.fig = figure.Figure(figsize=(12, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.toolbar = CustomToolbar(self.canvas, master)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=tk.YES)
        self.canvas.draw()

        # FRAME
        tk.Frame(master)
        master.title("HappyTools "+str(version.version) +
                     " (Build "+str(version.build)+")")
        iconbitmap = path.join(getcwd(), "HappyTools", "gui", "assets",
                               "Icon.ico")
        backgroundimage = path.join(getcwd(), "HappyTools", "gui",
                                    "assets", "UI.png")
        if path.isfile(iconbitmap):
            master.iconbitmap(default=iconbitmap)
        if path.isfile(backgroundimage):
            background_image = self.fig.add_subplot(111)
            img = image.imread(backgroundimage)
            background_image.axis('off')
            self.fig.set_tight_layout(True)
            background_image.imshow(img)
        self.progress = progressbar.SimpleProgressBar(self)
        self.progress.bar.pack(fill=tk.X)

        # QUIT
        master.protocol("WM_DELETE_WINDOW", self.close)

        # MENU
        menu = tk.Menu(master)
        master.config(menu=menu)

        filemenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open Chromatogram",
                             command=self.open_chromatogram_window)
        filemenu.add_command(label="Smooth Chromatogram",
                             command=self.smooth_chromatogram)
        filemenu.add_command(label="Baseline Correction",
                             command=self.baseline_correction)
        filemenu.add_command(label="Normalize chromatogram",
                             command=self.normalize_chromatogram)
        filemenu.add_command(label="Save Chromatogram",
                             command=self.save_chromatogram)

        processmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Process", menu=processmenu)
        processmenu.add_command(label="Calibrate Chromatogram",
                                command=self.calibrate_chromatogram)
        processmenu.add_command(label="Quantify Chromatogram",
                                command=self.quantify_chromatogram)
        processmenu.add_command(label="Select Outputs",
                                command=self.select_outputs)

        advancedmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Advanced", menu=advancedmenu)
        advancedmenu.add_command(label="Peak Detection",
                                 command=self.peak_detection)
        advancedmenu.add_command(label="Save Calibrants",
                                 command=self.save_calibrants)
        advancedmenu.add_command(label="Save Annotation",
                                 command=self.save_annotation)

        batchmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Batch", menu=batchmenu)
        batchmenu.add_command(label="Batch Process",
                              command=self.open_batch_window)

        settingsmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Settings", menu=settingsmenu)
        settingsmenu.add_command(label="Settings",
                                 command=self.open_settings_window)

        aboutmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="About", menu=aboutmenu)
        aboutmenu.add_command(label="About HappyTools",
                              command=self.open_about_window)

    @classmethod
    def open_about_window(cls):
        AboutWindow()

    def open_chromatogram_window(self):
        files = filedialog.askopenfilenames(title="Open Chromatogram File(s)")
        data = []
        if files:
            for file in files:
                foo = Chromatogram(file)
                data.append(foo)
            self.data = data
            #self.data[0].plot_data(self)
        for chrom in self.data:
            chrom.plot_data_new(self)

    def open_settings_window(self):
        self.settings.settings_popup(self.settings)

    def calibrate_chromatogram(self):
        try:
            self.cal_file = filedialog.askopenfilename(
                title="Select Calibration File")
            self.reference = self.functions.read_peak_list(self.cal_file)

            self.progress.reset_bar(self)
            for index, data in enumerate(self.data):

                progress = (float(index) / len(self.data))*100
                self.counter.set(progress)
                self.progress.update_progress_bar(self)

                self.time_pairs = self.functions.find_peak(self, data)
                self.function = self.functions.determine_calibration_function(self)
                data = self.functions.apply_calibration_function(self, data)

            self.progress.fill_bar(self)
            self.data[0].plot_data(self)
        except AttributeError:
            pass

    def close(self):
        self.master.destroy()
        self.master.quit()

    def open_batch_window(self):
        batchWindow(self)

    def normalize_chromatogram(self):
        try:
            for chrom in self.data:
                chrom.trace.norm_chrom(self)
                chrom.plot_data_new(self)
        except AttributeError:
            pass

    def quantify_chromatogram(self):
        try:
            self.results = []
            self.quant_file = filedialog.askopenfilename(
                title="Select Quantitation File")
            self.reference = self.functions.read_peak_list(self.quant_file)

            self.progress.reset_bar(self)
            for index, data in enumerate(self.data):

                progress = (float(index) / len(self.data))*100
                self.counter.set(progress)
                self.progress.update_progress_bar(self)

                self.results.append({'file': path.basename(data.filename),
                                     'results': self.functions.quantify_chrom(self, data)})

            self.output = Output(self)
            self.output.init_output_file(self)
            self.output.build_output_file(self)
            self.progress.fill_bar(self)
        except AttributeError:
            pass

    def peak_detection(self):
        self.functions.peak_detection(self)

    def save_annotation(self):
        self.functions.save_annotation(self)

    def save_calibrants(self):
        self.functions.save_calibrants(self)

    @classmethod
    def select_outputs(self):
        OutputWindow(self)

    def smooth_chromatogram(self):
        try:
            for chrom in self.data:
                chrom.trace.smooth_chrom(self)
                chrom.plot_data_new(self)
        except AttributeError:
            pass

    def save_chromatogram(self):
        try:
            Trace().save_chrom(self)
        except AttributeError:
            pass

    def baseline_correction(self):
        try:
            for chrom in self.data:
                chrom.trace.baseline_correction(self)
                chrom.plot_data_new(self)
        except AttributeError:
            pass

# Call the main app
if __name__ == "__main__":
    HappyToolsGui.run()
