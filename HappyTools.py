#! /usr/bin/env python
#
# Copyright 2017-2018 Bas C. Jansen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# You should have received a copy of the Apache 2.0 license along
# with this program; if not, see 
# http://www.apache.org/licenses/LICENSE-2.0

# General imports
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
try:
    # Python 2
    import Tkinter as tk
    import tkFileDialog as filedialog
except ImportError:
    # Python 3
    import tkinter as tk
    import tk.filedialog as filedialog
import glob
import matplotlib
import os
import sys
import tkMessageBox

# Custom libraries
import HappyTools.plugins as plugins
import HappyTools.libs.functions as functions
from HappyTools.util.Functions import Functions

# Gui elements
from HappyTools.gui.CustomToolbar import CustomToolbar
from HappyTools.gui.AboutWindow import AboutWindow
from HappyTools.gui.BatchWindow import batchWindow
import HappyTools.gui.Settings as settings
import HappyTools.gui.Version as version
import HappyTools.gui.Debug as debug

# Class imports
from HappyTools.bin.Chromatogram import Chromatogram
from HappyTools.bin.Trace import Trace

# Directories
directories = [
    os.path.join(os.getcwd(),"HappyTools","plugins"),
    os.path.join(os.getcwd(),"HappyTools","gui"),
    os.path.join(os.getcwd(),"HappyTools","bin"),
    os.path.join(os.getcwd(),"HappyTools","util")
]

# Function overwrites
def dynamic_update(foo):
    pass
matplotlib.backends.backend_tkagg.NavigationToolbar2TkAgg.dynamic_update = dynamic_update

# Applicatiom
class HappyToolsGui(object):
    @classmethod
    def run(cls):
        root = tk.Tk()
        app = HappyToolsGui(root)
        root.mainloop()

    def __init__(self, master):
        # ACCESS CHECK
        if not functions.checkAccess(directories):
            tkMessageBox.showinfo("Access Error", "HappyTools does not have sufficient disk access rights. Please close "+
                    "HappyTools and check if the current user has read/write access to all folders in the Happytools "+
                    "folder.")
        # SETTINGS
        #if os.path.isfile(settings):
        #    functions.getSettings()

        # CANVAS
        self.fig = matplotlib.figure.Figure(figsize=(12, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.toolbar = CustomToolbar(self.canvas, master)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=tk.YES)
        self.canvas.draw()

        # FRAME
        frame = tk.Frame(master)
        master.title("HappyTools "+str(version.version)+" (Build "+str(version.build)+")")
        iconbitmap = os.path.join(os.getcwd(),"HappyTools","gui","assets","Icon.ico")
        backgroundimage = os.path.join(os.getcwd(),"HappyTools","gui","assets","UI.png")
        if os.path.isfile(iconbitmap):
            master.iconbitmap(default=iconbitmap)
        if os.path.isfile(backgroundimage):
            background_image = self.fig.add_subplot(111)
            image = matplotlib.image.imread(backgroundimage)
            background_image.axis('off')
            self.fig.set_tight_layout(True)
            background_image.imshow(image)
        
        # QUIT
        def close():
            master.destroy()
            master.quit()
        master.protocol("WM_DELETE_WINDOW", lambda: close())

        # MENU
        menu = tk.Menu(master)
        master.config(menu=menu)

        filemenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open Chromatogram", command=self.open_chromatogram_window)
        filemenu.add_command(label="Smooth Chromatogram", command=self.smooth_chromatogram)
        filemenu.add_command(label="Baseline Correction", command=self.baseline_correction)
        filemenu.add_command(label="Normalize chromatogram", command=self.normalize_chromatogram)
        filemenu.add_command(label="Save Chromatogram", command=self.save_chromatogram)

        processmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Process", menu=processmenu)
        processmenu.add_command(label="Calibrate Chromatogram", command=self.calibrate_chromatogram)
        processmenu.add_command(label="Quantify Chromatogram", command=self.foo)
        
        advancedmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Advanced", menu=advancedmenu)
        advancedmenu.add_command(label="Peak Detection", command=self.foo)
        advancedmenu.add_command(label="Save Calibrants", command=self.foo)
        advancedmenu.add_command(label="Save Annotation", command=self.foo)

        batchmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Batch", menu=batchmenu)
        batchmenu.add_command(label="Batch Process", command=self.open_batch_window)

        settingsmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Settings", menu=settingsmenu)
        settingsmenu.add_command(label="Settings", command=self.foo)

        aboutmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="About", menu=aboutmenu)
        aboutmenu.add_command(label="About HappyTools", command= self.open_about_window)

        #if glob.glob(os.path.join(".","plugins","*.py")):
            #import importlib
            #pluginsmenu = tk.Menu(menu,tearoff=0)
            #menu.add_cascade(label="Plugins", menu=pluginsmenu)
            #for file in glob.glob(os.path.join(".","plugins","*.py")):
                #moduleName = os.path.split(file)[-1].split(".")[0]
                #module = importlib.import_module(moduleName)
                #pluginsmenu.add_command(label=moduleName, command=functions.makeFunc(module))

    def open_about_window(self):
        AboutWindow(tk.Toplevel())

    def open_chromatogram_window(self):
        files = filedialog.askopenfilenames(title="Open Chromatogram File(s)")
        data = []
        if files:
            for file in files:
                foo = Chromatogram(file)
                data.append(foo)
            self.data = data
            self.data[0].plot_data(data, self.fig, self.canvas)

    def calibrate_chromatogram(self):
        cal_file = tk.StringVar()
        cal_file = filedialog.askopenfilename(title="Select Calibration File")
        calibrants = Functions().read_peak_list(cal_file)
        for data in self.data:
            time_pairs = Functions().find_peak(calibrants, data)
            function = Functions().determine_calibration_function(time_pairs)
            data = Functions().apply_calibration_function(function, data)
        self.data[0].plot_data(self.data, self.fig, self.canvas)

    def open_batch_window(self):
        batchWindow(self)

    def normalize_chromatogram(self):
        try:
            self.data = Trace().norm_chrom(self.data)
            self.data[0].plot_data(self.data, self.fig, self.canvas)
        except AttributeError:
            pass

    def smooth_chromatogram(self):
        try:
            self.data = Trace().smooth_chrom(self.data)
            self.data[0].plot_data(self.data, self.fig, self.canvas)
        except AttributeError:
            pass

    def save_chromatogram(self):
        try:
            for data in self.data:
                Trace().save_chrom(data)
        except AttributeError:
            pass

    def baseline_correction(self):
        try:
            self.data = Trace().baseline_correction(self.data)
            self.data[0].plot_data(self.data, self.fig, self.canvas)
        except AttributeError:
            pass

    def foo(self):
        raise NotImplementedError("This feature is not implemented in the refactor yet.")

# Call the main app
if __name__ == "__main__":
    HappyToolsGui.run()
