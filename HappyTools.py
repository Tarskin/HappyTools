#! /usr/bin/env python
#
# Copyright 2017-2019 Bas C. Jansen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.
#
# See the License for the specific language governing permissions and
# limitations under the License.
#
# You should have received a copy of the Apache 2.0 license along
# with this program; if not, see
# http://www.apache.org/licenses/LICENSE-2.0

# Compatability check
import HappyTools.util.requirement_checker as req_check
req_check.check_requirements()

# General imports
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk
)
import importlib
import logging
import os
import tkinter as tk
import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox
import matplotlib
from matplotlib import image, figure
from pathlib import Path, PurePath

# Platform specific bits
if os.name == 'posix':
    import matplotlib
    matplotlib.use('TkAgg')

# Custom libraries
from HappyTools.util.peak_detection import PeakDetection
from HappyTools.util.functions import (Functions, check_disk_access,
                                       determine_calibration_function,
                                       apply_calibration_function)
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
from HappyTools.bin.chromatogram import Chromatogram, finalize_plot

# Directories
directories = [
    Path.cwd() / 'HappyTools' / 'plugins',
    Path.cwd() / 'HappyTools' / 'gui',
    Path.cwd() / 'HappyTools' / 'bin',
    Path.cwd() / 'HappyTools' / 'util'
]


# Function overwrites
def dynamic_update(foo):
    pass
NavigationToolbar2Tk.dynamic_update = dynamic_update


# Applicatiom
class HappyToolsGui(object):
    @classmethod
    def run(cls):
        root = tk.Tk()
        HappyToolsGui(root)
        root.mainloop()

    def __init__(self, master):
        self.output_window_open = tk.IntVar(value=0)
        self.batch_folder = tk.StringVar(value=Path.cwd())
        self.abs_int = tk.IntVar(value=1)
        self.rel_int = tk.IntVar(value=1)
        self.gauss_int = tk.IntVar(value=1)
        self.bck_sub = tk.IntVar(value=1)
        self.bck_noise = tk.IntVar(value=1)
        self.peak_qual = tk.IntVar(value=1)
        self.create_figure = 'True'

        self.master = master
        self.counter = tk.DoubleVar(value=0.0)
        logging.basicConfig(filename='HappyTools.log',
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %H:%M', filemode='a',
                            level=logging.WARNING)
        self.logger = logging.getLogger(__name__)
        self.functions = Functions(self)

        # ACCESS CHECK
        self.directories = directories
        if not check_disk_access(self):
            messagebox.showinfo(
                'Access Error', 'HappyTools does ' +
                'not have sufficient disk access rights. Please close ' +
                'HappyTools and check if the current user has read/' +
                'write access to all folders in the Happytools folder.')

        # SETTINGS
        self.settings = Settings(self)
        if (Path.cwd() / self.settings.settings).is_file():
            self.settings.read_settings()

        # CANVAS
        self.fig = figure.Figure(figsize=(12,6))
        #self.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
        self.axes = self.fig.add_subplot(111)
        self.axes.axis('off')
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.toolbar = CustomToolbar(self.canvas, master)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=tk.YES)
        self.canvas.draw()

        # FRAME
        tk.Frame(master)
        master.title('HappyTools '+str(version.version) +
                     ' (Build '+str(version.build)+')')
        iconbitmap = Path.cwd() / 'HappyTools' / 'gui' / 'assets' / 'Icon.ico'
        backgroundimage = Path.cwd() / 'HappyTools' / 'gui' / 'assets' / 'UI.png'
        try:
            master.iconbitmap(default=iconbitmap)
        except tk.TclError as e:
            self.logger.warn(e)
        if backgroundimage.is_file():
            img = image.imread(str(backgroundimage))
            self.axes.imshow(img)
            self.axes.set_aspect('auto')
        self.progress = progressbar.SimpleProgressBar(self)
        self.progress.bar.pack(fill=tk.X)

        # QUIT
        master.protocol('WM_DELETE_WINDOW', self.close)

        # MENU
        menu = tk.Menu(master)
        master.config(menu=menu)

        filemenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='File', menu=filemenu)
        filemenu.add_command(label='Open Chromatogram',
                             command=self.open_chromatogram_window)
        filemenu.add_command(label='Smooth Chromatogram',
                             command=self.smooth_chromatogram)
        filemenu.add_command(label='Baseline Correction',
                             command=self.baseline_correction)
        filemenu.add_command(label='Normalize chromatogram',
                             command=self.normalize_chromatogram)
        filemenu.add_command(label='Save Chromatogram',
                             command=self.save_chromatogram)

        processmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Process', menu=processmenu)
        processmenu.add_command(label='Calibrate Chromatogram',
                                command=self.calibrate_chromatogram)
        processmenu.add_command(label='Quantify Chromatogram',
                                command=self.quantify_chromatogram)
        processmenu.add_command(label='Select Outputs',
                                command=self.open_output_window)

        advancedmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Advanced', menu=advancedmenu)
        advancedmenu.add_command(label='Peak Detection',
                                 command=self.peak_detection)
        advancedmenu.add_command(label='Save Calibrants',
                                 command=self.save_calibrants)
        advancedmenu.add_command(label='Save Annotation',
                                 command=self.save_annotation)

        batchmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Batch', menu=batchmenu)
        batchmenu.add_command(label='Batch Process',
                              command=self.open_batch_window)

        settingsmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Settings', menu=settingsmenu)
        settingsmenu.add_command(label='Settings',
                                 command=self.open_settings_window)

        aboutmenu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='About', menu=aboutmenu)
        aboutmenu.add_command(label='About HappyTools',
                              command=self.open_about_window)

        if (Path.cwd() / 'HappyTools' / 'plugins').glob("*.py"):
            pluginsmenu = tk.Menu(menu,tearoff=0)
            menu.add_cascade(label="Plugins", menu=pluginsmenu)
            for file in(Path.cwd() / 'HappyTools' / 'plugins').glob('*.py'):
                if '__' in str(file):
                    continue
                module_name = PurePath(file).stem
                module = "HappyTools.plugins."+str(module_name)
                module = importlib.import_module(module)
                try:
                    module_name = module.name
                except Exception as e:
                    self.logger.error(e)
                pluginsmenu.add_command(label=module_name,
                                        command=self.make_function(module))

    def make_function(self, module):
        try:
            def x():
                return module.start(self)
        except AttributeError as e:
            self.logger.error('Problem with the plugin: '+str(e))
        return x

    @classmethod
    def open_about_window(cls):
        AboutWindow()

    def open_chromatogram_window(self):
        files = filedialog.askopenfilenames(title='Open Chromatogram File(s)')
        data = []
        if files:
            for file in files:
                foo = Chromatogram(Path(file))
                data.append(foo)
            self.data = data

        self.axes.clear()
        for chrom in self.data:
            chrom.plot_chrom(self)
        finalize_plot(self)

    def open_output_window(self):
        OutputWindow(self)

    def open_settings_window(self):
        self.settings.settings_popup(self.settings)

    def calibrate_chromatogram(self):
        try:
            self.cal_file = filedialog.askopenfilename(
                title='Select Calibration File')
            if not self.cal_file:
                self.cal_file = None
            self.reference = self.functions.read_peak_list(self.cal_file)

            self.progress.reset_bar()
            for index, self.chrom in enumerate(self.data):

                progress = (float(index) / len(self.data))*100
                self.counter.set(progress)
                self.progress.update_progress_bar()

                self.time_pairs = self.functions.find_peak(self)
                self.function = determine_calibration_function(self)
                apply_calibration_function(self)

        except Exception as e:
            self.logger.error(e)
        self.progress.fill_bar()

        self.axes.clear()
        for chrom in self.data:
            chrom.plot_chrom(self)
        finalize_plot(self)

    def close(self):
        self.master.destroy()
        self.master.quit()

    def open_batch_window(self):
        batchWindow(self)

    def normalize_chromatogram(self):
        try:
            self.axes.clear()
            for chrom in self.data:
                chrom.trace.norm_chrom(self)
                chrom.plot_chrom(self)
            finalize_plot(self)
        except Exception as e:
            self.logger.error(e)

    def quantify_chromatogram(self):
        try:
            self.results = []
            self.functions.batch_folder = self.batch_folder
            self.quant_file = filedialog.askopenfilename(
                title='Select Quantitation File')
            if not self.quant_file:
                self.quant_file = None
            self.reference = self.functions.read_peak_list(self.quant_file)

            self.progress.reset_bar()
            for index, self.chrom in enumerate(self.data):

                progress = (float(index) / len(self.data))*100
                self.counter.set(progress)
                self.progress.update_progress_bar()

                self.results.append({'file': Path(self.chrom.filename).name,
                                     'results': self.functions.quantify_chrom(self)})

            self.output = Output(self)
            self.output.init_output_file(self)
            self.output.build_output_file(self)
            self.progress.fill_bar()
        except Exception as e:
            self.logger.error(e)

    def peak_detection(self):
        try:
            self.axes.clear()
            for self.chrom in self.data:
                self.detected_peaks = PeakDetection(self)
                self.detected_peaks.detect_peaks(self)
                self.detected_peaks.plot_peaks(self)
                self.chrom.plot_chrom(self)
            finalize_plot(self)
        except Exception as e:
            self.logger.error(e)

    def save_annotation(self):
        try:
            for self.chrom in self.data:
                self.detected_peaks.write_peaks(self)
        except Exception as e:
            self.logger.error(e)

    def save_calibrants(self):
        self.functions.save_calibrants(self)

    def smooth_chromatogram(self):
        try:
            self.axes.clear()
            for chrom in self.data:
                chrom.trace.smooth_chrom(self)
                chrom.plot_chrom(self)
            finalize_plot(self)
        except Exception as e:
            self.logger.error(e)

    def save_chromatogram(self):
        try:
            for chrom in self.data:
                chrom.save_chrom(self)
        except Exception as e:
            self.logger.error(e)

    def baseline_correction(self):
        try:
            self.axes.clear()
            for chrom in self.data:
                chrom.trace.baseline_correction(self)
                chrom.plot_chrom(self)
            finalize_plot(self)
        except Exception as e:
            self.logger.error(e)

# Call the main app
if __name__ == '__main__':
    HappyToolsGui.run()
