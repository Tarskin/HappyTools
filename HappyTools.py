#! /usr/bin/env python
#
# Copyright 2017-2017 Bas C. Jansen
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
# You should have received a coyp of the Apache 2.0 license along
# with this program; if not, see 
# http://www.apache.org/licenses/LICENSE-2.0

# General imports
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from Tkinter import *
import glob
import matplotlib
import os

# Custom libraries
sys.path.append('libs')
if glob.glob(os.path.join(".","plugins","*.py")):
    sys.path.append('plugins')

import functions

# Innate variables
version = "0.0.2"
build = "170714a"

# General variables
output = "summary.results"
settings = "HappyTools.ini"

# Debug variables
logFile = "HappyTools.log"
logging = True
logLevel = 2

# Function overwrites
def dynamic_update(foo):
    pass
matplotlib.backends.backend_tkagg.NavigationToolbar2TkAgg.dynamic_update = dynamic_update

# Applicatiom
class App():

    def __init__(self, master):
        # SETTINGS
        if os.path.isfile(settings):
            functions.getSettings()

        # CANVAS
        self.fig = matplotlib.figure.Figure(figsize=(12, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, master)
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=YES)
        self.canvas.draw()

        # FRAME
        frame = Frame(master)
        master.title("HappyTools "+str(version))
        if os.path.isfile(os.path.join(".","ui","Icon.ico")):
            master.iconbitmap(default=os.path.join(".","ui","Icon.ico"))
        if os.path.isfile(os.path.join(".","ui","UI.png")):
            background_image = self.fig.add_subplot(111)
            image = matplotlib.image.imread(os.path.join(".","ui","UI.png"))
            background_image.axis('off')
            self.fig.set_tight_layout(True)
            background_image.imshow(image)
        
        # QUIT
        def close():
            root.destroy()
            root.quit()
        root.protocol("WM_DELETE_WINDOW", lambda: close())

        # MENU
        menu = Menu(master)
        master.config(menu=menu)

        filemenu = Menu(menu, tearoff=0)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open Chromatogram", command=lambda: functions.openFile(self.fig, self.canvas))
        filemenu.add_command(label="Smooth Chromatogram", command=lambda: functions.smoothChrom(self.fig, self.canvas))
        filemenu.add_command(label="Compare Chromatogram", command=lambda: functions.addFile(self.fig, self.canvas))
        filemenu.add_command(label="Baseline Correction", command=lambda: functions.baselineCorrection(self.fig, self.canvas))
        filemenu.add_command(label="Chromatogram Calibration", command=lambda: functions.chromCalibration(self.fig, self.canvas))
        filemenu.add_command(label="Save Chromatogram", command=functions.saveChrom)
        filemenu.add_command(label="Quantify Chromatogram", command=lambda: functions.quantifyChrom(self.fig, self.canvas))

        multimenu = Menu(menu, tearoff=0)
        menu.add_cascade(label="Multi File", menu=multimenu)
        multimenu.add_command(label="Batch Plot", command=lambda: functions.batchPlot(self.fig, self.canvas))

        advancedmenu = Menu(menu, tearoff=0)
        menu.add_cascade(label="Advanced Tools", menu=advancedmenu)
        advancedmenu.add_command(label="Peak Detection", command=lambda: functions.peakDetection(self.fig, self.canvas))
        advancedmenu.add_command(label="Save Annotation", command=lambda: functions.saveAnnotation(self.fig, self.canvas))

        menu.add_command(label="Batch Process", command=functions.batchPopup)

        menu.add_command(label="Settings", command=functions.settingsPopup)

        menu.add_command(label="About MassyTools", command=lambda: functions.infoPopup())

        if glob.glob(os.path.join(".","plugins","*.py")):
            import importlib
            pluginsmenu = Menu(menu,tearoff=0)
            menu.add_cascade(label="Plugins", menu=pluginsmenu)
            for file in glob.glob(os.path.join(".","plugins","*.py")):
                moduleName = os.path.split(file)[-1].split(".")[0]
                module = importlib.import_module(moduleName)
                pluginsmenu.add_command(label=moduleName, command=functions.makeFunc(module))

        # CLEANUP
        functions.fileCleanup()

# Call the main app
if __name__ == "__main__":
    root = Tk()
    app = App(root)
    root.mainloop()
