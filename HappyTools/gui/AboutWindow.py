try:
    # Python 2
    import Tkinter as tk
except ImportError:
    # Python 3
    import tkinter as tk
import HappyTools.gui.Version as version

class AboutWindow(object):
    # TODO: Get the version and build from HappyTools.py in here
    INFORMATION = ("HappyTools Version "+str(version.version)+" (build "+str(version.build)+") " +
                   "by Bas Cornelis Jansen (bas.c.jansen@gmail.com).\n\n" +
                   "This software is released under the Apache 2.0 License. " +
                   "Full details regarding this license can be found at " +
                   "the following URL:\n\n" +
                   "http://www.apache.org/licenses/LICENSE-2.0")

    def __init__(self, master):
        self.master = master
        self.master.title("HappyTools")
        self.master.resizable(width=False, height=False)
        self.master.protocol("WM_DELETE_WINDOW", self.close)

        self.frame = tk.Frame(self.master)
        about = tk.Label(self.frame, text=self.INFORMATION, justify=tk.LEFT, wraplength=250)
        about.pack()
        self.frame.pack()
        self.master.lift()

    def close(self):
        self.master.destroy()
