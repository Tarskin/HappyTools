import HappyTools.gui.version as version
try:
    # Python 2
    import Tkinter as tk
except ImportError:
    # Python 3
    import tkinter as tk


class AboutWindow(object):
    # TODO: Get the version and build from HappyTools.py in here
    INFORMATION = ("HappyTools Version "+str(version.version)+" (build " +
                   str(version.build)+") by Bas Cornelis Jansen (bas.c." +
                   "jansen@gmail.com).\n\n This software is released under " +
                   "the Apache 2.0 License. Full details regarding this " +
                   "license can be found at the following URL:\n\n" +
                   "http://www.apache.org/licenses/LICENSE-2.0")

    def __init__(self):
        root = tk.Toplevel()
        frame = tk.Frame(root)

        root.title("HappyTools")
        root.resizable(width=False, height=False)
        root.protocol("WM_DELETE_WINDOW", self.close)
        
        about = tk.Label(frame, text=self.INFORMATION, justify=tk.LEFT,
                         wraplength=250)

        about.pack()
        frame.pack()
        root.lift()

        self.root = root

    def close(self):
        self.root.destroy()
