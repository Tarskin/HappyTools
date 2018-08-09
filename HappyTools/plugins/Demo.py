#! /usr/bin/env python

# General imports
from Tkinter import *

def start():
    def close():
        top.destroy()

    top = Tk.top = Toplevel()
    top.protocol("WM_DELETE_WINDOW", lambda: close())
    top.title("Plugin Demo")
    information = ("This is a demonstration of HappyTools plugin "+
                   "functionality. Please add a single python file "+
                   "to the plugins directory, from which you can link "+
                   "to further files (in sub directories) where "+ 
                   "needed. HappyTools will parse the plugins "+
                   "directory and add a button for each python file "+
                   "that it finds. The program expects a function "+
                   "called start as the entry point for a plugin.")
    about = Label(top, text=information, justify=LEFT, wraplength=250)
    about.pack()
    top.lift()
