import tkinter as tk

def start(master):
    def close():
        top.destroy()

    top = tk.Tk()
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
    about = tk.Label(top, text=information, justify=tk.LEFT, wraplength=250)
    about.pack()
    top.lift()
