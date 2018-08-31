import tkinter as tk
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk


class CustomToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas_, parent_):
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move',
                'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            ('Axes', 'Define the X- and Y-axis limits', 'subplots',
                'plot_axes'),
            ('Subplots', 'Configure subplots', 'subplots',
                'configure_subplots'),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )
        NavigationToolbar2Tk.__init__(self, canvas_, parent_)

    def plot_axes(self):
        self.push_current()
        top = tk.top = tk.Toplevel()
        top.title("Configure Axes")
        top.protocol("WM_DELETE_WINDOW", self.close)

        x_label = tk.Label(top, text="X-axis", font="bold")
        x_label.grid(row=0, column=0, sticky=tk.W)
        x_min_window = tk.Entry(top)
        x_min_window.grid(row=0, column=1, sticky=tk.W)
        x_max_window = tk.Entry(top)
        x_max_window.grid(row=0, column=2, sticky=tk.W)
        y_label = tk.Label(top, text="Y-axis", font="bold")
        y_label.grid(row=1, column=0, sticky=tk.W)
        y_min_window = tk.Entry(top)
        y_min_window.grid(row=1, column=1, sticky=tk.W)
        y_max_window = tk.Entry(top)
        y_max_window.grid(row=1, column=2, sticky=tk.W)

        self.top = top
        self.x_min_window = x_min_window
        self.x_max_window = x_max_window
        self.y_min_window = y_min_window
        self.y_max_window = y_max_window

    def close(self):
        try:
            self.canvas.figure.axes[0].set_xlim([float(self.x_min_window.get()),
                float(self.x_max_window.get())])
        except ValueError:
            pass
        try:
            self.canvas.figure.axes[0].set_ylim([float(self.y_min_window.get()),
                float(self.y_max_window.get())])
        except ValueError:
            pass
        self.canvas.draw()
        self.top.destroy()
