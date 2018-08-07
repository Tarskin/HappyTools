import os

import Trace

class Chromatogram(object):
    def __init__(self, filename):
        self.filename = filename
        self.data = Trace.Trace().open_chrom(filename)

    def plot_data(self, data, fig, canvas):
        """Plot all chromatograms in data on the canvas.

        This function first clears the canvas and then draws all the
        chromatograms that are in the data list on the canvas.

        Keyword arguments:
        fig -- matplotlib figure object
        canvas -- tkinter canvas object
        data -- list of tuples, consisting of two numbers per tuple
        """
        """ TODO
        """
        fig.clear()
        axes = fig.add_subplot(111)
        for i in data:
            x_array = []
            y_array = []
            for j in i.data:
                x_array.append(j[0])
                y_array.append(j[1])
            label = os.path.splitext(os.path.basename(i.filename))[0]
            line, = axes.plot(x_array,y_array,label=str(label))
        axes.set_xlabel("Time [m]")
        axes.set_ylabel("Intensity [au]")
        handles, labels = axes.get_legend_handles_labels()
        fig.legend(handles,labels)
        axes.get_xaxis().get_major_formatter().set_useOffset(False)
        canvas.draw()
