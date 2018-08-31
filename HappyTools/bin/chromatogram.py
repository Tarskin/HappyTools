from HappyTools.bin.trace import Trace
from os import path


class Chromatogram(object):
    def __init__(self, filename):
        self.filename = filename
        self.trace = Trace(self)
        self.trace.open_chrom(self)

    def plot_data(self, master):
        axes = master.fig.add_subplot(111)
        label = path.splitext(path.basename(self.filename))[0]
        x_array, y_array = zip(*self.trace.chrom_data)
        axes.plot(x_array, y_array, label=str(label))
        axes.set_xlabel("Time [m]")
        axes.set_ylabel("Intensity [au]")
        handles, labels = axes.get_legend_handles_labels()
        master.fig.legend(handles, labels)
        axes.get_xaxis().get_major_formatter().set_useOffset(False)
