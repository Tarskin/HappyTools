from HappyTools.bin.trace import Trace
from os import path


class Chromatogram(object):
    def __init__(self, filename):
        self.filename = filename
        self.trace = Trace(self)
        self.trace.open_chrom(self)

    def plot_chrom(self, master):
        axes = master.fig.add_subplot(111)
        label = path.splitext(path.basename(self.filename))[0]
        x_array, y_array = zip(*self.trace.chrom_data)
        axes.plot(x_array, y_array, label=str(label))
        axes.set_xlabel("Time [m]")
        axes.set_ylabel("Intensity [au]")
        handles, labels = axes.get_legend_handles_labels()
        master.fig.legend(handles, labels)
        axes.get_xaxis().get_major_formatter().set_useOffset(False)

    def save_chrom(self, master):
        with open(self.filename, 'w') as fw:
            for data_point in self.trace.chrom_data:
                fw.write(
                    str(format(data_point[0], '0.' +
                        str(master.settings.decimal_numbers)+'f'))+"\t" +
                    str(format(data_point[1], '0.' +
                        str(master.settings.decimal_numbers)+'f'))+"\n")
