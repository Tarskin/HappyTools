from HappyTools.bin.trace import Trace
from pathlib import PurePath


class Chromatogram(object):
    def __init__(self, filename):
        self.filename = filename
        self.trace = Trace(self)
        self.trace.open_chrom(self)

    def plot_chrom(self, master):
        label = PurePath(self.filename).stem
        x_array, y_array = zip(*self.trace.chrom_data)
        master.axes.plot(x_array, y_array, label=str(label))
       
    def save_chrom(self, master):
        with open(self.filename, 'w') as fw:
            for data_point in self.trace.chrom_data:
                fw.write(
                    str(format(data_point[0], '0.' +
                        str(master.settings.decimal_numbers)+'f'))+'\t' +
                    str(format(data_point[1], '0.' +
                        str(master.settings.decimal_numbers)+'f'))+'\n')

def finalize_plot(master):
    master.axes.set_xlabel('Time [m]')
    master.axes.set_ylabel('Intensity [au]')
    handles, labels = master.axes.get_legend_handles_labels()
    master.fig.legend(handles, labels)
    master.axes.get_xaxis().get_major_formatter().set_useOffset(False)
    master.canvas.draw()
