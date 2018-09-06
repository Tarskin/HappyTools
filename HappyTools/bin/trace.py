from bisect import bisect_left, bisect_right
from scipy.signal import savgol_filter
import re
import operator
import numpy as np


class Trace(object):
    def __init__(self, master):
        """ TODO
        """
        self.chrom_data = None

    def open_chrom(self, master):
        """ TODO
        """
        with open(master.filename, 'r') as fr:
            chrom_data = []
            if master.filename.suffix.lower().endswith('txt'):
                for line in fr:
                    if line[0].isdigit() is True:
                        line_chunks = line.strip().split()

                        time_sep = re.sub(r'-?\d', '', line_chunks[0],
                                          flags=re.U)
                        for sep in time_sep[:-1]:
                            line_chunks[0] = line_chunks[0].replace(sep, '')
                        if time_sep:
                            line_chunks[0] = line_chunks[0].replace(
                                time_sep[-1], '.')

                        int_sep = re.sub(r'-?\d', '', line_chunks[-1],
                                         flags=re.U)
                        for sep in int_sep[:-1]:
                            line_chunks[-1] = line_chunks[-1].replace(
                                sep[-1], '')
                        if int_sep:
                            line_chunks[-1] = line_chunks[-1].replace(
                                int_sep[-1], '.')

                        try:
                            chrom_data.append((float(line_chunks[0]),
                                               float(line_chunks[-1])))
                        except UnicodeEncodeError:
                            print("Omitting line: "+str(line))
            elif master.filename.lower().endswith('arw'):
                for line in fr:
                    lines = line.split('\r')
                for line in lines:
                    try:
                        if line[0][0].isdigit() is False:
                            pass
                        else:
                            line_chunks = line.rstrip()
                            line_chunks = line_chunks.split()
                            chrom_data.append((float(line_chunks[0]),
                                               float(line_chunks[1])))
                    except IndexError:
                        pass
            else:
                print("Incorrect inputfile format, please upload a raw " +
                      "data 'txt' or 'arw' file.")
        self.chrom_data = chrom_data

    def baseline_correction(self, master):
        """ TODO
        """

        # Background determination
        background = []
        chunks = [self.chrom_data[x:x+master.settings.points] for x in range(
                  0, len(self.chrom_data), master.settings.points)]
        for j in chunks:
            buff1, buff2 = zip(*j)
            min_index, _ = min(enumerate(buff2), key=operator.itemgetter(1))
            if buff1[0] > master.settings.start and buff1[-1] < master.settings.end:
                background.append((buff1[min_index], buff2[min_index]))

        # Baseline function
        time, intensity = zip(*background)
        func = np.polyfit(time, intensity, master.settings.baseline_order)
        p = np.poly1d(func)

        # Transform
        time = [a for a, b in self.chrom_data]
        new_chrom_intensity = [b-p(a) for a, b in self.chrom_data]

        # Uplift
        low = bisect_left(time, master.settings.start)
        high = bisect_right(time, master.settings.end)
        offset = abs(min(min(new_chrom_intensity[low:high]), 0))

        self.chrom_data = list(zip(time,
                               [x+offset for x in new_chrom_intensity]))

    def smooth_chrom(self, master):
        """ TODO
        """
        # Apply Savitzky-Golay filter
        time, intensity = zip(*self.chrom_data)
        new = savgol_filter(intensity, 21, 3)
        self.chrom_data = list(zip(time, new))

    def norm_chrom(self, master):
        """ TODO
        """
        time, intensity = zip(*self.chrom_data)

        # Normalize to maximum intensity
        maximum = max(intensity[bisect_left(
            time, master.settings.start):bisect_right(
            time, master.settings.end)])
        normalized_intensity = [b/maximum for a, b, in self.chrom_data]

        self.chrom_data = list(zip(time, normalized_intensity))
