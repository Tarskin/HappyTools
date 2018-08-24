from bisect import bisect_left, bisect_right
from scipy.signal import savgol_filter
import re
import operator
import numpy as np


class Trace(object):

    def open_chrom(self, file):
        """Read a chromatogram and return the data.

        This function opens a chromatogram (txt or arw), interprets the
        local thousands/decimal seperators and creates a list of retention
        time and intensity tuples which is returned.

        Keyword arguments:
        file -- unicode string
        """
        with open(file, 'r') as fr:
            chrom_data = []
            if 'txt' in file:
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
            elif 'arw' in file:
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
        return chrom_data

    def baseline_correction(self, master):
        """Perform baseline correction and return the corrected data.

        This function determines the baseline of a chromatogram between
        two timepoints, specified with the start and end parameter. The
        chromatogram is split into segments of a length specified in the
        points parameter. The lowest intensity of each segment is used to
        determine a function of the order specified in the baseline_order
        using the numpy.polyfit function. The original chromatogram is then
        transformed by subtracting the function from the original data. The
        resulting chromatogram might have negative intensities between the
        start and end timepoints, the minimum intensity within that region
        is used to uplift the entire chromatogram. The transformed and
        uplifted chromatogram is returned to the calling function.

        Keyword arguments:
        data --
        """

        for i in master.data:
            # Background determination
            background = []
            chunks = [i.data[x:x+master.settings.points] for x in xrange(
                0, len(i.data), master.settings.points)]
            for j in chunks:
                buff1, buff2 = zip(*j)
                min_index, _ = min(enumerate(buff2),
                                           key=operator.itemgetter(1))
                if buff1[0] > master.settings.start and buff1[-1] < master.settings.end:
                    background.append((buff1[min_index], buff2[min_index]))

            # Baseline function
            time, intensity = zip(*background)
            func = np.polyfit(time, intensity, master.settings.baseline_order)
            p = np.poly1d(func)

            # Transform
            time = [a for a, b in i.data]
            new_chrom_intensity = [b-p(a) for a, b in i.data]

            # Uplift
            low = bisect_left(time, master.settings.start)
            high = bisect_right(time, master.settings.end)
            offset = abs(min(min(new_chrom_intensity[low:high]), 0))
            i.data = zip(time, [x+offset for x in new_chrom_intensity])

        # Return
        return master.data

    def smooth_chrom(self, master):
        """ TODO
        """
        # Apply Savitzky-Golay filter
        for i in master.data:
            time, intensity = zip(*i.data)
            new = savgol_filter(intensity, 21, 3)
            i.data = zip(time, new)

        # Return
        return master.data

    def norm_chrom(self, master):
        """ TODO
        """
        # Normalize to maximum intensity
        for i in master.data:
            time, intensity = zip(*i.data)
            maximum = max(intensity[bisect_left(
                time, master.settings.start):bisect_right(
                time, master.settings.end)])
            normalized_intensity = [b/maximum for a, b, in i.data]
            i.data = zip(time, normalized_intensity)

        # Return
        return master.data

    def save_chrom(self, master):
        """ TODO
        """
        for i in master.data:
            with open(i.filename, 'w') as fw:
                for data_point in i.data:
                    fw.write(
                        str(format(data_point[0], '0.' +
                            str(master.settings.decimal_numbers)+'f'))+"\t" +
                        str(format(data_point[1], '0.' +
                            str(master.settings.decimal_numbers)+'f'))+"\n")
