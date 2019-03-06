from bisect import bisect_left, bisect_right
from numpy import greater, less, linspace, poly1d, polyfit
from os import W_OK, access
from pathlib import Path
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal import argrelextrema
import logging

def check_disk_access(master):
    disk_access = True
    for directory in master.directories:
        if not access(directory, W_OK):
            disk_access = False
    return disk_access

def determine_breakpoints(master):
    time, intensity = zip(*master.chrom_data)
    low = bisect_left(time, master.time-master.window)
    high = bisect_right(time, master.time+master.window)

    f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])
    f_prime = f.derivative()

    new_x = linspace(time[low], time[high], 2500*(time[high]-time[low]))
    new_prime_y = f_prime(new_x)

    maxm = argrelextrema(new_prime_y, greater)
    minm = argrelextrema(new_prime_y, less)

    breaks = maxm[0].tolist() + minm[0].tolist()
    breaks = sorted(breaks)

    return breaks
    
def read_peak_list(file_name):
    '''Read and parse the peak file and return a list of peaks.

    This function opens the file that is specified in 'fileName', and
    reads it on a line per line basis. The function will split each
    line on '\t' prior to trying to append the parts to the 'peaks'
    list. The function will write a message to the logFile if logging
    is enabled if the previous mentioned try goes to the except clause.

    Keyword argments:
    fileName -- string
    '''
    logger = logging.getLogger(__name__)
    peaks = []
    try:
        with open(file_name, 'r') as fr:
            for line in fr:
                line = line.rstrip('\n').split('\t')
                try:
                    peaks.append((str(line[0]), float(line[1]),
                        float(line[2])))
                except ValueError:
                    logger.info('Ignoring line: '+str(line)+' from '+
                                     'file: '+str(file_name))
    except IOError:
        logger.error('The selected reference file '+str(file_name)+
                          ' could not be opened.')
    return peaks
        
def subset_data(master):

    max_point = 0
    time, intensity = zip(*master.chrom_data)
    low = bisect_left(time, master.time-master.window)
    high = bisect_right(time, master.time+master.window)

    f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])
    new_x = linspace(time[low], time[high], 2500*(time[high]-time[low]))
    new_y = f(new_x)

    x_data = new_x
    y_data = new_y

    # Subset the data
    # Region from newY[0] to breaks[0]
    try:
        if max(new_y[0:master.breaks[0]]) > max_point:
            max_point = max(new_y[0:master.breaks[0]])
            x_data = new_x[0:master.breaks[0]]
            y_data = new_y[0:master.breaks[0]]
    except IndexError:
        pass

    # Regions between breaks[x] and breaks[x+1]
    try:
        for index, _ in enumerate(master.breaks):
            if max(new_y[master.breaks[index]:master.breaks[index+1]]) > max_point:
                max_point = max(new_y[master.breaks[index]:
                    master.breaks[index+1]])
                x_data = new_x[master.breaks[index]:master.breaks[index+1]]
                y_data = new_y[master.breaks[index]:master.breaks[index+1]]
    except IndexError:
        pass

    # Region from break[-1] to newY[-1]
    try:
        if max(new_y[master.breaks[-1]:-1]) > max_point:
            max_point = max(new_y[master.breaks[-1]:-1])
            x_data = new_x[master.breaks[-1]:-1]
            y_data = new_y[master.breaks[-1]:-1]
    except IndexError:
        pass

    return list(zip(x_data, y_data))

def save_calibrants(master):
    # TODO: Remove this function once this is implemented elsewhere
    raise NotImplementedError('This feature is not implemented in the ' +
                              'refactor yet.')
