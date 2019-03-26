from os import W_OK, access
import logging

def check_disk_access(master):
    disk_access = True
    for directory in master.directories:
        if not access(directory, W_OK):
            disk_access = False
    return disk_access

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

def save_calibrants(master):
    # TODO: Remove this function once this is implemented elsewhere
    raise NotImplementedError('This feature is not implemented in the ' +
                              'refactor yet.')
