import re
import numpy as np

class Trace(object):
    def openChrom(self, file):
        """Read a chromatogram and return the data.

        This function opens a chromatogram (txt or arw), interprets the
        local thousands/decimal seperators and creates a list of retention
        time and intensity tuples which is returned.

        Keyword arguments:
        file -- unicode string
        """
        with open(file,'r') as fr:
            chromData = []
            if 'txt' in file:
                for line in fr:
                    if line[0].isdigit() == True:
                        lineChunks = line.strip().split()
                        # Number based regex splitting to get rid of thousand seperators
                        timeSep = re.sub(r'-?\d', '', lineChunks[0], flags=re.U)
                        for sep in timeSep[:-1]:
                            lineChunks[0] = lineChunks[0].replace(sep, '')
                        if timeSep:
                            lineChunks[0] = lineChunks[0].replace(timeSep[-1], '.')
                        intSep = re.sub(r'-?\d', '', lineChunks[-1], flags=re.U)
                        for sep in intSep[:-1]:
                            lineChunks[-1] = lineChunks[-1].replace(sep[-1], '')
                        if intSep:
                            lineChunks[-1] = lineChunks[-1].replace(intSep[-1], '.')
                        # End of regex based splitting
                        try:
                            chromData.append((float(lineChunks[0]),float(lineChunks[-1])))
                        except UnicodeEncodeError:
                            print("Omitting line: "+str(line))
            elif 'arw' in file:
                for line in fr:
                    lines = line.split('\r')
                for line in lines:
                    try:
                        if line[0][0].isdigit() == False:
                            pass
                        else:
                            chunks = line.rstrip()
                            chunks = chunks.split()
                            chromData.append((float(chunks[0]),float(chunks[1])))
                    except IndexError:
                        pass
            else:
                print("Incorrect inputfile format, please upload a raw data 'txt' or 'arw' file.")
        return chromData
