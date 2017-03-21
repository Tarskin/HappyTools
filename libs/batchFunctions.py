#! /usr/bin/env python

# General imports
import bisect
import glob
import numpy as np
import operator
import os

# Custom libraries
import functions

# Variables
points = 100
start = 5
end = 30
backgroundOrder = 1

# Defines
CALIBRATION_FILETYPES = ["*.txt","*.arw"]
INTEGRATION_FILETYPES = ["calibrated*.txt"]

# Functions
def batchBaselineCorrection(data):
    """ TODO
    """

    # Background determination
    background = []
    chunks = [data[x:x+points] for x in xrange(0, len(data), points)]
    for i in chunks:
        buff1, buff2 = zip(*i)
        min_index, min_value = min(enumerate(buff2), key=operator.itemgetter(1))
        if buff1[0] > start and buff1[-1] < end:
            background.append((buff1[min_index], buff2[min_index]))
    time, intensity = zip(*background)
    newX = np.linspace(min(time), max(time),100)
    func = np.polyfit(time, intensity, backgroundOrder)
    p = np.poly1d(func)

    # Transform
    x = []
    newChromIntensity = []
    for i in data:
        x.append(i[0])
        newChromIntensity.append(int(i[1]-p(i[0])))
    newData = zip(x,newChromIntensity)

    # Return corrected data
    return newData
    
def batchChromCalibration(data, calFile):
    """ TODO
    """
    # Get calibration values
    try:
        refPeaks = []
        with open(calFile.get(),'r') as fr:
            for line in fr:
                line = line.rstrip("\n").split("\t")
                refPeaks.append((str(line[0]), float(line[1]), float(line[2])))
    except IOError:
        tkMessageBox.showinfo("File Error","The selected reference file could not be opened.")

    # Get observed times
    time, intensity = zip(*data)
    timePairs = []
    for i in refPeaks:
        low = bisect.bisect_left(time, i[1]-i[2])
        high = bisect.bisect_right(time, i[1]+i[2])
        max_value = max(intensity[low:high])
        max_index = intensity[low:high].index(max_value)
        timePairs.append((i[1],time[low+max_index]))

    # Calibration
    observedTime = []
    expectedTime = []
    for i in timePairs:
        expectedTime.append(float(i[0]))
        observedTime.append(float(i[1]))
    z = np.polyfit(observedTime,expectedTime,2)
    f = np.poly1d(z)
    calibratedData = []
    for i in data:
        newX = f(i[0])
        calibratedData.append((newX,i[1]))

    # Return calibrated data
    return calibratedData

def batchProcess(calFile, analFile):
    """ TODO
    """
    # Calibration
    filesGrabbed = []
    for files in CALIBRATION_FILETYPES:
        for file in glob.glob(files):
            filesGrabbed.append(file)
    for file in filesGrabbed:
        data = {'Data':functions.openChrom(file),'Name':file}
        data['Data'] = batchBaselineCorrection(data['Data'])
        data['Data'] = batchChromCalibration(data['Data'], calFile)
        data['Name'] = "calibrated_"+str(data['Name'])
        batchWriteData(data)

    # Integration
    filesGrabbed = []
    for files in INTEGRATION_FILETYPES:
        for file in glob.glob(files):
            filesGrabbed.append(file)
    for file in filesGrabbed:
        batchQuantifyChrom(data, analFile)

def batchQuantifyChrom(data, analFile):
    """ TODO
    """
    peaks = []
    with open(analFile.get(),'r') as fr:
        for line in fr:
            line = line.rstrip("\n").split("\t")
            peaks.append((str(line[0]), float(line[1]), float(line[2])))
    time, intensity = zip(*data['Data'])
    results = []
    for i in peaks:
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])
        peakArea = 0
        for index,j in enumerate(intensity[low:high]):
            try:
                peakArea += j * (time[low+index]-time[low+index-1])
            except IndexError:
                continue
        results.append((i[0], i[1], peakArea))
    data['Name'] = str(data['Name'].split('.')[0])+".raw"
    with open(data['Name'],'w') as fw:
        for i in results:
            fw.write(str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\n")
            
def batchWriteData(data):
    """ TODO
    """
    with open(data['Name'],'w') as fw:
        for i in data['Data']:
            fw.write(str(i[0])+"\t"+str(i[1])+"\n")
