#! /usr/bin/env python

# General imports
from datetime import datetime
from scipy.optimize import curve_fit
import bisect
import glob
import math
import numpy as np
import operator
import os
import sys
import tkMessageBox

# Test imports
import matplotlib.pyplot as plt
# Custom libraries
sys.path.append('..')
import functions
import HappyTools

# Defines
CALIBRATION_FILETYPES = ["*.txt","*.arw"]
INTEGRATION_FILETYPES = ["calibrated*.txt"]

# Functions
def batchBaselineCorrection(data):
    """ TODO
    """

    # Background determination
    background = []
    chunks = [data[x:x+HappyTools.points] for x in xrange(0, len(data), HappyTools.points)]
    for i in chunks:
        buff1, buff2 = zip(*i)
        min_index, min_value = min(enumerate(buff2), key=operator.itemgetter(1))
        if buff1[0] > HappyTools.start and buff1[-1] < HappyTools.end:
            background.append((buff1[min_index], buff2[min_index]))
    time, intensity = zip(*background)
    newX = np.linspace(min(time), max(time),100)
    func = np.polyfit(time, intensity, HappyTools.backgroundOrder)
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
        data = {'Data':functions.openChrom(file),'Name':file}
        batchQuantifyChrom(data, analFile)
    combineResults()
    tkMessageBox.showinfo("Status Message", "Batch Process finished on "+str(datetime.now()))

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
        residual = 0
        signalNoise = "Nan"
        # Get peak Area
        for index,j in enumerate(intensity[low:high]):
            try:
                peakArea += j * (time[low+index]-time[low+index-1])
            except IndexError:
                continue
        # Get signal-to-noise
        lowBackground = bisect.bisect_left(time,max(i[1]-HappyTools.backgroundWindow,HappyTools.start))
        highBackground = bisect.bisect_right(time,min(i[1]+HappyTools.backgroundWindow,HappyTools.end))
        backgroundData = intensity[lowBackground:low]+intensity[high:highBackground]
        NOBAN = noban(backgroundData)
        signalNoise = (max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise']
        # Gaussian fit (to get residuals)
        x_data = np.array(time[low:high])
        y_data = np.array(intensity[low:high])
        newX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        p0 = [np.max(y_data), x_data[np.argmax(y_data)],0.1]
        try:
            coeff, var_matrix = curve_fit(gaussFunction, x_data, y_data, p0)
            newY = gaussFunction(newX, *coeff)
            # Generate plot
            fig =  plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(time[low:high], intensity[low:high], 'b*')
            plt.plot(newX,newY, 'b--')
            plt.legend(['Raw Data','Gaussian Fit'], loc='best')
            plt.savefig(str(data['Name'])+"-"+str(i[0])+".pdf",bbox_inches="tight")
            plt.close(fig)
            # Get residuals
            for index,j in enumerate(time[low:high]):
                residual += abs(intensity[index]-gaussFunction(j, *coeff))**2
            residual = math.sqrt(residual)
        except RuntimeError:
            residual = "Nan"
        results.append({'Peak':i[0], 'Time':i[1], 'Area':peakArea, 'Residual':residual, 'S/N':signalNoise})
    data['Name'] = str(data['Name'].split('.')[0])+".raw"
    with open(data['Name'],'w') as fw:
        fw.write("Name\tTime\tPeak Area\tS/N\tGaussian Residual RMS\n")
        for i in results:
            fw.write(str(i['Peak'])+"\t"+str(i['Time'])+"\t"+str(i['Area'])+"\t"+str(i['S/N'])+"\t"+str(i['Residual'])+"\n")
            
def batchWriteData(data):
    """ TODO
    """
    with open(data['Name'],'w') as fw:
        for i in data['Data']:
            fw.write(str(i[0])+"\t"+str(i[1])+"\n")

def combineResults():
    """ TODO
    """
    # Read the raw files and construct a data structure
    Results = []
    for file in glob.glob("*.raw"):
        Buffer = []
        with open(file,'r') as fr:
            fr.readline()
            for line in fr:
                chunks = line.rstrip('\n').split('\t')
                Buffer.append({'Peak':str(chunks[0]),'Time':float(chunks[1]),'Area':float(chunks[2]),'S/N':float(chunks[3]),'Residual':float(chunks[4])})
            Results.append({'File':str(file),'Data':Buffer})

    # Construct the filename for the output
    utc_datetime = datetime.utcnow()
    s = utc_datetime.strftime("%Y-%m-%d-%H%MZ")
    filename = s + "_" + HappyTools.output

    # Construct header
    header = ""
    for i in Results:
        for j in i['Data']:
            header = header + "\t"+str(j['Peak'])
        header = header + "\n"
        for j in i['Data']:
            header = header + "\t"+str(j['Time'])
        header = header + "\n"
        break

    # Write results
    with open(filename,'w') as fw:
        # Area
        fw.write("Peak Area")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['Area']))
            fw.write("\n")
        fw.write("\n")

        # S/N
        fw.write("Signal-to-Noise")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['S/N']))
            fw.write("\n")
        fw.write("\n")

        # Residuals
        fw.write("Residuals")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['Residual']))
            fw.write("\n")
        fw.write("\n")

def gaussFunction(x, *p):
    """ TODO
    """
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
def noban(data):
    """ NOBAN implementation based on Jansen et al, 2016.
    """
    sortedData = sorted(data)
    startSize = int(0.25 * float(len(sortedData)))
    currSize = startSize
    currAverage = np.average(sortedData[0:currSize])
    if HappyTools.noise == "MM":
        currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
    elif HappyTools.noise == "RMS":
        currNoise = np.std(sortedData[0:currSize])
    directionFlag = 0
    for k in range(0,len(sortedData)-(startSize+1)):
        if sortedData[currSize+1] < currAverage + 3 * currNoise:
            directionFlag == 1
            currSize += 1
            currAverage =  np.average(sortedData[0:currSize])
            if HappyTools.noise == "MM":
                currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
            elif HappyTools.noise == "RMS":
                currNoise = np.std(sortedData[0:currSize])
        else:
            if sortedData[currSize-1] > currAverage + 3 * currNoise and directionFlag == 0:
                currSize -= 1
                currAverage = np.average(sortedData[0:currSize])
                if HappyTools.noise == "MM":
                    currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                elif HappyTools.noise == "RMS":
                    currNoise = np.std(sortedData[0:currSize])
            else:
                break
    background = currAverage
    noise = currNoise
    return {'Background': background, 'Noise': noise}
