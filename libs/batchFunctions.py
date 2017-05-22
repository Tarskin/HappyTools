#! /usr/bin/env python

# General imports
from datetime import datetime
#from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from Tkinter import StringVar, Toplevel, Label
import bisect
import glob
import math
import numpy as np
import operator
import os
import sys
import ttk
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

    # Baseline determination
    background = []
    chunks = [data[x:x+functions.points] for x in xrange(0, len(data), functions.points)]
    for i in chunks:
        buff1, buff2 = zip(*i)
        min_index, min_value = min(enumerate(buff2), key=operator.itemgetter(1))
        if buff1[0] > functions.start and buff1[-1] < functions.end:
            background.append((buff1[min_index], buff2[min_index]))
    time, intensity = zip(*background)
    newX = np.linspace(min(time), max(time),100)
    func = np.polyfit(time, intensity, functions.baselineOrder)
    p = np.poly1d(func)

    # Transform
    x = []
    newChromIntensity = []
    for i in data:
        x.append(i[0])
        newChromIntensity.append(int(i[1]-p(i[0])))
    
    # Uplift
    foo1, foo2 = zip(*data)
    start = bisect.bisect_left(foo1, functions.start)
    end = bisect.bisect_right(foo1, functions.end)
    offset = abs(min(min(newChromIntensity[start:end]),0))
    newData = zip(foo1,[x+offset for x in newChromIntensity])

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
                try:
                    refPeaks.append((str(line[0]), float(line[1]), float(line[2])))
                except ValueError:
                    pass
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
        newX = float(f(i[0]))
        calibratedData.append((newX,i[1]))

    # Return calibrated data
    return calibratedData

def batchProcess(calFile, analFile):
    """ TODO
    """
    start = datetime.now()

    # Progress bar
    calPerc = StringVar()
    intPerc = StringVar()
    calPerc.set("0%")
    intPerc.set("0%")
    barWindow = Toplevel()
    barWindow.title("Progress Bar")
    cal = Label(barWindow, text="Calibration", padx=25)
    cal.grid(row=0, column=0, sticky="W")
    ft = ttk.Frame(barWindow)
    ft.grid(row=1, columnspan=2, sticky="")
    perc1 = Label(barWindow, textvariable=calPerc)
    perc1.grid(row=0, column=1, padx=25)
    progressbar = ttk.Progressbar(ft, length=100, mode='determinate')
    progressbar.grid(row=1, columnspan=2, sticky="")
    ext = Label(barWindow, text="Integration", padx=25)
    ext.grid(row=2, column=0, sticky="W")
    ft2 = ttk.Frame(barWindow)
    ft2.grid(row=3, columnspan=2, sticky="")
    perc2 = Label(barWindow, textvariable=intPerc)
    perc2.grid(row=2, column=1, padx=25)
    progressbar2 = ttk.Progressbar(ft2, length=100, mode='determinate')
    progressbar2.grid(row=3, columnspan=2, sticky="")

    # Calibration   
    if calFile.get() != "":
        filesGrabbed = []
        for files in CALIBRATION_FILETYPES:
            for file in glob.glob(files):
                filesGrabbed.append(file)
        for index,file in enumerate(filesGrabbed):
            functions.updateProgressBar(progressbar, calPerc, index, len(filesGrabbed))
            try:
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now())+"\tCalibrating file: "+str(file)+"\n")
                data = {'Data':functions.openChrom(file),'Name':file}
                data['Data'] = batchBaselineCorrection(data['Data'])
                if 'blank' in data['Name'] or 'blanc' in data['Name']:
                    pass
                else:
                    data['Data'] = batchChromCalibration(data['Data'], calFile)
                data['Name'] = "calibrated_"+str(data['Name'])
                batchWriteData(data)
            except ValueError:
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now())+"\tIgnoring file: "+str(file)+" for calibration\n")
                pass
    functions.updateProgressBar(progressbar, calPerc, len(filesGrabbed), len(filesGrabbed))
    # Integration
    if analFile.get() != "":
        try:
            filesGrabbed = []
            for files in INTEGRATION_FILETYPES:
                for file in glob.glob(files):
                    filesGrabbed.append(file)
            for index,file in enumerate(filesGrabbed):
                functions.updateProgressBar(progressbar2, intPerc, index, len(filesGrabbed))
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now())+"\tQuantifying file: "+str(file)+"\n")
                data = {'Data':functions.openChrom(file),'Name':file}
                batchQuantifyChrom(data, analFile)
        except ValueError:
            if HappyTools.logging == True and HappyTools.logLevel >= 1:
                with open(HappyTools.logFile,'a') as fw:
                    fw.write(str(datetime.now())+"\tIgnoring a file: <UNKNOWN> for quantitation\n")
            pass    
        functions.updateProgressBar(progressbar2, intPerc, len(filesGrabbed), len(filesGrabbed))
        combineResults()
    end = datetime.now()
    tkMessageBox.showinfo("Status Message", "Batch Process finished on "+str(end)+" and took a total time of "+str(end-start))

def batchQuantifyChrom(data, analFile):
    """ TODO
    """
    peaks = []
    with open(analFile.get(),'r') as fr:
        for line in fr:
            line = line.rstrip("\n").split("\t")
            try:
                peaks.append((str(line[0]), float(line[1]), float(line[2])))
            except ValueError:
                if HappyTools.logging == True and HappyTools.logLevel > 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now())+"\tIgnoring line: "+str(line)+" from file: "+str(analFile)+"\n")
                pass
    time, intensity = zip(*data['Data'])
    results = []
    for i in peaks:
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])
        peakArea = 0
        backgroundArea = 0
        residual = 0
        signalNoise = "Nan"
        # Get signal-to-noise
        lowBackground = bisect.bisect_left(time,max(i[1]-functions.backgroundWindow,functions.start))
        highBackground = bisect.bisect_right(time,min(i[1]+functions.backgroundWindow,functions.end))
        backgroundData = intensity[lowBackground:low]+intensity[high:highBackground]
        if functions.backgroundNoiseMethod == "NOBAN":
            NOBAN = functions.noban(backgroundData)
        elif functions.backgroundNoiseMethod == "MT":
            NOBAN = functions.backgroundNoise(backgroundData)
        signalNoise = (max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise']
        # Get peak Area
        for index,j in enumerate(intensity[low:high]):
            try:
                peakArea += max(j,0) * (time[low+index]-time[low+index-1])
                backgroundArea += max(NOBAN['Background'],0) * (time[low+index]-time[low+index-1])
            except IndexError:
                continue
        peakNoise = np.std(intensity[low:high])
        #################################################
        # Gaussian fit on main points(to get residuals) #
        #################################################
        # 1. fit a spline through all data points
        # 2. get first derivative of spline function
        # 3. get both data points (left/right) where x' = 0 around the data point with max intensity
        # 4. fit a gaussian through the data points returned from 3
        # 4a. Plot the gaussian (stopping at baseline) in the pdf
        # 5. Calculate area under gaussian curve
        # 6. Calculate area under raw data points
        # 7. Calculate percentage of total area explained by gaussian area
        if HappyTools.logging == True and HappyTools.logLevel > 1:
            with open(HappyTools.logFile,'a') as fw:
                fw.write(str(datetime.now())+"\tIdentifying local minima and maxima, using first derivative\n")
        x_data = np.array(time[low:high])
        y_data = np.array(intensity[low:high])
        newX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        f = InterpolatedUnivariateSpline(x_data, y_data)
        fPrime = f.derivative()
        newY = f(newX)
        newPrimeY = fPrime(newX)
        maxm = argrelextrema(newPrimeY, np.greater)
        minm = argrelextrema(newPrimeY, np.less)
        breaks = maxm[0].tolist() + minm[0].tolist()
        breaks = sorted(breaks)
        maxPoint = 0
        try:
            if max(newY[0:breaks[0]]) > maxPoint:
                maxPoint = max(newY[0:breaks[0]])
                xData = newX[0:breaks[0]]
                yData = [x - NOBAN['Background'] for x in newY[0:breaks[0]]]
            for index,j in enumerate(breaks):
                try:
                    if max(newY[breaks[index]:breaks[index+1]]) > maxPoint:
                        maxPoint = max(newY[breaks[index]:breaks[index+1]])
                        xData = newX[breaks[index]:breaks[index+1]]
                        yData = [x - NOBAN['Background'] for x in newY[breaks[index]:breaks[index+1]]]
                except:
                    if max(newY[breaks[index]:-1]) > maxPoint:
                        maxPoint = max(newY[breaks[index]:-1])
                        xData = newX[breaks[index]:-1]
                        yData = [x - NOBAN['Background'] for x in newY[breaks[index]:-1]]
                    pass
        except IndexError:
            pass
        # Gaussian fit on main points
        newGaussX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        p0 = [np.max(yData), xData[np.argmax(yData)],0.1]
        try:
            coeff, var_matrix = curve_fit(functions.gaussFunction, xData, yData, p0)
            newGaussY = functions.gaussFunction(newGaussX, *coeff)
            newGaussY = [x + NOBAN['Background'] for x in newGaussY]
            totalArea = 0
            gaussArea = 0
            for index,j in enumerate(intensity[low:high]):
                try:
                    totalArea += max(j-NOBAN['Background'],0) * (time[low+index]-time[low+index-1])
                    gaussArea += max(functions.gaussFunction(time[low+index],*coeff),0) * (time[low+index]-time[low+index-1])
                except IndexError:
                    if HappyTools.logging == True and HappyTools.logLevel > 1:
                        with open(HappyTools.logFile,'a') as fw:
                            fw.write(str(datetime.now())+"\t<PLACEHOLDER 1>"+str(low)+" - "+str(high)+"\n")
                    continue
            residual = gaussArea / totalArea
            if functions.createFigure == "True":
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                       fw.write(str(datetime.now())+"\tCreating figure for analyte: "+str(i[0])+"\n")
                # Generate plot
                fig =  plt.figure()
                ax = fig.add_subplot(111)
                plt.plot(time[low:high], intensity[low:high], 'b*')
                plt.plot((newX[0],newX[-1]),(NOBAN['Background'],NOBAN['Background']),'red')
                plt.plot((newX[0],newX[-1]),(NOBAN['Background']+NOBAN['Noise'],NOBAN['Background']+NOBAN['Noise']),color='green')
                plt.plot(newX,newY, color='blue',linestyle='dashed')
                plt.plot(newGaussX, newGaussY, color='green',linestyle='dashed')
                plt.plot((time[intensity[low:high].index(max(intensity[low:high]))+low],time[intensity[low:high].index(max(intensity[low:high]))+low]),(NOBAN['Background'],max(intensity[low:high])),color='orange',linestyle='dotted')
                plt.legend(['Raw Data','Background','Noise','Univariate Spline','Gaussian Fit ('+str(int(residual*100))+'%)','Signal (S/N '+str(round((max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise'],1))+")"], loc='best')
                plt.title(str(data['Name']+"-"+str(i[0])))
                plt.xlabel("rt [m]")
                plt.ylabel("intensity [au]")
                plt.savefig(str(data['Name'])+"-"+str(i[0])+".pdf",bbox_inches="tight")
                plt.close(fig)
        except:
            if HappyTools.logging == True and HappyTools.logLevel > 1:
                with open(HappyTools.logFile,'a') as fw:
                    fw.write(str(datetime.now())+"\tUnable to determine residuals for peak: "+str(i[1])+"\n")
            residual = "Nan"
            pass
        results.append({'Peak':i[0], 'Time':i[1], 'Area':peakArea, 'PeakNoise':peakNoise, 'Residual':residual, 'S/N':signalNoise,'Background':NOBAN['Background'],'Noise':NOBAN['Noise'],'BackgroundArea':backgroundArea})
    data['Name'] = str(data['Name'].split('.')[0])+".raw"
    with open(data['Name'],'w') as fw:
        fw.write("Name\tTime\tPeak Area\tS/N\tBackground\tNoise\tGaussian Residual RMS\tPeak Noise\tBackground Area\n")
        for i in results:
            fw.write(str(i['Peak'])+"\t"+str(i['Time'])+"\t"+str(i['Area'])+"\t"+str(i['S/N'])+"\t"+str(i['Background'])+"\t"+str(i['Noise'])+"\t"+str(i['Residual'])+"\t"+str(i['PeakNoise'])+"\t"+str(i['BackgroundArea'])+"\n")

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
                Buffer.append({'Peak':str(chunks[0]),'Time':float(chunks[1]),'Area':float(chunks[2]),'S/N':float(chunks[3]),'Background':float(chunks[4]),'Noise':float(chunks[5]),'Residual':float(chunks[6]),'PeakNoise':float(chunks[7]),'BackgroundArea':float(chunks[8])})
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

    # Write results, settings and version information
    with open(filename,'w') as fw:
        # Metadata
        fw.write("HappyTools Settings\n")
        fw.write("Version:\t"+str(HappyTools.version)+"\n")
        fw.write("Build:\t"+str(HappyTools.build)+"\n")
        fw.write("Start Time:\t"+str(functions.start)+"\n")
        fw.write("End Time:\t"+str(functions.end)+"\n")
        fw.write("Baseline Order:\t"+str(functions.baselineOrder)+"\n")
        fw.write("Background Window:\t"+str(functions.backgroundWindow)+"\n")
        fw.write("Background and noise method:\t"+str(functions.backgroundNoiseMethod)+"\n")
        if functions.backgroundNoiseMethod == "MT":
            fw.write("MT Slice Points:\t"+str(functions.slicepoints)+"\n")
        elif functions.backgroundNoiseMethod == "NOBAN":
            fw.write("NOBAN Initial Estimate:\t"+str(functions.nobanStart)+"\n")
        fw.write("Noise:\t"+str(functions.noise)+"\n")
        fw.write("\n")

        # Area (non background subtracted)
        fw.write("Peak Area")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['Area']))
            fw.write("\n")
        fw.write("\n")

        # Area (Background subtracted)
        fw.write("Peak Area (Background Subtracted)")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(max(j['Area']-j['BackgroundArea'],0)))
            fw.write("\n")
        fw.write("\n")

        # Relative Area
        fw.write("Relative Peak Area (TAN)")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            total = 0.
            for j in i['Data']:
                total += j['Area']
            for j in i['Data']:
                try:
                    fw.write("\t"+str(j['Area']/total))
                except ZeroDivisionError:
                    fw.write("\t"+str(0.0))
            fw.write("\n")
        fw.write("\n")

        # Relative Area (Background subtracted)
        fw.write("Relative Peak Area (TAN, Background Subtracted)")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            total = 0.
            for j in i['Data']:
                total += max(j['Area']-j['BackgroundArea'],0)
            for j in i['Data']:
                try:
                    fw.write("\t"+str(max(j['Area']-j['BackgroundArea'],0)/total))
                except ZeroDivisionError:
                    fw.write("\t"+str(0.0))
            fw.write("\n")
        fw.write("\n")

        # Peak Noise (standard deviation of the integration window)
        fw.write("Peak Noise (standard deviation of integration window)")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            total = 0.
            for j in i['Data']:
                fw.write("\t"+str(j['PeakNoise']))
            fw.write("\n")
        fw.write("\n")
        
        # Background
        fw.write("Background")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['Background']))
            fw.write("\n")
        fw.write("\n")

        # Noise
        fw.write("Noise")
        fw.write(header)
        for i in Results:
            fw.write(i['File'])
            for j in i['Data']:
                fw.write("\t"+str(j['Noise']))
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
