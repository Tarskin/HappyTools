#! /usr/bin/env python

# General imports
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import PchipInterpolator, Akima1DInterpolator
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
EXCLUSION_FILES = ["LICENSE.txt","CHANGELOG.txt"]
CALIBRATION_FILETYPES = ["*.txt","*.arw"]
INTEGRATION_FILETYPES = ["calibrated*.txt"]

# Functions
def baselineCorrection(data):
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
    x = [a for a,b in data]
    newChromIntensity = [b-p(a) for a,b in data]

    # Uplift
    foo1, foo2 = zip(*data)
    start = bisect.bisect_left(foo1, functions.start)
    end = bisect.bisect_right(foo1, functions.end)
    offset = abs(min(min(newChromIntensity[start:end]),0))
    newData = zip(foo1,[x+offset for x in newChromIntensity])

    # Return corrected data
    return newData



def batchCalibrationControl(data, calFile):
    """ TODO 
    """
    # Get calibration values
    refPeaks = functions.getPeakList(calFile.get())

    # Get observed times
    timePairs = determineTimepairs(refPeaks, data)

    # Check number calibrants
    calibratedData = performCalibration(timePairs, data)

    # Return calibrated data
    return calibratedData

def batchProcess(calFile, analFile, batchFolder):
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
            for file in glob.glob(os.path.join(batchFolder.get(),files)):
                if file not in EXCLUSION_FILES:
                    filesGrabbed.append(os.path.join(batchFolder.get(),file))
        for index,file in enumerate(filesGrabbed):
            functions.updateProgressBar(progressbar, calPerc, index, len(filesGrabbed))
            try:
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now().replace(microsecond=0))+"\tCalibrating file: "+str(file)+"\n")
                data = {'Data':functions.openChrom(file),'Name':file}
                data['Data'] = baselineCorrection(data['Data'])
                # This is some custom code to quantify blancs (if desired)
                #if 'blank' in data['Name'].lower() or 'blanc' in data['Name'].lower():
                    #pass
                #else:
                    #data['Data'] = batchChromCalibration(data['Data'], calFile)
                # Must be a better way of doing this
                foo = batchCalibrationControl(data['Data'], calFile)
                data['Data'] = foo['Data']
                data['Function'] = foo['Function']
                if data['Data'] == None:
                    continue
                data['Name'] = os.path.join(batchFolder.get(), "calibrated_"+os.path.basename(data['Name']))
                writeData(batchFolder, data)
            except ValueError:
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now().replace(microsecond=0))+"\tIgnoring file: "+str(file)+" for calibration\n")
                pass
    functions.updateProgressBar(progressbar, calPerc, 1, 1)
    # Integration
    if analFile.get() != "":
        try:
            filesGrabbed = []
            for files in INTEGRATION_FILETYPES:
                for file in glob.glob(os.path.join(batchFolder.get(),files)):
                    if file not in EXCLUSION_FILES:
                        filesGrabbed.append(os.path.join(batchFolder.get(),file))
            for index,file in enumerate(filesGrabbed):
                functions.updateProgressBar(progressbar2, intPerc, index, len(filesGrabbed))
                if HappyTools.logging == True and HappyTools.logLevel >= 1:
                    with open(HappyTools.logFile,'a') as fw:
                        fw.write(str(datetime.now().replace(microsecond=0))+"\tQuantifying file: "+str(file)+"\n")
                data = {'Data':functions.openChrom(file),'Name':file}
                batchQuantitationControl(data, analFile, batchFolder)
        except ValueError:
            if HappyTools.logging == True and HappyTools.logLevel >= 1:
                with open(HappyTools.logFile,'a') as fw:
                    fw.write(str(datetime.now().replace(microsecond=0))+"\tIgnoring file: "+str(file)+" for quantitation. " \
                        "The 'Start' or 'End' parameter do not match the specified analytes.\n")
            pass    
        functions.updateProgressBar(progressbar2, intPerc, 1, 1)

        if HappyTools.logging == True and HappyTools.logLevel >= 1:
            with open(HappyTools.logFile,'a') as fw:
               fw.write(str(datetime.now().replace(microsecond=0))+"\tCreating summary file\n")
        combineResults(batchFolder)
    end = datetime.now()
    tkMessageBox.showinfo("Status Message", "Batch Process finished on "+str(end)+" and took a total time of "+str(end-start))

def batchQuantitationControl(data, analFile, batchFolder):
    """Quantify the current chromatogram and write results to disk.

    This function will open the analyte file (analFile), read all lines
    and split the line on tabs. The individual segments (name, time and
    time window) are then appended as a tuple to the list peaks.
    Next, the function will iterate over all tuples in the list peaks 
    and isolate the relevant segment of the chromatogram using a binary
    search. The local background and noise is then determined using 
    either the NOBAN or MT method, prior to integrating the peak and
    background areas. The best fitting Gaussian (for the highest
    intensity datapoints) is determined and used to calculate the
    overlap between the Gaussian and observed pattern. Optionally, a
    figure is created showing the raw data, fitted Gaussian peak,
    background, noise and the overlap percentage, which is saved to the
    disk. Lastly, the function writes all results to the disk in a raw
    file.

    Keyword arguments:
    data -- list of (time,intensity) tuples
    analFile -- unicode string
    """
    peaks = functions.getPeakList(analFile.get())
    time, intensity = zip(*data['Data'])
    results = []

    # Plot chromatogram region of interest (check if X[0] and X[-1] can be found before start)
    if functions.createFigure == "True" and bisect.bisect_left(time,functions.start) and bisect.bisect_right(time,functions.end):
        pdf = PdfPages(os.path.join(batchFolder.get(),os.path.splitext(os.path.basename(data['Name']))[0]+".pdf"))
        plotOverview(pdf, peaks, data, time, intensity)

    for i in peaks:
        # Initialize values
        peakArea = 0
        backgroundArea = 0
        totalArea = 0
        gaussArea = 0
        height = 0
        signalNoise = "NAN"
        residual = "NAN"
        fwhm = {'fwhm':0, 'width':0, 'center':0}

        # Get time boundaries
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])

        # Get signal-to-noise
        lowBackground = bisect.bisect_left(time,max(i[1]-functions.backgroundWindow,functions.start))
        highBackground = bisect.bisect_right(time,min(i[1]+functions.backgroundWindow,functions.end))
        backgroundData = intensity[lowBackground:highBackground]
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

        # Get breakpoints (where f'(x) == 0)
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

        # Initialize maxPoint, xData and yData
        maxPoint = 0
        xData = newX
        yData = [x - NOBAN['Background'] for x in newY]

        # Subset the data
        # Region from newY[0] to breaks[0]
        try:
            if max(newY[0:breaks[0]]) > maxPoint:
                maxPoint = max(newY[0:breaks[0]])
                xData = newX[0:breaks[0]]
                yData = [x - NOBAN['Background'] for x in newY[0:breaks[0]]]
        except IndexError:
            pass
        # Regions between breaks[x] and breaks[x+1]
        try:
            for index,j in enumerate(breaks):
                if max(newY[breaks[index]:breaks[index+1]]) > maxPoint:
                    maxPoint = max(newY[breaks[index]:breaks[index+1]])
                    xData = newX[breaks[index]:breaks[index+1]]
                    yData = [x - max(NOBAN['Background'],0) for x in newY[breaks[index]:breaks[index+1]]]
        except IndexError:
            pass
        # Region from break[-1] to newY[-1]
        try:
            if max(newY[breaks[-1]:-1]) > maxPoint:
                maxPoint = max(newY[breaks[-1]:-1])
                xData = newX[breaks[-1]:-1]
                yData = [x - NOBAN['Background'] for x in newY[breaks[-1]:-1]]
        except IndexError:
            pass

        # Gaussian fit on main points
        peak = xData[yData > np.exp(-0.5)*max(yData)]
        guess_sigma = 0.5*(max(peak) - min(peak))
        newGaussX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        p0 = [np.max(yData), xData[np.argmax(yData)],guess_sigma]
        try:
            coeff, var_matrix = curve_fit(functions.gaussFunction, xData, yData, p0)
            newGaussY = functions.gaussFunction(newGaussX, *coeff)
            newGaussY = [x + NOBAN['Background'] for x in newGaussY]
            for index,j in enumerate(intensity[low:high]):
                gaussArea += max(functions.gaussFunction(time[low+index],*coeff),0) * (time[low+index]-time[low+index-1])
            fwhm = functions.fwhm(coeff)
            height = functions.gaussFunction(fwhm['center']+fwhm['width'], *coeff)+NOBAN['Background']
        except TypeError:
            if HappyTools.logging == True and HappyTools.logLevel > 1:
                with open(HappyTools.logFile,'a') as fw:
                    fw.write(str(datetime.now().replace(microsecond=0))+"\tNot enough data points to fit a Gaussian to peak: "+str(i[0])+"\n")
        except RuntimeError:
            if HappyTools.logging == True and HappyTools.logLevel > 1:
                with open(HappyTools.logFile,'a') as fw:
                    fw.write(str(datetime.now().replace(microsecond=0))+"\tUnable to determine residuals for peak: "+str(i[1])+"\n")

        # Determine Area
        for index,j in enumerate(intensity[low:high]):
            totalArea += max(j-NOBAN['Background'],0) * (time[low+index]-time[low+index-1])

        # Determine Residual
        try:
            if gaussArea != 0:
                residual = min(gaussArea / totalArea, 1.0)
        except ZeroDivisionError:
            pass

        # Generate plot
        if functions.createFigure == "True" and residual != "NAN":
            if HappyTools.logging == True and HappyTools.logLevel >= 1:
                with open(HappyTools.logFile,'a') as fw:
                   fw.write(str(datetime.now().replace(microsecond=0))+"\tCreating figure for analyte: "+str(i[0])+"\n")
            details = {'fwhm':fwhm, 'height':height, 'NOBAN':NOBAN, 'newData':zip(newX,newY), 'newGauss':zip(newGaussX,newGaussY), 
                    'data':zip(time,intensity), 'low':low, 'high':high, 'residual':residual, 'i':i}
            plotIndividual(pdf, details)

        results.append({'Peak':i[0], 'Time':i[1], 'Area':peakArea, 'PeakNoise':peakNoise, 'Residual':residual, 'S/N':signalNoise,
                'Background':NOBAN['Background'],'Noise':NOBAN['Noise'],'BackgroundArea':backgroundArea, 'fwhm':fwhm['fwhm'],
                'ActualTime':fwhm['center']})
    if functions.createFigure == "True":
        pdf.close()

    # Write results to disk
    data['Name'] = os.path.splitext(os.path.basename(data['Name']))[0]+".raw"
    with open(data['Name'],'w') as fw:
        fw.write("Name\tTime\tPeak Area\tS/N\tBackground\tNoise\tGaussian Residual RMS\tPeak Noise\tBackground Area\tPeak Time\tFWHM\n")
        for i in results:
            fw.write(str(i['Peak'])+"\t"+str(i['Time'])+"\t"+str(i['Area'])+"\t"+str(i['S/N'])+"\t"+str(i['Background'])+"\t"+
                    str(i['Noise'])+"\t"+str(i['Residual'])+"\t"+str(i['PeakNoise'])+"\t"+str(i['BackgroundArea'])+"\t"+
                    str(i['ActualTime'])+"\t"+str(i['fwhm'])+"\n")

def combineResults(batchFolder):
    """ TODO
    """
    # Read the raw files and construct a data structure
    Results = []
    for file in glob.glob(os.path.join(batchFolder.get(),"*.raw")):
        Buffer = []
        with open(file,'r') as fr:
            fr.readline()
            for line in fr:
                chunks = line.rstrip('\n').split('\t')
                Buffer.append({'Peak':str(chunks[0]),'Time':float(chunks[1]),'Area':float(chunks[2]),'S/N':float(chunks[3]),
                            'Background':float(chunks[4]),'Noise':float(chunks[5]),'Residual':float(chunks[6]),'PeakNoise':float(chunks[7]),
                            'BackgroundArea':float(chunks[8]),'ActualTime':float(chunks[9]),'fwhm':float(chunks[10])})
        with open(os.path.splitext(os.path.basename(file))[0]+".cal") as fr:
            formula = fr.readline()
        Results.append({'File':str(os.path.splitext(os.path.basename(file))[0]), 'Calibration':str(formula), 'Data':Buffer})

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
    with open(os.path.join(batchFolder.get(),filename),'w') as fw:
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
        if functions.absInt.get() == 1 and functions.bckSub.get() == 0:
            fw.write("Peak Area")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['Area']))
                fw.write("\n")
            fw.write("\n")

        # Area (Background subtracted)
        if functions.absInt.get() == 1 and functions.bckSub.get() == 1:
            fw.write("Peak Area (Background Subtracted)")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(max(j['Area']-j['BackgroundArea'],0)))
                fw.write("\n")
            fw.write("\n")

        # Relative Area
        if functions.relInt.get() == 1 and functions.bckSub.get() == 0:
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
        if functions.relInt.get() == 1 and functions.bckSub.get() == 1:
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
        if functions.bckNoise.get() == 1:
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
        if functions.bckNoise.get() == 1:
            fw.write("Background")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['Background']))
                fw.write("\n")
            fw.write("\n")

        # Noise
        if functions.bckNoise.get() == 1:
            fw.write("Noise")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['Noise']))
                fw.write("\n")
            fw.write("\n")

        # S/N
        if functions.peakQual.get() == 1:
            fw.write("Signal-to-Noise")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['S/N']))
                fw.write("\n")
            fw.write("\n")

        # GPQ
        if functions.peakQual.get() == 1:
            fw.write("GPQ (Gaussian Peak Quality)")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['Residual']))
                fw.write("\n")
            fw.write("\n")

        # FWHM
        if functions.peakQual.get() == 1:
            fw.write("FWHM")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(j['fwhm']))
                fw.write("\n")
            fw.write("\n")

        # Tr residual
        if functions.peakQual.get() == 1:
            fw.write("Retention Time Residual")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                if i['Calibration']:
                    fw.write(" ["+str(i['Calibration'])+"]")
                for j in i['Data']:
                    residualTime = abs(float(j['ActualTime']) - float(j['Time']))
                    fw.write("\t"+str(residualTime))
                fw.write("\n")
            fw.write("\n")

        # Peak Tr
        if functions.peakQual.get() == 1:
            fw.write("Retention Time")
            fw.write(header)
            for i in Results:
                fw.write(i['File'])
                for j in i['Data']:
                    fw.write("\t"+str(float(j['ActualTime'])))
                fw.write("\n")
            fw.write("\n")

def determineTimepairs(refPeaks, data):
    """ TODO
    """
    time, intensity = zip(*data)
    timePairs = []
    for i in refPeaks:
        low = bisect.bisect_left(time, i[1]-i[2])
        high = bisect.bisect_right(time, i[1]+i[2])
        lowBackground = bisect.bisect_left(time, max(i[1]-functions.backgroundWindow,functions.start))
        highBackground = bisect.bisect_right(time, min(i[1]+functions.backgroundWindow,functions.end))
        if functions.backgroundNoiseMethod == "NOBAN":
            NOBAN = functions.noban(intensity[lowBackground:highBackground])
        elif functions.backgroundNoiseMethod == "MT":
            NOBAN = functions.backgroundNoise(intensity[lowBackground:highBackground])
        max_value = max(intensity[low:high])
        max_index = intensity[low:high].index(max_value)
        if ((max_value - NOBAN['Background'])/NOBAN['Noise']) >= functions.minPeakSN:
            timePairs.append((i[1],time[low+max_index]))
    return timePairs

def performCalibration(timePairs, data):
    """ TODO
    """
    time, intensity = zip(*data)
    try:
        if len(timePairs) >= functions.minPeaks:
            expectedTime, observedTime = zip(*timePairs)
            #z = np.polyfit(observedTime,expectedTime,2)
            #f = np.poly1d(z)
            f = functions.ultraPerformanceCalibration(observedTime,expectedTime,time[0], time[-1])
            calibratedData = zip(f(time),intensity)
        else:
            calibratedData = None
            if HappyTools.logging == True and HappyTools.logLevel >= 1:
                with open(HappyTools.logFile,'a') as fw:
                     fw.write(str(datetime.now().replace(microsecond=0))+"\tFile not calibrated due to lack of features, "+
                            str(len(timePairs))+" passed the minimum S/N ("+str(functions.minPeakSN)+") while "+str(functions.minPeaks)+
                            " were needed\n")
    except NameError:
        calibratedData = None
        if HappyTools.logging == True and HappyTools.logLevel >= 1:
            with open(HappyTools.logFile,'a') as fw:
                 fw.write(str(datetime.now().replace(microsecond=0))+"\tFile not calibrated due to lack of features, "+
                        str(len(timePairs))+" passed the minimum S/N ("+str(functions.minPeakSN)+") while "+
                        str(functions.minPeaks)+" were needed\n")  
    return {"Data":calibratedData, "Function":f}

def plotIndividual(pdf, details):
    """ TODO
    """
    # Unpack details
    low = details['low']
    high = details['high']
    fwhm = details['fwhm']
    NOBAN = details['NOBAN']
    height = details['height']
    residual = details['residual']
    i = details['i']
    newX, newY = zip(*details['newData'])
    newGaussX, newGaussY = zip(*details['newGauss'])
    time, intensity = zip(*details['data'])

    # Plot
    fig =  plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.plot(time[low:high], intensity[low:high], 'b*')
    plt.plot((newX[0],newX[-1]),(NOBAN['Background'],NOBAN['Background']),'red')
    plt.plot((newX[0],newX[-1]),(NOBAN['Background']+NOBAN['Noise'],NOBAN['Background']+NOBAN['Noise']),color='green')
    plt.plot(newX,newY, color='blue',linestyle='dashed')
    plt.plot(newGaussX, newGaussY, color='green',linestyle='dashed')
    plt.plot((time[intensity[low:high].index(max(intensity[low:high]))+low],time[intensity[low:high].index(max(intensity[low:high]))+low]),
            (NOBAN['Background'],max(intensity[low:high])),color='orange',linestyle='dotted')
    plt.plot((min(max(fwhm['center']-fwhm['width'],newX[0]),newX[-1]),max(min(fwhm['center']+fwhm['width'],newX[-1]),newX[0])),
            (height,height),color='red',linestyle='dashed')
    plt.legend(['Raw Data','Background','Noise','Univariate Spline','Gaussian Fit ('+str(int(residual*100))+
            '%)','Signal (S/N '+str(round((max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise'],1))+")",
            "FWHM:"+"{0:.2f}".format(fwhm['fwhm'])], loc='best')
    plt.title("Detail view: "+str(i[0]))
    plt.xlabel("Retention Time [m]")
    plt.ylabel("Intensity [au]")
    pdf.savefig(fig)
    plt.close(fig)

def plotOverview(pdf, peaks, data, time, intensity):
    """ TODO
    """
    d = pdf.infodict()
    d['Title'] = 'PDF Report for: '+str(os.path.splitext(os.path.basename(data['Name']))[0])
    d['Author'] = 'HappyTools version: '+str(HappyTools.version)+" build: "+str(HappyTools.build)
    d['CreationDate'] = datetime.now()
    low = bisect.bisect_left(time,functions.start)
    high = bisect.bisect_right(time,functions.end)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.plot(time[low:high], intensity[low:high], 'b-')
    plt.legend(['Raw Data'], loc='best')
    plt.title(str(os.path.splitext(os.path.basename(data['Name']))[0]))
    plt.xlabel("Retention Time [m]")
    plt.ylabel("Intensity [au]")
    for i in peaks:
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])
        newTime = np.linspace(time[low], time[high], len(time[low:high]))
        f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])
        newIntensity = f(newTime)
        ax.fill_between(time[low:high], 0, newIntensity, alpha=0.5)
        ax.text(i[1], max(intensity[low:high]), i[0])
    pdf.savefig(fig)
    plt.close(fig)

def writeData(batchFolder, data):
    """ TODO
    """
    with open(os.path.join(batchFolder.get(),os.path.splitext(data['Name'])[0]+'.cal'),'w') as fw:
        if isinstance(data['Function'],functions.powerLawCall):
            formula = data['Function'].describe()
        elif isinstance(data['Function'],np.lib.polynomial.poly1d):
            formula = ""
            for index,i in enumerate(data['Function']):
                if index < len(data['Function']):
                    formula += "{0:.2e}".format(i)+"x^"+str(len(data['Function'])-index)+" + "
                else:
                    formula += "{0:.2e}".format(i)
        elif isinstance(data['Function'],Akima1DInterpolator):
            formula = "Akima 1D Interpolation"
        elif isinstance(data['Function'],PchipInterpolator):
            formula = "Monotonic Piecewise Cubic Hermite Interpolating Polynomial"
        else:
            formula = "Unknown"
        fw.write(formula)
    with open(os.path.join(batchFolder.get(),data['Name']),'w') as fw:
        for i in data['Data']:
            fw.write(str(format(i[0],'0.'+str(functions.decimalNumbers)+'f'))+"\t"+str(format(i[1],'0.'+str(functions.decimalNumbers)+'f'))+"\n")
