#! /usr/bin/env python

# General imports
from matplotlib.pyplot import gca
from scipy.signal import savgol_filter
from Tkinter import *
import bisect
import glob
import numpy as np
import operator
import os
import sys
import tkFileDialog
import tkMessageBox

# Custom libraries
sys.path.append('..')
import batchFunctions
import HappyTools

# Functions
def addFile(fig,canvas):
    """ TODO
    """
    data = readData()
    file_path = tkFileDialog.askopenfilename()
    if not file_path:
        pass
    else:
        data.append((file_path,openChrom(file_path)))
    fig.clear()
    axes = fig.add_subplot(111)
    for i in data:
        x_array, y_array = zip(*i[1])
        axes.plot(x_array,y_array,label=str(os.path.split(i[0])[-1]))
    axes.legend()
    canvas.draw()

def baselineCorrection(fig,canvas):
    """ TODO
    """
    data = readData()

    # Background determination
    background = []
    chunks = [data[0][1][x:x+HappyTools.points] for x in xrange(0, len(data[0][1]), HappyTools.points)]
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
    for i in data[0][1]:
        x.append(i[0])
        newChromIntensity.append(int(i[1]-p(i[0])))
    newData = zip(x,newChromIntensity)

    # Plot & Write Data to Disk  
    multiData = [(os.path.split(data[0][0])[-1], data[0][1]),(os.path.split(data[0][0])[-1]+" (BC)",newData)]
    plotMultiData(fig,canvas,multiData)
    writeData(newData,os.path.split(data[0][0])[-1]+" (BC)")

def plotMultiData(fig,canvas,data):
    """ TODO
    """
    fig.clear()
    axes = fig.add_subplot(111)
    for i in data:
        xd = []
        yd = []
        for j in i[1]:
            xd.append(j[0])
            yd.append(j[1])
        axes.plot(xd,yd,label=os.path.split(i[0])[-1])
    axes.set_xlabel("Time [m]")
    axes.set_ylabel("Intensity [au]")
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles,labels)
    canvas.draw()

def batchPlot(fig,canvas):
    """ TODO
    """
    folder_path = tkFileDialog.askdirectory()
    data = []
    for file in glob.glob(str(os.path.join(folder_path,"*.txt"))):
        data.append((str(file),openChrom(file)))
    fig.clear()
    axes = fig.add_subplot(111)
    for i in data:
        x_array, y_array = zip(*i[1])
        axes.plot(x_array,y_array,label=str(os.path.split(i[0])[-1]))
    axes.legend()
    canvas.draw()

def batchPopup():
    """ TODO
    """

    calFile = StringVar()
    analFile = StringVar()

    def close():
        """ This function closes the settings popup and applies
        all the entered values to the parameters.
        """
        top.destroy()

    def calibrationFile():
        """ TODO
        """
        calFile.set(tkFileDialog.askopenfilename(title="Calibration File"))

    def analyteFile():
        """ TODO
        """
        analFile.set(tkFileDialog.askopenfilename(title="Analyte File"))

    def run():
        """ TODO
        """
        batchFunctions.batchProcess(calFile, analFile)

    top = Tk.top = Toplevel()
    top.protocol("WM_DELETE_WINDOW", lambda: close())
    calibrationButton = Button(top, text="Calibration File", command=lambda: calibrationFile())
    calibrationButton.grid(row=1, column=0, sticky=W)
    calibrationLabel = Label(top, textvariable=calFile)
    calibrationLabel.grid(row=1, column=1, sticky=W)

    analyteButton = Button(top, text="Analyte File", command=lambda: analyteFile())
    analyteButton.grid(row=2, column=0, sticky=W)
    analyteLabel = Label(top, textvariable=analFile)
    analyteLabel.grid(row=2, column=1, sticky=W)

    runButton = Button(top, text="Run", command=lambda: run())
    runButton.grid(row=3, column=0, sticky=W)
    closeButton = Button(top, text="Close", command=lambda: close())
    closeButton.grid(row=3, column=1, sticky=E)

def chromCalibration(fig,canvas):
    """ TODO
    """
    refFile = tkFileDialog.askopenfilename(title="Reference File")
    try:
        refPeaks = []
        with open(refFile,'r') as fr:
            for line in fr:
                line = line.rstrip("\n").split("\t")
                refPeaks.append((str(line[0]), float(line[1]), float(line[2])))
    except IOError:
        tkMessageBox.showinfo("File Error","The selected reference file could not be opened.")

    print "Rawr"
    # Get observed times
    data = readData()
    time, intensity = zip(*data[0][1])
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
    for i in data[0][1]:
        newX = f(i[0])
        calibratedData.append((newX,i[1]))

    # Plot & Write Data to Disk  
    multiData = [(os.path.split(data[0][0])[-1], data[0][1]),(os.path.split(data[0][0])[-1]+" (Cal)",calibratedData)]
    plotMultiData(fig,canvas,multiData)
    writeData(calibratedData,os.path.split(data[0][0])[-1]+" (Cal)")

def fileCleanup():
    """
    """
    search = os.path.join(str(os.getcwd()),"temp","*.*")
    files = glob.glob(search)
    for f in files:
        try:
            os.remove(f)
        except OSError:
            tkMessageBox.showinfo("File Error", "A temporary file is open in another program. Please exit HappyTools and close all temporary files before running HappyTools.")

def openChrom(file):
    """ TODO
    """
    with open(file,'r') as fr:
        chromData = []
        if 'txt' in file:
            for line in fr:
                if line[0].isdigit() == True:
                    lineChunks = line.rstrip().split()
                    # Check for non american thousand/decimal seperators
                    if ',' in lineChunks[0]:
                        time = lineChunks[0].replace('.','')
                        time = time.replace(',','.')
                        intensity = lineChunks[-1].replace('.','')
                        intensity = intensity.replace(',','.')
                        chromData.append((float(time),float(intensity)))
                    else:
                        chromData.append((float(lineChunks[0]),float(lineChunks[1])))
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
                    # Skipping empty lines
                    pass
        else:
            print "No file"
    return chromData

def openFile(fig,canvas):
    """ TODO
    """
    global inputFile
    file_path = tkFileDialog.askopenfilename()
    inputFile = file_path
    if not file_path:
        pass
    else:
        data = openChrom(file_path)
        try:
            data = data.tolist()
        except AttributeError:
            pass
        if data:
            writeData(data,file_path)
            plotData(data,fig,canvas,file_path)

def plotData(data,fig,canvas,file_path):
    """ TODO
    """
    x_array = []
    y_array = []
    for i in data:
        x_array.append(i[0])
        y_array.append(i[1])
    fig.clear()
    axes = fig.add_subplot(111)
    line, = axes.plot(x_array,y_array,label=os.path.split(file_path)[-1])
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles,labels)
    axes.get_xaxis().get_major_formatter().set_useOffset(False)
    axes.set_xlabel("Time [m]")
    axes.set_ylabel("Intensity [au]")
    canvas.draw()

def quantifyChrom(fig, canvas):
    """ TODO
    
    This is super prelimenary, should/will produce a lot more values
    """
    peakList = tkFileDialog.askopenfilename()
    peaks = []
    with open(peakList,'r') as fr:
        for line in fr:
            line = line.rstrip("\n").split("\t")
            peaks.append((str(line[0]), float(line[1]), float(line[2])))
    data = readData()
    time, intensity = zip(*data[0][1])
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
    tkMessageBox("Finished (?)","The selected chromatogram has been quantified, the results have been saved in 'output.txt'. Please note that this is a very early build that only performs a very basic quantitation.")
    with open('output.txt','w') as fw:
        for i in results:
            fw.write(str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\n")

def readData():
    """ TODO
    """
    data = []
    with open('temp/tracebuffer.txt','r') as fr:
        for line in fr:
            if line[0] == ">" and len(data) > 0:
                data.append(name, dataBuffer)
                name = line[2:].rstrip("\n")
                dataBuffer = []
            elif line[0] == ">" and len(data) == 0:
                name = line[2:].rstrip("\n")
                dataBuffer = []
            else:
                chunks = line.rstrip().split("\t")
                dataBuffer.append((float(chunks[0]),int(chunks[1])))
        data.append((name, dataBuffer))
    return data

def saveChrom():
    """ TODO
    """
    data = readData()
    saveFile = tkFileDialog.asksaveasfilename()
    with open(saveFile,'w') as fw:
        for i in data[0][1]:
            fw.write(str(i[0])+"\t"+str(i[1])+"\n")

def smoothChrom(fig, canvas):
    """ TODO
    """
    data = readData()
    time, intensity = zip(*data[0][1])
    new = savgol_filter(intensity,21,3)
    newData = zip(time,new)
    
    # Plot & Write Data to Disk  
    multiData = [(os.path.split(data[0][0])[-1], data[0][1]),(os.path.split(data[0][0])[-1]+" (Smoothed)",newData)]
    plotMultiData(fig,canvas,multiData)
    writeData(newData,os.path.split(data[0][0])[-1]+" (Smoothed)")

def writeData(data,file_path):
    """ TODO
    """
    with open('temp/tracebuffer.txt','w') as fw:
        fw.write(">>"+str(file_path)+"\n")
        for i in data:
            fw.write(str(i[0])+"\t"+str(int(i[1]))+"\n")
