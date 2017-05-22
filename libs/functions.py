#! /usr/bin/env python

# General imports
from datetime import datetime
from matplotlib.pyplot import gca
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter
from Tkinter import *
import bisect
import glob
import math
import numpy as np
import operator
import os
import re
import sys
import tkFileDialog
import tkMessageBox

# Custom libraries
sys.path.append('..')
import batchFunctions
import HappyTools

# General variables
points = 100
start = 10
end = 60
baselineOrder = 1
backgroundWindow = 1
nobanStart = 0.25
slicepoints = 5
peakDetectionMin = 0.01
createFigure = "True"

# Advanced variables
noise = "RMS"                     # Accepts: RMS, MM
backgroundNoiseMethod = "MT"      # Accepts: NOBAN, MT

###########
# Classes #
###########

################################################################################################
# Tooltip code - Taken from http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml #
################################################################################################
class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 27
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        try:
            # For Mac OS
            tw.tk.call("::tk::unsupported::MacWindowStyle",
                       "style", tw._w,
                       "help", "noActivates")
        except TclError:
            pass
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      wraplength=500, font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

#############
# Functions #
#############
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

def backgroundNoise(data):
    """ Background and noise determination based on MassyTools
        method.
    """
    background = sys.maxint
    for index,i in enumerate(data[:-slicepoints]):
        buffer = data[index:index+slicepoints]
        if np.mean(buffer) < background:
            background = np.mean(buffer)
            if noise == "MM":
                currNoise = max(buffer)-min(buffer)
            elif noise == "RMS":
                currNoise = np.std(buffer)
    return {'Background': background, 'Noise': currNoise}

def baselineCorrection(fig,canvas):
    """ TODO
    """
    data = readData()

    # Background determination
    background = []
    chunks = [data[0][1][x:x+points] for x in xrange(0, len(data[0][1]), points)]
    for i in chunks:
        buff1, buff2 = zip(*i)
        min_index, min_value = min(enumerate(buff2), key=operator.itemgetter(1))
        if buff1[0] > start and buff1[-1] < end:
            background.append((buff1[min_index], buff2[min_index]))
    time, intensity = zip(*background)
    newX = np.linspace(min(time), max(time),100)
    func = np.polyfit(time, intensity, baselineOrder)
    p = np.poly1d(func)

    # Transform
    x = []
    newChromIntensity = []
    for i in data[0][1]:
        x.append(i[0])
        newChromIntensity.append(int(i[1]-p(i[0])))
        
    # Uplift
    foo1, foo2 = zip(*data[0][1])
    low = bisect.bisect_left(foo1, start)
    high = bisect.bisect_right(foo1, end)
    offset = abs(min(min(newChromIntensity[low:high]),0))
    newData = zip(foo1,[x+offset for x in newChromIntensity])

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
        """ TODO
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
                try:
                    refPeaks.append((str(line[0]), float(line[1]), float(line[2])))
                except ValueError:
                    pass
    except IOError:
        tkMessageBox.showinfo("File Error","The selected reference file could not be opened.")

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

def createToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

def gaussFunction(x, *p):
    """ TODO
    """
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def getSettings():
    """ TODO
    """
    with open(HappyTools.settings,'r') as fr:
        for line in fr:
            chunks = line.strip('\n').split('\t')
            if chunks[0] == "points:":
                global points
                points = int(chunks[1])
            elif chunks[0] == "start:":
                global start
                start = float(chunks[1])
            elif chunks[0] == "end:":
                global end
                end = float(chunks[1])
            elif chunks[0] == "baselineOrder:":
                global baselineOrder
                baselineOrder = int(chunks[1])
            elif chunks[0] == "backgroundWindow:":
                global backgroundWindow
                backgroundWindow = float(chunks[1])
            elif chunks[0] == "noise:":
                global noise
                noise = str(chunks[1])
            elif chunks[0] == "nobanStart:":
                global nobanStart
                nobanStart = float(chunks[1])
            elif chunks[0] == "slicepoints:":
                global slicepoints
                slicepoints = int(chunks[1])
            elif chunks[0] == "createFigure:":
                global createFigure
                createFigure = str(chunks[1])

def noban(data):
    """ NOBAN implementation based on Jansen et al, 2016.
    """
    def calcValues(sortedData,currSize,currAverage,currNoise,increment):
        currSize += increment
        currAverage =  np.average(sortedData[0:currSize])
        if noise == "MM":
            currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
        elif noise == "RMS":
            currNoise = np.std(sortedData[0:currSize])
        return currSize,currAverage,currNoise

    sortedData = sorted(data)
    startSize = int(nobanStart * float(len(sortedData)))
    currSize = startSize
    currAverage = np.average(sortedData[0:currSize])
    if noise == "MM":
        currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
    elif noise == "RMS":
        currNoise = np.std(sortedData[0:currSize])
    directionFlag = 0
    # <<NOBAN-V2>>
    # This algorithm now includes faster convergence (by starting with a large step
    # and taking smaller steps the closer we get to the minimum).
    for k in range(0,len(sortedData)-(startSize+1)):
        remainder = len(sortedData) - currSize
        try:
            if sortedData[currSize+int(math.ceil(0.1*remainder))] < currAverage + 3 * currNoise:
                directionFlag == 1
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, int(math.ceil(0.1*remainder)))
            elif sortedData[currSize+int(math.ceil(0.05*remainder))] < currAverage + 3 * currNoise:
                directionFlag == 1
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, int(math.ceil(0.05*remainder)))
            elif sortedData[currSize+1] < currAverage + 3 * currNoise:
                directionFlag == 1
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, 1)
            elif sortedData[currSize-int(math.ceil(0.1*remainder))] > currAverage + 3 * currNoise and directionFlag == 0:
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, -int(math.ceil(0.1*remainder)))
            elif sortedData[currSize-int(math.ceil(0.05*remainder))] > currAverage + 3 * currNoise and directionFlag == 0:
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, -int(math.ceil(0.05*remainder)))
            elif sortedData[currSize-1] > currAverage + 3 * currNoise and directionFlag == 0:
                currSize, currAverage, currNoise = calcValues(sortedData, currSize, currAverage, currNoise, -1)
            else:
                break
        except IndexError:
            break
    return {'Background': currAverage, 'Noise': currNoise}

def openChrom(file):
    """ TODO
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
                        print "Omitting line: "+str(line)
                        
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
            print "Incorrect inputfile format, please upload a raw data 'txt' or 'arw' file."
    return chromData

def openFile(fig,canvas):
    """ TODO
    """
    #global inputFile
    file_path = tkFileDialog.askopenfilename()
    #inputFile = file_path
    if not file_path:
        pass
    else:
        data = openChrom(file_path)
        try:
            data = data.tolist()
        except AttributeError:
            pass
        if data:
            writeData(data,os.path.split(file_path)[-1])
            plotData(data,fig,canvas,file_path)

def peakDetection(fig,canvas):
    """ TODO
    """
    data = readData()
    
    x_data,y_data = zip(*data[0][1])
    orig_x = x_data
    orig_y = y_data
    low = bisect.bisect_left(x_data,start)
    high = bisect.bisect_right(x_data,end)
    x_data = x_data[low:high]
    y_data = y_data[low:high] 
    if backgroundNoiseMethod == "NOBAN":
        NOBAN = noban(y_data)
    elif backgroundNoiseMethod == "MT":
        NOBAN = backgroundNoise(y_data)

    functions = []
    cutoff = peakDetectionMin*(max(y_data)-NOBAN['Background'])

    newX = np.linspace(x_data[0], x_data[-1], 25000*(x_data[-1]-x_data[0]))
    f = InterpolatedUnivariateSpline(x_data, y_data)
    fPrime = f.derivative()
    newY = f(newX)
    newPrimeY = fPrime(newX)
    maxm = argrelextrema(newPrimeY, np.greater)
    minm = argrelextrema(newPrimeY, np.less)
    breaks = maxm[0].tolist() + minm[0].tolist()
    breaks = sorted(breaks)

    counter = 0
    while max(y_data)-NOBAN['Background'] > cutoff:
        counter += 1
        print "Fitting peak: "+str(counter)
        f = InterpolatedUnivariateSpline(x_data, y_data)
        fPrime = f.derivative()
        fPP = fPrime.derivative()
        newY = f(newX)
        maxPoint = 0
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
            except IndexError:
                if max(newY[breaks[index]:-1]) > maxPoint:
                    maxPoint = max(newY[breaks[index]:-1])
                    xData = newX[breaks[index]:-1]
                    yData = [x - NOBAN['Background'] for x in newY[breaks[index]:-1]]
                pass
        # Gaussian fit on main points
        newGaussX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        p0 = [np.max(yData), xData[np.argmax(yData)],0.1]
        try:
            coeff, var_matrix = curve_fit(gaussFunction, xData, yData, p0)
            newGaussY = gaussFunction(newGaussX, *coeff)
            newGaussY = [x + backGround for x in newGaussY]
        except:
            pass
        x_buff = []
        y_buff = []
        clac = 0.01 * max(newGaussY)
        for bla,ka in enumerate(newGaussY):
            #if ka > NOBAN['Background']+NOBAN['Noise']:
            if ka > clac:
                x_buff.append(newGaussX[bla])
                y_buff.append(ka)
        newGaussX = x_buff
        newGaussY = y_buff
        functions.append(("Peak: "+str(xData[np.argmax(yData)]),zip(newGaussX,newGaussY)))
        # Adjust data
        new_y = []
        for index,j in enumerate(x_data):
            new_y.append(y_data[index] - gaussFunction(j,*coeff))
        if max(new_y)-NOBAN['Background'] == max(y_data)-NOBAN['Background']:
            break
        y_data = new_y
            
    # Plotting
    fig.clear()
    axes = fig.add_subplot(111)
    axes.plot(orig_x,orig_y, 'b',alpha=0.5)
    for i in functions:
        xd = []
        yd = []
        for j in i[1]:
            xd.append(j[0])
            yd.append(j[1])
        axes.plot(xd,yd,label=os.path.split(i[0])[-1])
        axes.fill_between(xd, 0, yd,alpha=0.2)
        #ax.fill_between(x, y1, y2, where=y2 <= y1, facecolor='red'
    axes.set_xlabel("Time [m]")
    axes.set_ylabel("Intensity [au]")
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles,labels)
    canvas.draw()
    """
    # Gaussian fit on all points
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plt.plot(orig_x, orig_y, 'b*')
    #for i in functions:
    #    x_data, y_data = zip(*i)
    #    plt.plot(x_data,y_data)
    #plt.plot((newX[0],newX[-1]),(backGround,backGround),'red',linestyle='dashed',alpha=0.5)
    #plt.show()"""
        

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
            try:
                peaks.append((str(line[0]), float(line[1]), float(line[2])))
            except ValueError:
                pass
    data = {'Name':readData()[0][0],'Data':readData()[0][1]}
    time, intensity = zip(*data['Data'])
    results = []
    for i in peaks:
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])
        peakArea = 0
        residual = 0
        signalNoise = "Nan"
        # Get signal-to-noise
        lowBackground = bisect.bisect_left(time,max(i[1]-backgroundWindow,start))
        highBackground = bisect.bisect_right(time,min(i[1]+backgroundWindow,end))
        backgroundData = intensity[lowBackground:low]+intensity[high:highBackground]
        if backgroundNoiseMethod == "NOBAN":
            NOBAN = noban(backgroundData)
        elif backgroundNoiseMethod == "MT":
            NOBAN = backgroundNoise(backgroundData)
        signalNoise = (max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise']
        # Get peak Area
        for index,j in enumerate(intensity[low:high]):
            try:
                #peakArea += j * (time[low+index]-time[low+index-1])
                peakArea += max(j-NOBAN['Background'],0) * (time[low+index]-time[low+index-1])
            except IndexError:
                continue
        # Gaussian fit (to get residuals)
        x_data = np.array(time[low:high])
        y_data = np.array(intensity[low:high])
        newX = np.linspace(x_data[0], x_data[-1], 2500*(x_data[-1]-x_data[0]))
        p0 = [np.max(y_data), x_data[np.argmax(y_data)],0.1]
        try:
            coeff, var_matrix = curve_fit(gaussFunction, x_data, y_data, p0)
            newY = gaussFunction(newX, *coeff)
            # Get residuals
            for index,j in enumerate(time[low:high]):
                residual += abs(intensity[index]-gaussFunction(j, *coeff))**2
            residual = math.sqrt(residual)
        except RuntimeError:
            residual = "Nan"
    data['Name'] = str(data['Name'].split('.')[0])+".raw"
    with open(data['Name'],'w') as fw:
        fw.write("Name\tTime\tPeak Area\tS/N\tGaussian Residual RMS\n")
        for i in results:
            fw.write(str(i['Peak'])+"\t"+str(i['Time'])+"\t"+str(i['Area'])+"\t"+str(i['S/N'])+"\t"+str(i['Residual'])+"\n")
    tkMessageBox.showinfo("Status Message", "Quantitation finished on "+str(datetime.now()))

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

def settingsPopup():
    """ TODO
    """

    figureVariable = StringVar()
    figureVariable.set(createFigure)

    def close():
        """ TODO
        """
        global points
        global start
        global end
        global baselineOrder
        global backgroundWindow
        global backgroundNoise
        global nobanStart
        global slicepoints
        global createFigure
        points = int(pointsWindow.get())
        start = float(startWindow.get())
        end = float(endWindow.get())
        baselineOrder = int(baselineOrderWindow.get())
        backgroundWindow = float(backgroundWindowWindow.get())
        nobanStart = float(nobanWindow.get())
        slicepoints = int(slicepointsWindow.get())
        createFigure = str(figureVariable.get())
        top.destroy()
    
    def save():
        """ TODO
        """
        with open(HappyTools.settings,'w') as fw:
            fw.write("points:\t"+str(int(pointsWindow.get()))+"\n")
            fw.write("start:\t"+str(float(startWindow.get()))+"\n")
            fw.write("end:\t"+str(float(endWindow.get()))+"\n")
            fw.write("baselineOrder:\t"+str(int(baselineOrderWindow.get()))+"\n")
            fw.write("backgroundWindow:\t"+str(float(backgroundWindowWindow.get()))+"\n")
            fw.write("nobanStart:\t"+str(float(nobanWindow.get()))+"\n")
            fw.write("noise:\t"+str(noise)+"\n")
            fw.write("slicepoints:\t"+str(slicepoints)+"\n")
            fw.write("createFigure:\t"+str(figureVariable.get())+"\n")
        
    top = Tk.top = Toplevel()
    top.protocol( "WM_DELETE_WINDOW", lambda: close())
    
    pointsLabel = Label(top, text="Datapoints", font="bold")
    pointsLabel.grid(row=0, column=0, sticky=W)
    pointsWindow = Entry(top)
    pointsWindow.insert(0, points)
    pointsWindow.grid(row=0, column=1, sticky=W)
    
    startLabel = Label(top, text="Start", font="bold")
    startLabel.grid(row=1, column=0, sticky=W)
    startWindow = Entry(top)
    startWindow.insert(0, start)
    startWindow.grid(row=1, column=1, sticky=W)
    
    endLabel = Label(top, text="End", font="bold")
    endLabel.grid(row=2, column=0, sticky=W)
    endWindow = Entry(top)
    endWindow.insert(0, end)
    endWindow.grid(row=2, column=1, sticky=W)

    baselineOrderLabel = Label(top, text="Baseline Order", font="bold")
    baselineOrderLabel.grid(row=3, column=0, sticky=W)
    baselineOrderWindow = Entry(top)
    baselineOrderWindow.insert(0, baselineOrder)
    baselineOrderWindow.grid(row=3, column=1, sticky=W)
    
    backgroundWindowLabel = Label(top, text="Background Window", font="bold")
    backgroundWindowLabel.grid(row=4, column=0, sticky=W)
    backgroundWindowWindow = Entry(top)
    backgroundWindowWindow.insert(0, backgroundWindow)
    backgroundWindowWindow.grid(row=4, column=1, sticky=W)

    nobanLabel = Label(top, text="NOBAN Start", font="bold")
    nobanLabel.grid(row=5, column=0, sticky=W)
    nobanWindow = Entry(top)
    nobanWindow.insert(0, nobanStart)
    nobanWindow.grid(row=5, column=1, sticky=W)

    slicepointsLabel = Label(top, text="MT Slice points", font="bold")
    slicepointsLabel.grid(row=6, column=0, sticky=W)
    slicepointsWindow = Entry(top)
    slicepointsWindow.insert(0, slicepoints)
    slicepointsWindow.grid(row=6, column=1, sticky=W)

    figureLabel = Label(top, text="Create figure for each analyte", font="bold")
    figureLabel.grid(row=7, column=0, sticky=W)
    options = ["True", "False"]
    figureWindow = OptionMenu(top, figureVariable, *options)
    figureWindow.grid(row=7, column=1, sticky=W)

    saveButton = Button(top, text="Save", command=lambda: save())
    saveButton.grid(row=8, column=0, sticky=W)
    closeButton = Button(top, text="Close", command=lambda: close())
    closeButton.grid(row=8, column=1, sticky=E)

    # Tooltips
    createToolTip(pointsLabel,"The number of data points that is used to determine the baseline. Specifically, this setting specifies how large each segment of the whole chromatogram will be to identify the lowest data point per window, i.e. a setting of 100 means that the chromatogram is split into segments of 100 data points per segment.")
    createToolTip(startLabel,"This setting tells the program from which time point it is supposed to begin processing, this setting should be set in such a way that it is before the analytes of interest but after any potential big increases or decrease in intensity.")
    createToolTip(endLabel,"This setting tells the program until which time point it is supposed to begin processing, this setting should be set in such a way that it is after the analytes of interest but before any potential big increases or decrease in intensity.")
    createToolTip(baselineOrderLabel,"This setting tells the program what sort of function should be used to correct the baseline. A value of 1 refers to a linear function, while a value of 2 refers to a quadratic function. We advise to use a linear function as the need for any higher order function indicates an unexpected event in the chromatography.")
    createToolTip(backgroundWindowLabel,"This setting tells the program the size of the region that will be examined to determine the background. A value of 1 means that the program will look from 20.0 to 21.0 minutes and 21.4 to 22.4 for an analyte that elutes from 21.0 to 21.4 minutes.")
    createToolTip(nobanLabel,"This setting specifies the initial estimate for the NOBAN algorithm, specifically a value of 0.25 means that the lowest 25% of all data points will be used as an initial estimate for the background. This value should be changed depending on how many signals there are in the chromatogram, e.g. in a crowded chromatogram this value should be low.")
    createToolTip(slicepointsLabel,"The number of consecutive data points that will be used to determine the background and noise using the MT method. The MT method will scan all datapoints that fall within the background window (specified above) to find the here specified number of consecutive data points that yield the lowest average intensity, the average of these data points is then used as background while the standard deviation of these data points is used as the noise.")

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

def updateProgressBar(bar, variable, index, length):
    """ TODO
    """
    variable.set(str(int((float(index)/float(length))*100))+"%")
    bar["value"] = int((float(index)/float(length))*100)
    bar.update()

def writeData(data,file_path):
    """ TODO
    """
    with open('temp/tracebuffer.txt','w') as fw:     
        fw.write(">>"+str(file_path)+"\n")
        for i in data:
            fw.write(str(i[0])+"\t"+str(int(i[1]))+"\n")
