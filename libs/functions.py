#! /usr/bin/env python

###################
# General imports #
###################
from datetime import datetime
from matplotlib.pyplot import gca
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import PchipInterpolator, Akima1DInterpolator
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
import shutil
import sys
import tkFileDialog
import tkMessageBox

####################
# Custom libraries #
####################
sys.path.append('..')
import batchFunctions
import HappyTools

##################### 
# General variables #
#####################
points = 100
start = 10
end = 60
baselineOrder = 1
backgroundWindow = 1
nobanStart = 0.25
slicepoints = 5
peakDetectionMin = 0.05
peakDetectionEdge = "Sigma"
peakDetectionEdgeValue = 2.0
createFigure = "True"
minPeaks = 4
minPeakSN = 27
useUPC = "False"

######################
# Advanced variables #
######################
# decimalNumbers specifies the number of decimal numbers to be used
# in the files created by HappyTools
decimalNumbers = 6
# min_improvement specifies the minimum improvement in RMS that a 
# calibration function must yield, to accept the more complex 
# calibration function as an acceptable function.
min_improvement = 0.05
# use_interpolation specifies if calibration is allowed to use
# interpolation methods as well, which will always return an RMS of 0
# but tend to be less smooth than regular functions. However, the
# allowed interpolation functions are guaranteed to be monotone.
use_interpolation = False
# noise specifies the method that is used to determine the noise,
# acceptable values are 'RMS' for root-mean-square which is based on
# the standard deviation of a signal (much like radio noise) or 'MM'
# for the difference between the maximum and mimum, much like classical
# chromatography
noise = "RMS"
# backgroundNoiseMethod specifies the how the region is determined
# from which both the background and noise values are extracted. The
# options included 'MT' for the MassyTools method (5 sequential lowest
# datapoints near the peak of interest) or 'NOBAN' for the normal
# distrubtion based background and noise determined first described
# in the LaCyTools paper (assumption that background and noise are
# normally distributed and any datapoint that is not in the normal
# distribution is a signal
backgroundNoiseMethod = "MT"

####################
# Output variables #
####################
root = Tk()
root.withdraw()
outputWindow = IntVar()
absInt = IntVar()
relInt = IntVar()
bckSub = IntVar()
bckNoise = IntVar()
peakQual = IntVar()

###########
# Classes #
###########
class powerLawCall:
    def __init__(self,x,a,b,c):
        self.a = a
        self.b = b
        self.c = c
        self.f = lambda x: a*x**b+c
    def __call__(self,x):
        return self.f(x)
    def describe(self):
        return str(self.a)+"*X^"+str(self.b)+"+"+str(self.c)

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
    """Ask for a file and draw it on the existing canvas."""
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
    """Return the background and noise.

    This function determines the average and the standard deviation or
    the maximum difference of all segments of data, where each segment
    has the length specified in the slicepoints parameter.

    Keyword arguments:
    data -- list of intensities
    """
    background = sys.maxint
    currNoise = 0
    for index,i in enumerate(data[:-slicepoints]):
        buffer = data[index:index+slicepoints]
        if np.mean(buffer) < background:
            background = np.mean(buffer)
            if noise == "MM":
                currNoise = max(buffer)-min(buffer)
            elif noise == "RMS":
                currNoise = np.std(buffer)
    if currNoise == 0:
        currNoise = 1
    return {'Background': background, 'Noise': currNoise}

def baselineCorrection(fig,canvas):
    """Perform baseline correction and draw the corrected chromatogram.

    This function determines the baseline of a chromatogram between
    two timepoints, specified with the start and end parameter. The
    chromatogram is split into segments of a length specified in the
    points parameter. The lowest intensity of each segment is used to
    determine a function of the order specified in the baselineOrder
    using the numpy.polyfit function. The original chromatogram is then
    transformed by subtracting the function from the original data. The
    resulting chromatogram might have negative intensities between the
    start and end timepoints, the minimum intensity within that region
    is used to uplift the entire chromatogram. The transformed and
    uplifted chromatogram is written to disk and plotted together with
    the original chromatogram on the canvas.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
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
    time = [a for a,b in data[0][1]]
    newChromIntensity = [b-p(a) for a,b in data[0][1]]
        
    # Uplift
    low = bisect.bisect_left(time, start)
    high = bisect.bisect_right(time, end)
    offset = abs(min(min(newChromIntensity[low:high]),0))
    newData = zip(time,[x+offset for x in newChromIntensity])

    # Plot & Write Data to Disk  
    multiData = [(os.path.split(data[0][0])[-1], data[0][1]),(os.path.split(data[0][0])[-1]+" (BC)",newData)]
    plotMultiData(fig,canvas,multiData)
    writeData(newData,os.path.split(data[0][0])[-1]+" (BC)")

def plotMultiData(fig,canvas,data):
    """Plot all chromatograms in data on the canvas.

    This function first clears the canvas and then draws all the
    chromatograms that are in the data list on the canvas.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
    data -- list of tuples, consisting of two numbers per tuple
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

def powerLaw(x,a,b,c):
    """ TODO
    """
    penalty = 0
    if b > 2.:
        penalty = abs(b-1.)*10000
    if b < 0.:
        penalty = abs(2.-b)*10000
    return a*x**b + c + penalty

def batchPlot(fig,canvas):
    """Read and plot all chromatograms in a directory.

    This function asks the user to select a directory from which the
    function will read all the files that are specified in the 
    CALIBRATION_FILETYPES paramater of the batchFunctions file and plot
    them to the canvas.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
    """
    folder_path = tkFileDialog.askdirectory()
    if folder_path:
        filesGrabbed = []
        for files in batchFunctions.CALIBRATION_FILETYPES:
            for file in glob.glob(str(os.path.join(folder_path,files))):
                if os.path.basename(file) not in batchFunctions.EXCLUSION_FILES:
                    if openChrom(file):
                        filesGrabbed.append(file)

        data = []
        for file in filesGrabbed:
            data.append((str(file),openChrom(file)))

        if data:
            fig.clear()
            axes = fig.add_subplot(111)
            for i in data:
                x_array, y_array = zip(*i[1])
                axes.plot(x_array,y_array,label=str(os.path.split(i[0])[-1]))
            axes.legend()
        canvas.draw()

def batchPlotNorm(fig,canvas):
    """Read and plot all chromatograms in a directory.

    This function asks the user to select a directory from which the
    function will read all the files that are specified in the 
    CALIBRATION_FILETYPES paramater of the batchFunctions file. The
    function will then find the lowest and maximum intensities between
    the start and end variable, normalize all chromatograms and plot
    them to the canvas.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
    """
    folder_path = tkFileDialog.askdirectory()
    if folder_path:
        filesGrabbed = []
        for files in batchFunctions.CALIBRATION_FILETYPES:
            for file in glob.glob(str(os.path.join(folder_path,files))):
                if os.path.basename(file) not in batchFunctions.EXCLUSION_FILES:
                    if openChrom(file):
                        filesGrabbed.append(file)

        data = []
        for file in filesGrabbed:
            chromData = openChrom(file)

            # Background determination
            background = []
            chunks = [chromData[x:x+points] for x in xrange(0, len(chromData), points)]
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
            time = [a for a,b in chromData]
            newChromIntensity = [b-p(a) for a,b in chromData]

            # Uplift
            low = bisect.bisect_left(time, start)
            high = bisect.bisect_right(time, end)
            offset = abs(min(min(newChromIntensity[low:high]),0))
            newIntensity = [x+offset for x in newChromIntensity]

            # Normalize
            correction = max(newIntensity[low:high])
            normIntensity = [x/correction for x in newIntensity]
            newData = zip(time,normIntensity)
            data.append((str(file),newData))

        # Plot
        if data:
            fig.clear()
            axes = fig.add_subplot(111)
            for i in data:
                x_array, y_array = zip(*i[1])
                axes.plot(x_array,y_array,label=str(os.path.split(i[0])[-1]))
            axes.legend()
        canvas.draw()

def batchPopup():
    """Create a batch processing pop-up.

    This function creates a new tkinter window that is used to control
    the batch processing. Specifically, it allows the user to select a
    calibration file, an analyte file, select the desired outputs (by
    calling the outputPopup function) and starting the batch process.

    Keyword arguments:
    none
    """

    calFile = StringVar()
    analFile = StringVar()
    batchFolder = StringVar()

    def close():
        """Close the batch processing pop-up.
        """
        top.destroy()

    def setCalibrationFile():
        """Ask for the calibration file.
        """
        calFile.set(tkFileDialog.askopenfilename(title="Calibration File"))

    def setAnalyteFile():
        """Ask for the analyte file.
        """
        analFile.set(tkFileDialog.askopenfilename(title="Analyte File"))

    def setBatchFolder():
        """Ask for the batch folder.
        """
        batchFolder.set(tkFileDialog.askdirectory(title="Batch Folder"))

    def run():
        """Start the batch process.
        """
        batchFunctions.batchProcess(calFile, analFile, batchFolder)

    top = Tk.top = Toplevel()
    top.title("HappyTools "+str(HappyTools.version)+" Batch Process")
    top.protocol("WM_DELETE_WINDOW", lambda: close())

    calibrationButton = Button(top, text="Calibration File", width=20, command=lambda: setCalibrationFile())
    calibrationButton.grid(row=1, column=0, sticky=W)
    calibrationLabel = Label(top, textvariable=calFile, width=20)
    calibrationLabel.grid(row=1, column=1)

    analyteButton = Button(top, text="Analyte File", width=20, command=lambda: setAnalyteFile())
    analyteButton.grid(row=2, column=0, sticky=W)
    analyteLabel = Label(top, textvariable=analFile, width=20)
    analyteLabel.grid(row=2, column=1)

    batchButton = Button(top, text="Batch Directory", width=20, command=lambda: setBatchFolder())
    batchButton.grid(row=3, column=0, sticky=W)
    batchLabel = Label(top, textvariable=batchFolder, width=20)
    batchLabel.grid(row=3, column=1, sticky=W)
    
    outputButton = Button(top, text="Output Options", command=lambda: outputPopup())
    outputButton.grid(row=4, column=0, columnspan=2, sticky=E+W)

    runButton = Button(top, text="Run", width=20, command=lambda: run())
    runButton.grid(row=5, column=0, sticky=W)
    closeButton = Button(top, text="Close", width=20, command=lambda: close())
    closeButton.grid(row=5, column=1, sticky=E)

    # Tooltips
    createToolTip(calibrationButton,"This button will allow you to select your calibration file, the program expects a "+
            "tab separated text file where each line consists of a peak ID, peak RT and a RT window.")
    createToolTip(analyteButton,"This button will allow you to select your analyte file, the program expects a tab separated "+
            "text file where each line consists of a peak ID, peak RT and a RT window.")
    createToolTip(outputButton,"This button will open another window in which you can select which outputs you want "+
            "HappyTools to show in the final summary.")

def chromCalibration(fig,canvas):
    """Ask for a reference file and calibrate the current chromatogram.

    This function will first ask the user to select a reference file 
    using a tkinter filedialog.  The function will then find the highest
    intensity timepoint for each calibrant window. The actual
    calibration is achieved by fitting a second degree polynomial
    through the observed and expected retention times and applying the
    formula on the original chromatogram, with the new retention time
    being cast into a float of a user defined number of decimal numbers.
    The function finishes by plotting the calibrated chromatograom on
    top of the original chromatogram and writing the calibrated
    chromatogram to the temporary file.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
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
    if useUPC == False:
        z = np.polyfit(observedTime,expectedTime,2)
        f = np.poly1d(z)
    elif useUPC == True:
        f = ultraPerformanceCalibration(observedTime,expectedTime, time[0], time[-1])
    calibratedData = zip(f(time),intensity)

    # Plot & Write Data to Disk  
    multiData = [(os.path.split(data[0][0])[-1], data[0][1]),(os.path.split(data[0][0])[-1]+" (Cal)",calibratedData)]
    plotMultiData(fig,canvas,multiData)
    writeData(calibratedData,os.path.split(data[0][0])[-1]+" (Cal)")

def fileCleanup():
    """Clean up the temporary files.

    This function will search open the temporary directory and attempt
    to delete all files in that directory. A warning will be thrown
    when a temporary file was open in another program.

    Keyword arguments:
    none
    """
    search = os.path.join(str(os.getcwd()),"temp","*.*")
    files = glob.glob(search)
    for f in files:
        try:
            os.remove(f)
        except OSError:
            tkMessageBox.showinfo("File Error", "A temporary file is open in another program. Please exit HappyTools "+
                    "and close all temporary files before running HappyTools.")

def fwhm(coeff):
    """Calculate the FWHM.
    
    This function will calculate the FWHM based on the following formula
    FWHM = 2*sigma*sqrt(2*ln(2)). The function will return a dictionary
    with the fwhm ('fwhm'), the Gaussian peak center ('center') and the
    +/- width, from the peak center ('width').
    
    Keyword arguments:
    coeff -- coefficients as calculated by SciPy curve_fit
    """
    fwhm = abs(2*coeff[2]*math.sqrt(2*math.log(2)))
    width = 0.5*fwhm
    center = coeff[1]
    return {'fwhm':fwhm, 'width':width, 'center':center}

def checkAccess(all_directories):
    """Check disk access
    
    This function will check if HappyTools has read/write access to
    all the folders belonging to HappyTools.
    
    Keyword arguments:
    None
    """
    access = True
    for directory in all_directories:
        if not os.access(directory, os.W_OK):
            access = False
    return access

def createToolTip(widget, text):
    """Create a tooltip.

    This function will create a tooltip and assign it to the widget that
    was handed to this function. The widget will then show the provided
    text upon a mouseover.

    Keyword arguments:
    widget -- tkinter object
    text -- string
    """ 
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

def determineCalibrants(functions):
    """ Automatically determine suitable calibrant peaks.
    
    This function is part of the automated peak detection
    functionality, where it attempts to determine the most suitable
    n number of calibrations (n is the user defined setting of minimum
    number of peaks for calibration). The function splits the time
    between the first and last detected peaks into n chunks, determines
    the highest peak in each chunk and classifies these local maxima
    as potential calibrants.
    
    Keyword Arguments:
    functions -- A list of dicts, containing Peak, Data and FWHM
    """
    calibrants = []
    timeChunks = []
    timeRange = functions[-1]['Data'][-1][0] - functions[0]['Data'][0][0]
    timeChunk = timeRange / minPeaks
    for i in range(0, minPeaks):
        timeChunks.append((functions[0]['Data'][0][0]+i*float(timeChunk),functions[0]['Data'][0][0]+(i+1)*float(timeChunk)))
    for i in timeChunks:
        maxIntensity = 0
        calBuffer = None
        for j in functions:
            try:
                X, Y = zip(*j['Data'])
            except:
                print(j)
            peakCenter = X[Y.index(max(Y))]
            peakIntensity = max(Y)
            if peakIntensity > maxIntensity and peakCenter > i[0] and peakCenter < i[1]:
                maxIntensity = peakIntensity
                calBuffer = j
        if calBuffer:
            calibrants.append(calBuffer)  
    return calibrants

def gaussFunction(x, *p):
    """Define and return a Gaussian function.

    This function returns the value of a Gaussian function, using the
    A, mu and sigma value that is provided as *p.

    Keyword arguments:
    x -- number
    p -- A, mu and sigma numbers
    """
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def getPeakList(fileName):
    """Read and parse the peak file and return a list of peaks.
    
    This function opens the file that is specified in 'fileName', and 
    reads it on a line per line basis. The function will split each 
    line on '\t' prior to trying to append the parts to the 'peaks'
    list. The function will write a message to the logFile if logging
    is enabled if the previous mentioned try goes to the except clause.
    
    Keyword argments:
    fileName -- string
    """
    peaks = []
    try:
        with open(fileName,'r') as fr:
            for line in fr:
                line = line.rstrip("\n").split("\t")
                try:
                    peaks.append((str(line[0]), float(line[1]), float(line[2])))
                except ValueError:
                    if HappyTools.logging == True and HappyTools.logLevel > 1:
                        with open(HappyTools.logFile,'a') as fw:
                            fw.write(str(datetime.now().replace(microsecond=0))+"\tIgnoring line: "+str(line)+" from file: "+str(fileName)+"\n")
                    pass
    except IOError:
        tkMessageBox.showinfo("File Error","The selected reference file could not be opened.")
    return peaks

def getSettings():
    """Read the settings file.

    This function opens the settings file (default is HappyTools.ini),
    parses the lines of the settings file and takes the value from the
    settings file as a value for the changeable settings (e.g. the start
    variable can be read from the settings file).

    Keyword arguments:
    none
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
            #elif chunks[0] == "nobanStart:":
                #global nobanStart
                #nobanStart = float(chunks[1])
            elif chunks[0] == "slicepoints:":
                global slicepoints
                slicepoints = int(chunks[1])
            elif chunks[0] == "createFigure:":
                global createFigure
                createFigure = str(chunks[1])
            elif chunks[0] == "minPeaks:":
                global minPeaks
                minPeaks = int(chunks[1])
            elif chunks[0] == "minPeakSN:":
                global minPeakSN
                minPeakSN = int(chunks[1])
            elif chunks[0] == "peakDetectionMin:":
                global peakDetectionMin
                peakDetectionMin = float(chunks[1])
            elif chunks[0] == "peakDetectionEdge:":
                global peakDetectionEdge
                peakDetectionEdge = str(chunks[1])
            elif chunks[0] == "peakDetectionEdgeValue:":
                global peakDetectionEdgeValue
                peakDetectionEdgeValue = float(chunks[1])

def makeFunc(module):
    """This function returns a lambda function for the plugins."""
    return lambda: module.start()

def noban(data):
    """Determine background and noise using the NOBAN algorithm.

    This function is based on the NOBAN algorith, published by Jansen et
    al, in 2016. The function sorts the data by increasing intensity,
    takes an initiale estimate (defined in nobanStart) and calculates
    the background (average) and noise (root-mean-square or the maximum
    difference) from the initial estimate. The algorithm will then loop
    over the subsequent datapoints until the next datapoint falls
    outside of the current average plus three times the standard
    definition (as any point that is > 3SD is considered a signal).
    Alternatively, the function can also shrink the initial region if it
    appears that the initial estimate was too greedy. A major difference
    between the original implementation and this implementation is that
    this function should converge faster by allowing to take different
    step sizes.

    Keyword arguments:
    data -- list of numbers
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
    if currNoise == 0:
        currNoise = 1
    return {'Background': currAverage, 'Noise': currNoise}

def openChrom(file):
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

def openFile(fig,canvas):
    """Open a file and show it on the canvas.

    This function first asks the user to select a file via a file 
    dialog. The function will then call two other functions to write a
    temporary file to the disk (writeData) and to plot the selected file
    on the canvas (plotData).

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
    """
    file_path = tkFileDialog.askopenfilename()
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

def outputPopup():
    """Create a pop-up enabling output selection.

    This function creates a pop up box that allows the user to specify 
    what output should be shown in the final summary. The default value 
    for all variables is off (0) and by ticking a box it is set to on
    (1).

    Keyword arguments:
    none
    """
    if outputWindow.get() == 1:
        return
    outputWindow.set(1)

    def select_all():
        """Set all variables to on (1).
        """
        absInt.set(1)
        relInt.set(1)
        bckSub.set(1)
        bckNoise.set(1)
        peakQual.set(1)

    def select_none():
        """Set all variables to off (0).
        """
        absInt.set(0)
        relInt.set(0)
        bckSub.set(0)
        bckNoise.set(0)
        peakQual.set(0)

    def close():
        """Close the output pop-up.
        """
        outputWindow.set(0)
        top.destroy()

    top = Toplevel()
    top.protocol("WM_DELETE_WINDOW", lambda: close())
    top.title("HappyTools "+str(HappyTools.version)+" Output Options")
    selAll = Button(top, text="Select All", command=lambda: select_all())
    selAll.grid(row=0, column=0, sticky=W)
    none = Button(top, text="Select None", command=lambda: select_none())
    none.grid(row=0, column=1, sticky=E)
    text1 = Label(top, text="Base Outputs", font="bold")
    text1.grid(row=1, column=0, sticky=W)
    text2 = Label(top, text="Output Modifiers", font="bold")
    text2.grid(row=1, column=1, sticky=W)
    ai = Checkbutton(top, text=u"Analyte Intensity\u00B9", variable=absInt, onvalue=1, offvalue=0)
    ai.grid(row=2, column=0, sticky=W)
    ri = Checkbutton(top, text=u"Relative Intensity\u00B9", variable=relInt, onvalue=1, offvalue=0)
    ri.grid(row=3, column=0, sticky=W)
    pq = Checkbutton(top, text="Peak Quality Criteria", variable=peakQual, onvalue=1, offvalue=0)
    pq.grid(row=4, column=0, sticky=W)
    bn = Checkbutton(top, text="Background and Noise", variable=bckNoise, onvalue=1, offvalue=0)
    bn.grid(row=5, column=0, sticky=W)
    bck = Checkbutton(top, text=u"\u00B9Background subtracted Intensities", variable=bckSub, onvalue=1, offvalue=0)
    bck.grid(row=2, column=1, sticky=W)
    button = Button(top,text='Ok',command = lambda: close())
    button.grid(row = 6, column = 0, columnspan = 2)
    top.lift()
    return

def overlayQuantitationWindows(fig,canvas):
    """ TODO
    """
    # Prompt for peaklist
    peakList = tkFileDialog.askopenfilename()
    peaks = []
    with open(peakList,'r') as fr:
        for line in fr:
            line = line.rstrip("\n").split("\t")
            try:
                peaks.append((str(line[0]), float(line[1]), float(line[2])))
            except ValueError:
                pass

    # Read data currently on canvas
    data = {'Name':readData()[0][0],'Data':readData()[0][1]}
    time, intensity = zip(*data['Data'])

    # Plot the original data
    fig.clear()
    axes = fig.add_subplot(111)
    line, = axes.plot(time,intensity,label=data['Name'])
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles,labels)
    axes.get_xaxis().get_major_formatter().set_useOffset(False)
    axes.set_xlabel("Time [m]")
    axes.set_ylabel("Intensity [au]")
    canvas.draw()

    # Plot the quantitation windows
    for i in peaks:
        low = bisect.bisect_left(time,i[1]-i[2])
        high = bisect.bisect_right(time,i[1]+i[2])
        newTime = np.linspace(time[low], time[high],len(time[low:high]))
        f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])
        newIntensity = f(newTime)
        axes.fill_between(time[low:high], 0, newIntensity, alpha=0.5)
        axes.text(i[1], max(intensity[low:high]), i[0])
        canvas.draw()

def peakDetection(fig,canvas):
    """Detect all peaks in the currently active chromatogram.

    This function performs peak detection by fitting a Gaussian function
    through the highest data points in a chromatogram. The fitted
    function is then subtracted from the original data to yield a 
    chromatogram without the removed analyte, after which this process
    is repeated until the highest datapoint falls below the specified
    cut-off (determined by comparing the intensity of the most intense
    analyte in the original data with the intensity of the most intense
    analyte in the current (residual) data).
    
    The peak detection is based on the assumption that the first
    derivative of the data is 0 at a local maxima or minima. 
    
    <<TODO>>
    the current implementation is overly complex and can be optimized 
    and the code has to be cleaned up.

    Keyword arguments:
    fig -- matplotlib figure object
    canvas -- tkinter canvas object
    """
    data = readData()

    # Retrieve subset of data and determine the background
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

    # Determine the local maxima and minima, using first order derivative
    newX = np.linspace(x_data[0], x_data[-1], 25000*(x_data[-1]-x_data[0]))
    f = InterpolatedUnivariateSpline(x_data, y_data)
    fPrime = f.derivative()
    newY = f(newX)
    newPrimeY = fPrime(newX)
    maxm = argrelextrema(newPrimeY, np.greater)
    minm = argrelextrema(newPrimeY, np.less)
    breaks = maxm[0].tolist() + minm[0].tolist()
    breaks = sorted(breaks)

    # Determine the maximum full peak within the specified window, for the cut-off
    maxIntensity = 0
    for i in range(0, len(breaks)-2):
        maxIntensity = max(max(newY[breaks[i]:breaks[i+1]]),maxIntensity)
    cutoff = peakDetectionMin*(maxIntensity - max(NOBAN['Background'],0))

    # Detect peaks
    functions = []
    counter = 0
    while max(y_data)-NOBAN['Background'] > cutoff:
        counter += 1
        print("Fitting peak: "+str(counter))
        f = InterpolatedUnivariateSpline(x_data, y_data)
        fPrime = f.derivative()
        newY = f(newX)
        newPrimeY = fPrime(newX)
        maxm = argrelextrema(newPrimeY, np.greater)
        minm = argrelextrema(newPrimeY, np.less)
        breaks = maxm[0].tolist() + minm[0].tolist()
        breaks = sorted(breaks)
        maxPoint = 0

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
            coeff, var_matrix = curve_fit(gaussFunction, xData, yData, p0)
            newGaussY = gaussFunction(newGaussX, *coeff)
        except:
            pass

        # Limit the peak to either FWHM or a user specified Sigma value
        FWHM = fwhm(coeff)
        if peakDetectionEdge == "FWHM":
            low = bisect.bisect_left(newGaussX,coeff[1]-FWHM['width'])
            high = bisect.bisect_right(newGaussX,coeff[1]+FWHM['width'])
            try:
                newGaussX = newGaussX[low:high]
                newGaussY = newGaussY[low:high]
            except:
                pass
        elif peakDetectionEdge == "Sigma":
            low = bisect.bisect_left(newGaussX,coeff[1]-peakDetectionEdgeValue*abs(coeff[2]))
            high = bisect.bisect_right(newGaussX,coeff[1]+peakDetectionEdgeValue*abs(coeff[2]))
            try:
                newGaussX = newGaussX[low:high]
                newGaussY = newGaussY[low:high]
            except:
                pass

        # Ignore breaks (f'(x) == 0) that did not match any data (reword this)
        if newGaussX.any():
            functions.append({'Peak':newGaussX[np.argmax(newGaussY)],'Data':zip(newGaussX,newGaussY),'FWHM':FWHM})

        # Subtract the fitted Gaussian from the raw or intermediate data and repeat
        # the peak detection step.
        GaussY = gaussFunction(x_data,*coeff)
        new_y = map(operator.sub, y_data, GaussY)
        if max(new_y) == max(y_data):
            break
        y_data = new_y
    functions = sorted(functions, key=lambda tup: tup['Peak'])

    # iterate over all peaks and remove overlap
    overlapDetected = False
    for index, i in enumerate(functions):
        try:
            if i['Data'][-1][0] > functions[index+1]['Data'][0][0]:
                overlapDetected = True
                overlap = abs(functions[index+1]['Data'][0][0] - i['Data'][-1][0])
                peak1 = max([x[1] for x in i['Data']])
                peak2 = max([x[1] for x in functions[index+1]['Data']])
                peak1fraction = (peak1 / (peak1 + peak2))*overlap
                peak2fraction = (peak2 / (peak1 + peak2))*overlap
                low = bisect.bisect_right([x[0] for x in i['Data']],i['Data'][-1][0]-peak2fraction)
                high = bisect.bisect_left([x[0] for x in functions[index+1]['Data']],functions[index+1]['Data'][0][0]+peak1fraction)
                i['Data'] = i['Data'][0:low]
                functions[index+1]['Data'] = functions[index+1]['Data'][high:-1]
        except IndexError:
            pass

    # Determine calibrants
    calibrants = determineCalibrants(functions) 

    # Writing to temp folder
    with open('temp/annotation.ref','w') as fw:
        fw.write("Peak\tRT\tWindow\n")
        for index, analyte in enumerate(functions):
            try:
                window = 0.5*(float(analyte['Data'][-1][0])-float(analyte['Data'][0][0]))
                center = float(analyte['Data'][0][0])+0.5*window
                fw.write(str("%.2f" % analyte['Peak'])+"\t"+str("%.2f" % center)+"\t"+str("%.2f" % window)+"\n")
            except:
                pass
    with open('temp/calibrants.ref','w') as fw:
        fw.write("Peak\tRT\tWindow\n")
        for index, analyte in enumerate(calibrants):
            try:
                window = 0.5*(float(analyte['Data'][-1][0])-float(analyte['Data'][0][0]))
                center = float(analyte['Data'][0][0])+0.5*window
                fw.write(str("%.2f" % analyte['Peak'])+"\t"+str("%.2f" % center)+"\t"+str("%.2f" % window)+"\n")
            except:
                pass

    # Plotting
    fig.clear()
    axes = fig.add_subplot(111)
    axes.plot(orig_x,orig_y, 'b',alpha=0.5)
    for index,i in enumerate(functions):
        try:
            xd,yd = zip(*i['Data'])
            axes.plot(xd,yd,label=str(index+1)+": "+str("%.2f" % i['Peak']))
            axes.fill_between(xd, 0, yd,alpha=0.2)
        except ValueError:
            pass
    for index,i in enumerate(calibrants):
        try:
            xd,yd = zip(*i['Data'])
            axes.annotate('Cal: '+str(index), xy=(xd[yd.index(max(yd))], max(yd)), xytext=(xd[yd.index(max(yd))], max(yd)), 
            arrowprops=dict(facecolor='black', shrink=0.05))
        except ValueError:
            pass
    axes.set_xlabel("Time [m]")
    axes.set_ylabel("Intensity [au]")
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles,labels)
    canvas.draw()

    # Warn (if needed)
    if overlapDetected == True:
        tkMessageBox.showinfo("Peak Overlap","HappyTools detected overlap between several automatically "+
        "detected peaks. HappyTools has attempted to automatically re-adjust the borders to capture the "+
        "largest possible portion of the analytes, based on their signal intensities. However, please feel "+
        "free to manually re-adjust the signals if desired in the peak list.")

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
    peaks = getPeakList(peakList)
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
        if backgroundNoiseMethod == "NOBAN":
            NOBAN = noban(intensity[lowBackground:highBackground])
        elif backgroundNoiseMethod == "MT":
            NOBAN = noban(intensity[lowBackground:highBackground])
        signalNoise = (max(intensity[low:high])-NOBAN['Background'])/NOBAN['Noise']
        # Get background subtracted peak area
        for index,j in enumerate(intensity[low:high]):
            try:
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
                dataBuffer.append((float(chunks[0]),float(chunks[1])))
        data.append((name, dataBuffer))
    return data

def saveCalibrants(fig, canvas):
    """ TODO
        Add correct try/except handling
    """
    origin = os.path.join(str(os.getcwd()),"temp","calibrants.ref")
    target = tkFileDialog.asksaveasfile(mode='w', defaultextension=".ref")
    shutil.copyfile(origin,target.name)

def saveAnnotation(fig, canvas):
    """ TODO
        Add correct try/except handling
    """
    origin = os.path.join(str(os.getcwd()),"temp","annotation.ref")
    target = tkFileDialog.asksaveasfile(mode='w', defaultextension=".ref")
    shutil.copyfile(origin,target.name)

def saveChrom():
    """ TODO
    """
    data = readData()
    saveFile = tkFileDialog.asksaveasfilename()
    with open(saveFile,'w') as fw:
        for i in data[0][1]:
            fw.write(str(format(i[0],'0.'+str(decimalNumbers)+'f'))+"\t"+str(format(i[1],'0.'+str(decimalNumbers)+'f'))+"\n")

def settingsPopup():
    """ TODO: Redesign settings window (it's fugly)
    """

    figureVariable = StringVar()
    figureVariable.set(createFigure)
    peakDetectionEdgeVar = StringVar()
    peakDetectionEdgeVar.set(peakDetectionEdge)
    ultraPerformanceCalibrationVar = StringVar()
    ultraPerformanceCalibrationVar.set(useUPC)
    options = ["True", "False"]

    def close():
        """ TODO
        """
        global points
        global start
        global end
        global baselineOrder
        global backgroundWindow
        global backgroundNoise
        #global nobanStart
        global slicepoints
        global createFigure
        global minPeaks
        global minPeakSN
        global peakDetectionMin
        global peakDetectionEdge
        global peakDetectionEdgeValue
        global useUPC
        points = int(pointsWindow.get())
        start = float(startWindow.get())
        end = float(endWindow.get())
        baselineOrder = int(baselineOrderWindow.get())
        backgroundWindow = float(backgroundWindowWindow.get())
        #nobanStart = float(nobanWindow.get())
        slicepoints = int(slicepointsWindow.get())
        createFigure = str(figureVariable.get())
        minPeaks = int(minPeakWindow.get())
        minPeakSN = int(minPeakSNWindow.get())
        peakDetectionMin = float(peakDetection.get())
        peakDetectionEdge = str(peakDetectionEdgeVar.get())
        peakDetectionEdgeValue = float(peakDetectionEdgeValueWindow.get())
        useUPC = str(ultraPerformanceCalibrationVar.get())
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
            #fw.write("nobanStart:\t"+str(float(nobanWindow.get()))+"\n")
            fw.write("noise:\t"+str(noise)+"\n")
            fw.write("slicepoints:\t"+str(slicepoints)+"\n")
            fw.write("createFigure:\t"+str(figureVariable.get())+"\n")
            fw.write("minPeaks:\t"+str(minPeakWindow.get())+"\n")
            fw.write("minPeakSN:\t"+str(minPeakSNWindow.get())+"\n")
            fw.write("peakDetectionMin:\t"+str(peakDetection.get())+"\n")
            fw.write("peakDetectionEdge:\t"+str(peakDetectionEdgeVar.get())+"\n")
            fw.write("peakDetectionEdgeValue:\t"+str(peakDetectionEdgeValueWindow.get())+"\n")
            fw.write("useUPC:\t"+str(ultraPerformanceCalibrationVar.get())+"\n")

    top = Tk.top = Toplevel()
    top.title("HappyTools "+str(HappyTools.version)+" Settings")
    top.protocol( "WM_DELETE_WINDOW", lambda: close())

    # General Settings
    generalLabel = Label(top, text="General Settings", font=("Helvetica", 16))
    generalLabel.grid(row=0, columnspan=2, sticky=W)
    
    startLabel = Label(top, text="Start Time", font=("Helvetica", 12))
    startLabel.grid(row=1, column=0, sticky=W)
    startWindow = Entry(top)
    startWindow.insert(0, start)
    startWindow.grid(row=1, column=1, sticky=W)
    
    endLabel = Label(top, text="End Time", font=("Helvetica", 12))
    endLabel.grid(row=2, column=0, sticky=W)
    endWindow = Entry(top)
    endWindow.insert(0, end)
    endWindow.grid(row=2, column=1, sticky=W)
    
    # Peak Detection Settings
    peakDetectionLabel = Label(top, text="Peak Detection Settings", font=("Helvetica", 16))
    peakDetectionLabel.grid(row=3, columnspan=2, sticky=W)

    peakDetectionLabel = Label(top, text="Minimum Intensity", font=("Helvetica", 12))
    peakDetectionLabel.grid(row=4, column=0, sticky=W)
    peakDetection = Entry(top)
    peakDetection.insert(0, peakDetectionMin)
    peakDetection.grid(row=4, column=1, sticky=W)

    peakDetectionEdgeLabel = Label(top, text="Edge Method", font=("Helvetica", 12))
    peakDetectionEdgeLabel.grid(row=5, column=0, sticky=W)
    peakDetectionEdgeWindow = OptionMenu(top, peakDetectionEdgeVar, "Sigma", "FWHM")
    peakDetectionEdgeWindow.grid(row=5, column=1, sticky=W)

    peakDetectionEdgeValueLabel = Label(top, text="Sigma Value", font=("Helvetica", 12))
    peakDetectionEdgeValueLabel.grid(row=6, column=0, sticky=W)
    peakDetectionEdgeValueWindow = Entry(top)
    peakDetectionEdgeValueWindow.insert(0, peakDetectionEdgeValue)
    peakDetectionEdgeValueWindow.grid(row=6, column=1, sticky=W)

    # Calibration Settings
    calibrationLabel = Label(top, text="Calibration Settings", font=("Helvetica", 16))
    calibrationLabel.grid(row=7, columnspan=2, sticky=W)

    minPeakLabel = Label(top, text="Minimum Peaks", font=("Helvetica", 12))
    minPeakLabel.grid(row=8, column=0, sticky=W)
    minPeakWindow = Entry(top)
    minPeakWindow.insert(0, minPeaks)
    minPeakWindow.grid(row=8, column=1, sticky=W)

    minPeakSNLabel = Label(top, text="Minimum S/N", font=("Helvetica", 12))
    minPeakSNLabel.grid(row=9, column=0, sticky=W)
    minPeakSNWindow = Entry(top)
    minPeakSNWindow.insert(0, minPeakSN)
    minPeakSNWindow.grid(row=9, column=1, sticky=W)
    
    ultraPerformanceCalibrationLabel = Label(top, text="UltraPerformanceCalibration", font=("Helvetica", 12))
    ultraPerformanceCalibrationLabel.grid(row=10, column=0, sticky=W)
    ultraPerformanceCalibrationWindow = OptionMenu(top, ultraPerformanceCalibrationVar, *options)
    ultraPerformanceCalibrationWindow.grid(row=10, column=1, sticky=W)
    
    # Quantitation Settings
    quantitationLabel = Label(top, text="Quantitation Settings", font=("Helvetica", 16))
    quantitationLabel.grid(row=11, columnspan=2, sticky=W)
    
    pointsLabel = Label(top, text="Datapoints", font=("Helvetica", 12))
    pointsLabel.grid(row=12, column=0, sticky=W)
    pointsWindow = Entry(top)
    pointsWindow.insert(0, points)
    pointsWindow.grid(row=12, column=1, sticky=W)

    baselineOrderLabel = Label(top, text="Baseline Order", font=("Helvetica", 12))
    baselineOrderLabel.grid(row=13, column=0, sticky=W)
    baselineOrderWindow = Entry(top)
    baselineOrderWindow.insert(0, baselineOrder)
    baselineOrderWindow.grid(row=13, column=1, sticky=W)
    
    backgroundWindowLabel = Label(top, text="Background Window", font=("Helvetica", 12))
    backgroundWindowLabel.grid(row=14, column=0, sticky=W)
    backgroundWindowWindow = Entry(top)
    backgroundWindowWindow.insert(0, backgroundWindow)
    backgroundWindowWindow.grid(row=14, column=1, sticky=W)

    #nobanLabel = Label(top, text="NOBAN Start", font=("Helvetica", 12))
    #nobanLabel.grid(row=14, column=0, sticky=W)
    #nobanWindow = Entry(top)
    #nobanWindow.insert(0, nobanStart)
    #nobanWindow.grid(row=14, column=1, sticky=W)

    slicepointsLabel = Label(top, text="MT Slice points", font=("Helvetica", 12))
    slicepointsLabel.grid(row=15, column=0, sticky=W)
    slicepointsWindow = Entry(top)
    slicepointsWindow.insert(0, slicepoints)
    slicepointsWindow.grid(row=15, column=1, sticky=W)

    figureLabel = Label(top, text="Create figure for each analyte", font=("Helvetica", 12))
    figureLabel.grid(row=16, column=0, sticky=W)
    #options = ["True", "False"]
    figureWindow = OptionMenu(top, figureVariable, *options)
    figureWindow.grid(row=16, column=1, sticky=W)

    # Close/Save Buttons
    saveButton = Button(top, text="Save", command=lambda: save())
    saveButton.grid(row=17, column=0, sticky=W)
    closeButton = Button(top, text="Close", command=lambda: close())
    closeButton.grid(row=17, column=1, sticky=E)

    # Tooltips
    createToolTip(pointsLabel,"The number of data points that is used to determine the baseline. Specifically, "+
            "this setting specifies how large each segment of the whole chromatogram will be to identify the lowest "+
            "data point per window, i.e. a setting of 100 means that the chromatogram is split into segments of 100 "+
            "data points per segment.")
    createToolTip(startLabel,"This setting tells the program from which time point it is supposed to begin "+
            "processing, this setting should be set in such a way that it is before the analytes of interest but "+
            "after any potential big increases or decrease in intensity.")
    createToolTip(endLabel,"This setting tells the program until which time point it is supposed to begin processing, "+
            "this setting should be set in such a way that it is after the analytes of interest but before any "+
            "potential big increases or decrease in intensity.")
    createToolTip(baselineOrderLabel,"This setting tells the program what sort of function should be used to correct "+
            "the baseline. A value of 1 refers to a linear function, while a value of 2 refers to a quadratic "+
            "function. We advise to use a linear function as the need for any higher order function indicates an "+
            "unexpected event in the chromatography.")
    createToolTip(backgroundWindowLabel,"This setting tells the program the size of the region that will be examined "+
            "to determine the background. A value of 1 means that the program will look from 20.0 to 21.0 minutes "+
            "and 21.4 to 22.4 for an analyte that elutes from 21.0 to 21.4 minutes.")
    #createToolTip(nobanLabel,"This setting specifies the initial estimate for the NOBAN algorithm, specifically a "+
            #"value of 0.25 means that the lowest 25% of all data points will be used as an initial estimate for the "+
            #"background. This value should be changed depending on how many signals there are in the chromatogram, "+
            #"e.g. in a crowded chromatogram this value should be low.")
    createToolTip(slicepointsLabel,"The number of conscutive data points that will be used to determine the "+
            "background and noise using the MT method. The MT method will scan all datapoints that fall within the "+
            "background window (specified above) to find the here specified number of consecutive data points that "+
            "yield the lowest average intensity, the average of these data points is then used as background while "+
            "the standard deviation of these data points is used as the noise.")
    createToolTip(figureLabel,"This setting specifies if HappyTools should create a figure for each integrated peak, "+
            "showing the raw datapoints, background, noise, S/N and GPQ values. This is a very performance intensive "+
            "option and it is recommended to only use this on a subset of your samples (e.g. less than 25 samples).")
    createToolTip(minPeakLabel,"This setting specifies the minimum number of calibrant peaks that have to pass the "+
            "specified S/N value that must be present in a chromatogram. A chromatogram for which there are not "+
            "enough calibrant peaks passing the specified criteria will not be calibrated and excluded from further "+
            "quantitation.")
    createToolTip(minPeakSNLabel,"This setting specifies the minimum S/N value a calibrant peak must surpass to be "+
            "included in the calibration. The actual S/N value that is determined by HappyTools depends heavily on "+
            "which method to determine signal and noise is used, the default method being rather conservative.")
    createToolTip(peakDetectionLabel, "This setting specifies the minimum intensity, relative to the main peak in "+
            "a chromatogram, that the peak detection algorithm will try to annotate. For example, a value of 0.01 "+
            "means that the program will attempt to annotate peaks until the next highest peak is below 1% of the "+
            "intensity of the main peak in the chromatogram.")
    createToolTip(peakDetectionEdgeLabel, "This setting specifies if HappyTools will determine the integration window "+
            "using either the full width at half maximum (FWHM) or a specified sigma value. The Sigma value has to be "+
            "specified in the field below if Sigma is the selected method.")
    createToolTip(peakDetectionEdgeValueLabel, "This setting specifies the Sigma value that will be used for "+
            "determining the border of the integration window. A value of 1.0 means that HappyTools will set the "+
            "integration window so that 68.3% of the Gaussian peak will be quantified (2 Sigma = 95.5% and 3 sigma "+
            "= 99.7%). Please note that this value should depend on how complex the chromatogram is, for instance "+
            "a low sigma will yield better results in a complex chromatogram.") 
    createToolTip(ultraPerformanceCalibrationLabel, "This setting specifies if HappyTools will use the standard "+
            "calibration method using a second degree polynomial (False) or if HappyTools is allowed to use the "+
            "experimental ultraPerformanceCalibration (True) that attempts to identify the optimal function for "+
            "calibration based on the residual time error of the calibrants.\n\nNOTE: The ultraPerformanceCalibration "+
            "is not ready to be used in a production environment yet!") 

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

def ultraPerformanceCalibration(measured, expected, minimum, maximum):
    """ This function tries various calibration methods, starting with 
    polynomials, followed by a power law function and lastly tries two
    interpolation methods (Pchip and Akima1D). The latter methods should
    always return an RMS of 0, which is why the function will only use
    those if they give a significant improvement in RMS (defined in
    min_improvement) and if the user has selected to use interpolation
    methods as well (defined in use_interpolation).
    
    INPUT1: List of measured data points
    INPUT2: List of expected data points
    OUTPUT: Function object
    """
    RMS = sys.maxint
    func = None

    # Polynomials between 1 and len(expected)
    for i in range(1,len(expected)):
        z = np.polyfit(measured, expected, i)
        f = np.poly1d(z)
        
        # Check if the fitted polynomial is monotone
        X_range = np.linspace(minimum, maximum, 10000)
        dx = np.diff(f(X_range))
        if not np.all(dx > 0):
            break

        # Calculate RMS
        RMS_buffer = []
        for index, j in enumerate(measured):
            RMS_buffer.append((f(j) - expected[index])**2)
        RMS_buffer = np.mean(RMS_buffer)
        RMS_buffer = math.sqrt(RMS_buffer)
        if RMS_buffer < RMS - min_improvement*RMS:
            RMS = RMS_buffer
            func = f


    # Power Law
    z = curve_fit(powerLaw, measured, expected)

    # Check if the fitted power law function is monotone
    X_range = np.linspace(minimum, maximum, 10000)
    dx = np.diff(powerLaw(X_range, *z[0]))
    if np.all(dx > 0):
        # Calculate RMS
        RMS_buffer = []
        for index, i in enumerate(measured):
            RMS_buffer.append((powerLaw(i,*z[0])-expected[index])**2)
        RMS_buffer = np.mean(RMS_buffer)
        RMS_buffer = math.sqrt(RMS_buffer)
        if RMS_buffer < RMS - min_improvement*RMS:
            RMS = RMS_buffer
            func = powerLawCall(0,*z[0])

    if use_interpolation == True:
        # Monotonic Piecewise Cubic Hermite Interpolating Polynomial
        RMS_buffer = []
        f = PchipInterpolator(measured, expected)
        for index, j in enumerate(measured):
            RMS_buffer.append((f(j)-expected[index])**2)
        RMS_buffer = np.mean(RMS_buffer)
        RMS_buffer = math.sqrt(RMS_buffer)
        if RMS_buffer < RMS - min_improvement*RMS:
            RMS = RMS_buffer
            func = f

        # Akima 1D Interpolator
        RMS_buffer = []
        f = Akima1DInterpolator(measured, expected)
        for index, j in enumerate(measured):
            RMS_buffer.append((f(j)-expected[index])**2)
        RMS_buffer = np.mean(RMS_buffer)
        RMS_buffer = math.sqrt(RMS_buffer)
        if RMS_buffer < RMS - min_improvement*RMS:
            RMS = RMS_buffer
            func = f

    return func

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
            fw.write(str(format(i[0],'0.'+str(decimalNumbers)+'f'))+"\t"+str(format(i[1],'0.'+str(decimalNumbers)+'f'))+"\n")
