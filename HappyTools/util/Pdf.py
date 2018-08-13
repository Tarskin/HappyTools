import HappyTools.gui.Version as version
import matplotlib.pyplot as plt
from HappyTools.util.Functions import Functions

from matplotlib.backends.backend_pdf import PdfPages
from bisect import bisect_left, bisect_right
from numpy import linspace
from datetime import datetime
from scipy.interpolate import InterpolatedUnivariateSpline
from os import path

class Pdf(object):
    def __init__(self, master):
        self.pdf = PdfPages(path.join(master.master.batch_folder.get(),path.splitext(path.basename(master.data.filename))[0]+".pdf"))
    
    def plot_overview(self, master):
        time, intensity = zip(*master.data.data)
        d = self.pdf.infodict()
        d['Title'] = 'PDF Report for: '+str(path.splitext(path.basename(master.data.filename))[0])
        d['Author'] = 'HappyTools version: '+str(version.version)+" build: "+str(version.build)
        d['CreationDate'] = datetime.now()
        low = bisect_left(time, master.settings.start)
        high = bisect_right(time, master.settings.end)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        plt.plot(time[low:high], intensity[low:high], 'b-')
        plt.legend(['Raw Data'], loc='best')
        plt.title(str(path.splitext(path.basename(master.data.filename))[0]))
        plt.xlabel("Retention Time [m]")
        plt.ylabel("Intensity [au]")
        for i in master.master.reference:
            low = bisect_left(time, i[1]-i[2])
            high = bisect_right(time, i[1]+i[2])
            new_time = linspace(time[low], time[high], len(time[low:high]))
            f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])
            new_intensity = f(new_time)
            ax.fill_between(time[low:high], 0, new_intensity, alpha=0.5)
            ax.text(i[1], max(intensity[low:high]), i[0])
        self.pdf.savefig(fig)
        plt.close(fig)

    def plot_individual(self, master):
        time, intensity = zip(*master.data.data)
        f = InterpolatedUnivariateSpline(time[low:high], intensity[low:high])

        new_x = np.linspace(time[master.low], time[master.high], 2500*(time[master.high]-time[master.low]))
        new_y = f(new_x)

        new_gauss_x = np.linspace(time[master.low], time[master.high], 2500*(time[master.high]-time[master.low]))
        new_gauss_y = Functions().gauss_function(new_gauss_x, *master.coeff)

        fig =  plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        plt.plot(time[master.low:master.high], intensity[master.low:master.high], 'b*')
        plt.plot((new_x[0],new_x[-1]),(master.NOBAN['Background'],master.NOBAN['Background']),'red')
        plt.plot((new_x[0],new_x[-1]),(master.NOBAN['Background']+master.NOBAN['Noise'],master.NOBAN['Background']+master.NOBAN['Noise']),color='green')
        plt.plot(new_x,new_y, color='blue',linestyle='dashed')
        plt.plot(new_gauss_x, new_gauss_y, color='green',linestyle='dashed')
        plt.plot((time[intensity[master.low:master.high].index(max(intensity[master.low:master.high]))+master.low],
                time[intensity[master.low:master.high].index(max(intensity[master.low:master.high]))+master.low]),
                (master.NOBAN['Background'],max(intensity[master.low:master.high])),
                color='orange',linestyle='dotted')
        plt.plot((min(max(master.fwhm['center']-master.fwhm['width'],new_x[0]),new_x[-1]),max(min(master.fwhm['center']+master.fwhm['width'],new_x[-1]),new_x[0])),
                (height,height),color='red',linestyle='dashed')
        plt.legend(['Raw Data','Background','Noise','Univariate Spline','Gaussian Fit ('+str(int(master.residual*100))+
                '%)','Signal (S/N '+str(master.signal_noise)+")","FWHM:"+"{0:.2f}".format(master.fwhm['fwhm'])], loc='best')
        plt.title("Detail view: "+str(master.peak))
        plt.xlabel("Retention Time [m]")
        plt.ylabel("Intensity [au]")
        self.pdf.savefig(fig)
        plt.close(fig)

    def close(self, master):
        self.pdf.close()
