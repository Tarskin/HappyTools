from HappyTools.util.fitting import gauss_function
from matplotlib.figure import Figure
import HappyTools.gui.version as version

from bisect import bisect_left, bisect_right
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
from numpy import linspace
from pathlib import Path, PurePath
from scipy.interpolate import InterpolatedUnivariateSpline


class Pdf(object):
    def __init__(self, master):
        pdf_file = str(PurePath(master.filename).stem)+'.pdf'
        pdf = PdfPages(master.master.process_parameters.data_folder /
                       Path(pdf_file))

        fig = Figure(figsize=(8, 6))
        axes = fig.add_subplot(111)
        axes.set_xlabel('Retention Time [m]')
        axes.set_ylabel('Intensity [au]')

        self.settings = master.settings
        self.master = master
        self.fig = fig
        self.axes = axes
        self.pdf = pdf

    def plot_overview(self):
        time, intensity = zip(*self.master.chrom_data)
        d = self.pdf.infodict()
        d['Title'] = 'PDF Report for: '+str(
                PurePath(self.master.filename).stem)
        d['Author'] = ('HappyTools version: '+str(version.version)+
                       ' build: '+str(version.build))
        d['CreationDate'] = datetime.now()
        low = bisect_left(time, self.settings.start)
        high = bisect_right(time, self.settings.end)

        self.axes.clear()
        self.axes.plot(time[low:high], intensity[low:high], 'b-')
        self.axes.legend(['Raw Data'], loc='best')
        self.axes.set_title(str(
                PurePath(self.master.filename).stem))
        for i in self.master.master.reference:
            low = bisect_left(time, i[1]-i[2])
            high = bisect_right(time, i[1]+i[2])
            new_time = linspace(time[low], time[high],
                                len(time[low:high]))
            f = InterpolatedUnivariateSpline(time[low:high],
                                             intensity[low:high])
            new_intensity = f(new_time)
            self.axes.fill_between(time[low:high], 0, new_intensity,
                                   alpha=0.5)
            self.axes.text(i[1], max(intensity[low:high]), i[0],
                           fontsize=6, rotation=90, ha='left',
                           va='bottom')
        self.pdf.savefig(self.fig)

    def plot_individual(self):
        time, intensity = zip(*self.master.chrom_data)
        low = bisect_left(time, self.master.peak.time-
                          self.master.peak.window)
        high = bisect_right(time, self.master.peak.time+
                            self.master.peak.window)

        f = InterpolatedUnivariateSpline(time[low:high],
                                         intensity[low:high])

        new_x = linspace(time[low], time[high],
                         2500*(time[high]-time[low]))
        new_y = f(new_x)

        if self.master.peak.coeff.size > 0:
            new_gauss_x = linspace(time[low], time[high], 2500*(
                                   time[high]-time[low]))
            new_gauss_y = gauss_function(new_gauss_x,
                                         *self.master.peak.coeff)

        self.axes.clear()
        self.axes.plot(time[low:high], intensity[low:high], 'b*')
        self.axes.plot(
                (new_x[0], new_x[-1]),(self.master.peak.background,
                self.master.peak.background), 'red')
        self.axes.plot(
                (new_x[0],new_x[-1]),(self.master.peak.background+
                self.master.peak.noise,self.master.peak.background+
                self.master.peak.noise), color='green')
        self.axes.plot(new_x,new_y, color='blue',linestyle='dashed')
        if self.master.peak.coeff.size > 0:
            self.axes.plot(new_gauss_x, new_gauss_y, color='green',
                           linestyle='dashed')
        self.axes.plot(
                (time[intensity[low:high].index(max(intensity[low:high]))+low],
                time[intensity[low:high].index(max(intensity[low:high]))+low]),
                (self.master.peak.background,max(intensity[low:high])),
                color='orange',linestyle='dotted')
        self.axes.plot(
                (min(max(self.master.peak.center-self.master.peak.width,new_x[0]),
                new_x[-1]),max(min(self.master.peak.center+self.master.peak.width,
                new_x[-1]),new_x[0])), (self.master.peak.height,
                self.master.peak.height),color='red',linestyle='dashed')
        self.axes.legend(
                ['Raw Data','Background','Noise','Univariate Spline',
                'Gaussian Fit ('+str(int(self.master.peak.residual*100))+
                '%)','Signal (S/N '+'{0:.2f}'.format(
                self.master.peak.signal_noise)+')','FWHM: '+'{0:.2f}'.format(
                self.master.peak.fwhm)], loc='best')
        self.axes.set_title('Detail view: '+str(self.master.peak.peak))
        self.pdf.savefig(self.fig)

    def close(self):
        self.pdf.close()
