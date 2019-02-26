import HappyTools.gui.version as version
from datetime import datetime
from pathlib import Path

class Output(object):
    def __init__(self, master):
        self.master = master
        self.settings = master.settings
        self.results = master.results
        self.output_parameters = master.output_parameters
        self.batch_folder = master.batch_folder

        utc_datetime = datetime.utcnow()
        s = utc_datetime.strftime(master.settings.date_format)
        self.filename = s + '_' + master.settings.output

    def build_header(self):
        header = ''
        for i in self.results:

            for j in i['results']:
                header = header + '\t'+str(j['peak'])
            header = header + '\n'

            header = header + 'RT'
            for j in i['results']:
                header = header + '\t'+str(j['time'])
            header = header + '\n'

            break
        self.header = header

    def build_output_file(self):
        self.build_header()

        if self.output_parameters.absolute_intensity.get() == 1:
            if self.output_parameters.background_subtraction.get() == 0:
                self.write_non_back_sub_abs_peak_area()
            elif self.output_parameters.background_subtraction.get() == 1:
                self.write_back_sub_abs_peak_area()

        if self.output_parameters.relative_intensity.get() == 1:
            if self.output_parameters.background_subtraction.get() == 0:
                self.write_non_back_sub_rel_peak_area()
            elif self.output_parameters.background_subtraction.get() == 1:
                self.write_back_sub_rel_peak_area()

        if self.output_parameters.gaussian_intensity.get() == 1:
            if self.output_parameters.background_subtraction.get() == 0:
                self.write_non_back_sub_gauss_peak_area()
            elif self.output_parameters.background_subtraction.get() == 1:
                self.write_back_sub_gauss_peak_area()

        if self.output_parameters.background_and_noise.get() == 1:
            self.write_peak_noise()
            self.write_background()
            self.write_noise()

        if self.output_parameters.analyte_quality_criteria.get() == 1:
            self.write_gaussian_pattern_quality()
            self.write_fwhm()
            self.write_signal_noise()
            self.write_retention_time()
            self.write_retention_time_residual()

    def init_output_file(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('w') as fw:
            fw.write('HappyTools Settings\n')
            fw.write('Version:\t'+str(version.version)+'\n')
            fw.write('Build:\t'+str(version.build)+'\n')
            fw.write('Start Time:\t'+str(self.settings.start)+'\n')
            fw.write('End Time:\t'+str(self.settings.end)+'\n')
            fw.write('Baseline Order:\t'+str(self.settings.baseline_order)+'\n')
            fw.write('Background Window:\t'+str(self.settings.background_window)+'\n')
            fw.write('Background and noise method:\t'+str(self.settings.background_noise_method)+'\n')
            if self.settings.background_noise_method == 'MT':
                fw.write('MT Slice Points:\t'+str(self.settings.slicepoints)+'\n')
            elif self.settings.background_noise_method == 'NOBAN':
                fw.write('NOBAN Initial Estimate:\t'+str(self.settings.noban_start)+'\n')
            fw.write('Noise:\t'+str(self.settings.noise)+'\n')
            fw.write('\n')

    def write_back_sub_abs_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Peak Area (Background Subtracted)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['peak_area']-j['background_area']))
                fw.write('\n')
            fw.write('\n')

    def write_back_sub_gauss_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Gaussian Area (Background Subtracted)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['gaussian_area']-j['background_area']))
                fw.write('\n')
            fw.write('\n')

    def write_back_sub_rel_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Relative Peak Area (TAN, Background Subtracted)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                total = 0.
                for j in i['results']:
                    total += max(j['peak_area']-j['background_area'],0)
                for j in i['results']:
                    try:
                        fw.write('\t'+str(max(j['peak_area']-j['background_area'],0)/total))
                    except ZeroDivisionError:
                        fw.write('\t'+str(0.0))
                fw.write('\n')
            fw.write('\n')

    def write_background(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Background')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['background']))
                fw.write('\n')
            fw.write('\n')

    def write_fwhm(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('FWHM')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['fwhm']))
                fw.write('\n')
            fw.write('\n')

    def write_gaussian_pattern_quality(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('GPQ (Gaussian Peak Quality)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['residual']))
                fw.write('\n')
            fw.write('\n')

    def write_noise(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Noise')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['noise']))
                fw.write('\n')
            fw.write('\n')

    def write_non_back_sub_abs_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Peak Area')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['peak_area']))
                fw.write('\n')
            fw.write('\n')

    def write_non_back_sub_gauss_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Gaussian Area')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['gaussian_area']))
                fw.write('\n')
            fw.write('\n')

    def write_non_back_sub_rel_peak_area(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Relative Peak Area (TAN)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                total = 0.
                for j in i['results']:
                    total += j['peak_area']
                for j in i['results']:
                    try:
                        fw.write('\t'+str(j['peak_area']/total))
                    except ZeroDivisionError:
                        fw.write('\t'+str(0.0))
                fw.write('\n')
            fw.write('\n')

    def write_peak_noise(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Peak Noise (Standard deviation of integration window)')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['peak_noise']))
                fw.write('\n')
            fw.write('\n')

    def write_retention_time(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Retention Time [min.]')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(float(j['actual_time'])))
                fw.write('\n')
            fw.write('\n')

    def write_retention_time_residual(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Retention Time Residual [min.]')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(abs(j['actual_time'] - j['time'])))
                fw.write('\n')
            fw.write('\n')

    def write_signal_noise(self):
        with Path(self.batch_folder.get() /
                  Path(self.filename)).open('a') as fw:
            fw.write('Signal-to-Noise')
            fw.write(self.header)
            for i in self.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write('\t'+str(j['signal_noise']))
                fw.write('\n')
            fw.write('\n')
