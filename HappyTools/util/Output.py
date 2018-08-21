import HappyTools.gui.version as version
from datetime import datetime
from os import path

class Output(object):
    def __init__(self, master):
        self.abs_int = master.abs_int
        self.rel_int = master.rel_int
        self.gauss_int = master.gauss_int
        self.bck_sub = master.rel_int
        self.bck_noise = master.rel_int
        self.peak_qual = master.rel_int
        self.settings = master.settings
        self.results = master.results
        self.batch_folder = master.batch_folder

        utc_datetime = datetime.utcnow()
        s = utc_datetime.strftime(master.settings.date_format)
        self.filename = s + "_" + master.settings.output

    def build_output_file(self, master):
        header = ""
        for i in master.results:
            for j in i['results']:
                header = header + "\t"+str(j['peak'])
            header = header + "\n"
            for j in i['results']:
                header = header + "\t"+str(j['time'])
            header = header + "\n"
            break
        self.header = header

        if master.abs_int.get() == 1:
            if master.bck_sub.get() == 0:
                self.write_non_back_sub_abs_peak_area(self)
            elif master.bck_sub.get() == 1:
                self.write_back_sub_abs_peak_area(self)

        elif master.rel_int.get() == 1:
            if master.bck_sub.get() == 0:
                self.write_non_back_sub_rel_peak_area(self)
            elif master.bck_sub.get() == 1:
                self.write_back_sub_rel_peak_area(self)

        if master.gauss_int.get() == 1:
            if master.bck_sub.get() == 0:
                self.write_non_back_sub_gauss_peak_area(self)
            elif master.bck_sub.get() == 1:
                self.write_back_sub_gauss_peak_area(self)

        if master.bck_noise.get() == 1:
            self.write_peak_noise(self)
            self.write_background(self)
            self.write_noise(self)

        if master.peak_qual.get() == 1:
            self.write_gaussian_pattern_quality(self)
            self.write_fwhm(self)
            self.write_signal_noise(self)
            self.write_retention_time(self)
            self.write_retention_time_residual(self)

    def init_output_file(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'w') as fw:
            fw.write("HappyTools Settings\n")
            fw.write("Version:\t"+str(version.version)+"\n")
            fw.write("Build:\t"+str(version.build)+"\n")
            fw.write("Start Time:\t"+str(master.settings.start)+"\n")
            fw.write("End Time:\t"+str(master.settings.end)+"\n")
            fw.write("Baseline Order:\t"+str(master.settings.baseline_order)+"\n")
            fw.write("Background Window:\t"+str(master.settings.background_window)+"\n")
            fw.write("Background and noise method:\t"+str(master.settings.background_noise_method)+"\n")
            if master.settings.background_noise_method == "MT":
                fw.write("MT Slice Points:\t"+str(master.settings.slicepoints)+"\n")
            elif master.settings.background_noise_method == "NOBAN":
                fw.write("NOBAN Initial Estimate:\t"+str(master.settings.noban_start)+"\n")
            fw.write("Noise:\t"+str(master.settings.noise)+"\n")
            fw.write("\n")

    def write_back_sub_abs_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Peak Area (Background Subtracted)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['peak_area']-j['background_area']))
                fw.write("\n")
            fw.write("\n")

    def write_back_sub_gauss_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Gaussian Area (Background Subtracted)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['gaussian_area']-j['background_area']))
                fw.write("\n")
            fw.write("\n")

    def write_back_sub_rel_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Relative Peak Area (TAN, Background Subtracted)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                total = 0.
                for j in i['results']:
                    total += max(j['peak_area']-j['background_area'],0)
                for j in i['results']:
                    try:
                        fw.write("\t"+str(max(j['peak_area']-j['background_area'],0)/total))
                    except ZeroDivisionError:
                        fw.write("\t"+str(0.0))
                fw.write("\n")
            fw.write("\n")

    def write_background(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Background")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['background']))
                fw.write("\n")
            fw.write("\n")

    def write_fwhm(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("FWHM")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['fwhm']))
                fw.write("\n")
            fw.write("\n")

    def write_gaussian_pattern_quality(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("GPQ (Gaussian Peak Quality)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['residual']))
                fw.write("\n")
            fw.write("\n")

    def write_noise(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Noise")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['noise']))
                fw.write("\n")
            fw.write("\n")

    def write_non_back_sub_abs_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Peak Area")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['peak_area']))
                fw.write("\n")
            fw.write("\n")

    def write_non_back_sub_gauss_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Gaussian Area")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['gaussian_area']))
                fw.write("\n")
            fw.write("\n")

    def write_non_back_sub_rel_peak_area(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Relative Peak Area (TAN)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                total = 0.
                for j in i['results']:
                    total += j['peak_area']
                for j in i['results']:
                    try:
                        fw.write("\t"+str(j['peak_area']/total))
                    except ZeroDivisionError:
                        fw.write("\t"+str(0.0))
                fw.write("\n")
            fw.write("\n")

    def write_peak_noise(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Peak Noise (Standard deviation of integration window)")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['peak_noise']))
                fw.write("\n")
            fw.write("\n")

    def write_retention_time(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Retention Time")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(float(j['actual_time'])))
                fw.write("\n")
            fw.write("\n")

    def write_retention_time_residual(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Retention Time Residual")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(abs(j['actual_time'] - j['time'])))
                fw.write("\n")
            fw.write("\n")

    def write_signal_noise(self, master):
        with open(path.join(master.batch_folder.get(), self.filename), 'a') as fw:
            fw.write("Signal-to-Noise")
            fw.write(self.header)
            for i in master.results:
                fw.write(i['file'])
                for j in i['results']:
                    fw.write("\t"+str(j['signal_noise']))
                fw.write("\n")
            fw.write("\n")
