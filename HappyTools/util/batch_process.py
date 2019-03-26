from HappyTools.util.functions import read_peak_list
from HappyTools.gui.batch_process_progress_window import \
        BatchProcessProgressWindow
from HappyTools.bin.chromatogram import Chromatogram
from HappyTools.util.output import Output
from pathlib import Path


class BatchProcess(object):
    def __init__(self, master):
        self.master = master
        self.process_parameters = master.process_parameters
        self.output_parameters = master.output_parameters
        self.settings = master.settings
        self.axes = master.axes
        self.logger = master.logger

    def batch_process(self):
        if self.process_parameters.data_folder:
            self.process_window = BatchProcessProgressWindow(self)
            self.process_window.create_window()

            # Calibration
            if self.process_parameters.calibration_file:
                self.reference = read_peak_list(
                        self.process_parameters.calibration_file)
                self.files = self.get_calibration_files()

                # Read data
                self.read_data()

                # Perform calibration
                progress = self.process_window.calibration_progress_bar
                for index, self.chrom in enumerate(self.data):
                    progress.counter.set(
                            (float(index) / len(self.data))*100)
                    progress.update_progress_bar()
                    self.chrom.determine_calibration_timepairs()
                    # Still need to include a check against number of calibrants
                    self.chrom.determine_calibration_function()
                    self.chrom.calibrate_chromatogram()
                progress.fill_bar()

            # Quantitation
            if self.process_parameters.quantitation_file:
                self.reference = read_peak_list(
                        self.process_parameters.quantitation_file)
                self.files = self.get_quantitation_files()

                # Read data
                if not self.data:
                    self.data = self.read_data()

                # Perform quantitation
                progress = self.process_window.quantitation_progress_bar
                for index, self.chrom in enumerate(self.data):
                    progress.counter.set(
                            (float(index) / len(self.data))*100)
                    progress.update_progress_bar()
                    self.chrom.quantify_chromatogram()
                progress.fill_bar()

                # Generate summary file
                output = Output(self)
                output.init_output_file()
                output.build_output_file()

            # Report generation
            if self.output_parameters.pdf_report.get() == True:
                progress = self.process_window.report_progress_bar
                for index, self.chrom in enumerate(self.data):
                    progress.counter.set(
                            (float(index) / len(self.data))*100)
                    progress.update_progress_bar()
                    self.chrom.generate_pdf_report()
                progress.fill_bar()

    def read_data(self):
        data = []
        progress = self.process_window.reading_progress_bar
        for index, filename in enumerate(self.files):
            progress.counter.set((float(index) / len(self.files))*100)
            progress.update_progress_bar()
            self.filename = Path(filename)
            chromatogram = Chromatogram(self)
            chromatogram.open_chromatogram()
            data.append(chromatogram)
        progress.fill_bar()
        self.data = data

    def get_calibration_files(self):
        calibration_files = []
        for files in self.settings.calibration_filetypes:
            for file in Path(
                    self.process_parameters.data_folder).glob(files):
                if file not in self.settings.exclusion_files:
                    calibration_files.append(Path(
                            self.process_parameters.data_folder) /
                            file)
        return calibration_files

    def get_quantitation_files(self):
        quantitation_files = []
        for files in self.settings.quantitation_filetypes:
            for file in Path(self.process_parameters.data_folder).glob(files):
                if file not in self.settings.exclusion_files:
                    quantitation_files.append(Path(
                            self.process_parameters.data_folder) /
                            file)
        return quantitation_files
