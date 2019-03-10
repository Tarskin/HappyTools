from HappyTools.util.functions import read_peak_list
from HappyTools.gui.progress_bar import ProgressBar
from HappyTools.bin.chromatogram import Chromatogram
from HappyTools.util.output import Output
from pathlib import Path
import logging

class BatchProcess(object):
    def __init__(self, master):
        self.process_parameters = master.process_parameters
        self.output_parameters = master.output_parameters
        self.settings = master.settings
        self.axes = master.axes
        self.logger = master.logger

    def batch_process(self):
        if self.process_parameters.data_folder:
            process_window = ProgressBar(self)

            # Calibration
            if self.process_parameters.calibration_file:
                self.reference = read_peak_list(
                        self.process_parameters.calibration_file)
                self.files = self.get_calibration_files()

                # Read data
                self.read_data()

                # Perform calibration
                for index, self.chrom in enumerate(self.data):
                    self.chrom.determine_calibration_timepairs()
                    # Still need to include a check against number of calibrants
                    self.chrom.determine_calibration_function()
                    self.chrom.calibrate_chromatogram()
                    process_window.update_progress_bar(
                            process_window.progressbar,
                            process_window.calibration_percentage,
                            index, len(self.data))
                process_window.update_progress_bar(
                        process_window.progressbar,
                        process_window.calibration_percentage, 1, 1)

            # Quantitation
            if self.process_parameters.quantitation_file:
                self.reference = read_peak_list(
                        self.process_parameters.quantitation_file)
                self.files = self.get_quantitation_files()

                # Read data
                if not self.data:
                    self.data = self.read_data()

                # Perform quantitation
                for index, self.chrom in enumerate(self.data):
                    self.chrom.quantify_chromatogram()
                    process_window.update_progress_bar(
                            process_window.progressbar2,
                            process_window.quantitation_percentage,
                            index, len(self.data))

                # Generate summary file
                output = Output(self)
                output.init_output_file()
                output.build_output_file()

                process_window.update_progress_bar(
                        process_window.progressbar2,
                        process_window.quantitation_percentage, 1, 1)

    def read_data(self):
        data = []
        for index, filename in enumerate(self.files):
            self.filename = Path(filename)
            chromatogram = Chromatogram(self)
            chromatogram.open_chromatogram()
            data.append(chromatogram)
        self.data = data
        
    def get_calibration_files(self):
        calibration_files = []
        for files in self.settings.calibration_filetypes:
            for file in Path(self.process_parameters.data_folder).glob(files):
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
