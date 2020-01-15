import configparser
import logging
from pathlib import Path


class Settings(object):
    def __init__(self, master):
        self.config = configparser.ConfigParser()
        self.logger = logging.getLogger(__name__)

        # General Settings
        self.start = 10
        self.end = 60
        self.points = 100
        self.baseline_order = 1
        self.background_window = 1

        # Peak Detection Settings
        self.peak_detection_min = 0.05
        self.peak_detection_edge = 'Sigma'
        self.peak_detection_edge_value = 2.0

        # Calibration Settings
        self.min_peaks = 4
        self.min_peak_SN = 27

        # Advanced Settings (non-GUI)
        self.encoding='utf-8_sig'
        self.noban_start = 0.25
        self.slicepoints = 5
        self.decimal_numbers = 6
        self.min_improvement = 0.05
        self.use_interpolation = False
        self.use_UPC = 'False'
        self.noise = 'RMS'
        self.background_noise_method = 'MT'

        self.read_from_disk()

        # These should go somewhere else, processing parameters?
        self.output = 'summary.results'
        self.date_format = '%Y-%m-%d-%H%MZ'
        self.settings = 'HappyTools.ini'

        self.exclusion_files = ['LICENSE.txt', 'CHANGELOG.txt']
        self.calibration_filetypes = ['*.txt', '*.arw', '*.csv']
        self.quantitation_filetypes = ['calibrated*.txt']

    def read_from_disk(self):
        if Path(Path.cwd() / 'HappyTools.ini').is_file():
            self.config.read('HappyTools.ini')

            # General Settings
            try:
                self.start = float(self.config.get(
                        'General', 'Start time'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.end = float(self.config.get(
                        'General', 'End time'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.points = int(self.config.get(
                        'General', 'Points for baseline'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.baseline_order = int(self.config.get(
                        'General', 'Baseline function order'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.background_window = float(self.config.get(
                        'General', 'Background window'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            # Peak Detection Settings
            try:
                self.peak_detection_min = float(self.config.get(
                        'Peak Detection', 'Minimum peak intensity'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.peak_detection_edge = self.config.get(
                        'Peak Detection', 'Edge type')
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.peak_detection_edge_value = float(self.config.get(
                        'Peak Detection', 'Edge value'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            # Calibration Settings
            try:
                self.min_peaks = int(self.config.get(
                        'Calibration', 'Minimum number of peaks'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.min_peak_SN = float(self.config.get(
                        'Calibration', 'Minimum sn for calibration'))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

    def save_to_disk(self):
        print ('Trying Boss!')
        try:
            self.config.add_section('General')
        except configparser.DuplicateSectionError:
            pass

        try:
            self.config.add_section('Peak Detection')
        except configparser.DuplicateSectionError:
            pass

        try:
            self.config.add_section('Calibration')
        except configparser.DuplicateSectionError:
            pass

        try:
            self.config.add_section('Quantitation')
        except configparser.DuplicateSectionError:
            pass

        # General Settings
        self.config.set('General', 'Points for baseline',
                        str(self.points))
        self.config.set('General', 'Start time',
                        str(self.start))
        self.config.set('General', 'End time',
                        str(self.end))
        self.config.set('General', 'Baseline function order',
                        str(self.baseline_order))
        self.config.set('General', 'Background window',
                        str(self.background_window))

        # Peak Detection Settings
        self.config.set('Peak Detection', 'Minimum peak intensity',
                        str(self.peak_detection_min))
        self.config.set('Peak Detection', 'Edge type',
                        str(self.peak_detection_edge))
        self.config.set('Peak Detection', 'Edge value',
                        str(self.peak_detection_edge_value))

        # Calibration Settings
        self.config.set('Calibration', 'Minimum number of peaks',
                        str(self.min_peaks))
        self.config.set('Calibration', 'Minimum sn for calibration',
                        str(self.min_peak_SN))

        # Save
        with Path(Path.cwd() / 'HappyTools.ini').open('w') as config_file:
            self.config.write(config_file)
