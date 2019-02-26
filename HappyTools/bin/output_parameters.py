import configparser
import logging
import tkinter as tk
from pathlib import Path


class OutputParameters(object):
    def __init__(self, master):
        self.config = configparser.ConfigParser()
        self.logger = logging.getLogger(__name__)

        self.absolute_intensity = tk.IntVar()
        self.relative_intensity = tk.IntVar()
        self.gaussian_intensity = tk.IntVar()
        self.background_subtraction = tk.IntVar()
        self.background_and_noise = tk.IntVar()
        self.analyte_quality_criteria = tk.IntVar()
        self.pdf_report = tk.IntVar()

        self.read_from_disk()

    def save_to_disk(self):
        try:
            self.config.add_section('Output')
        except configparser.DuplicateSectionError as e:
            self.logger.error(e)
        self.config.set('Output', 'Absolute Intensity',
                        str(self.absolute_intensity.get()))
        self.config.set('Output', 'Relative Intensity',
                        str(self.relative_intensity.get()))
        self.config.set('Output', 'Gaussian Intensity',
                        str(self.gaussian_intensity.get()))
        self.config.set('Output', 'Background Subtraction',
                        str(self.background_subtraction.get()))
        self.config.set('Output', 'Background and Noise',
                        str(self.background_and_noise.get()))
        self.config.set('Output', 'Analyte Quality Criteria',
                        str(self.analyte_quality_criteria.get()))
        self.config.set('Output', 'PDF Reports',
                        str(self.pdf_report.get()))
        with Path(Path.cwd() / 'HappyTools.ini').open('w') as config_file:
            self.config.write(config_file)

    def read_from_disk(self):
        if Path(Path.cwd() / 'HappyTools.ini').is_file():
            self.config.read('HappyTools.ini')

            try:
                self.absolute_intensity.set(int(
                    self.config.get('Output', 'Absolute Intensity')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.relative_intensity.set(int(
                    self.config.get('Output', 'Relative Intensity')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.gaussian_intensity.set(int(
                    self.config.get('Output', 'Gaussian Intensity')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.background_subtraction.set(int(
                    self.config.get('Output', 'Background Subtraction')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.background_and_noise.set(int(
                    self.config.get('Output', 'Background and Noise')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.analyte_quality_criteria.set(int(
                    self.config.get('Output', 'Analyte Quality Criteria')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass

            try:
                self.pdf_report.set(int(
                    self.config.get('Output', 'PDF Reports')))
            except (configparser.NoSectionError,
                    configparser.NoOptionError):
                pass
