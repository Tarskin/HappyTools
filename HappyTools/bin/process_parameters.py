from pathlib import Path


class ProcessParameters(object):
    def __init__(self, master):
        self.calibration_file = None
        self.calibration = False
        self.quantitation_file = None
        self.quantitation = False
        self.data_folder = Path.cwd()
