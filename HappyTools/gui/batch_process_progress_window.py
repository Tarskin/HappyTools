from HappyTools.gui.progress_bar import ProgressBar
import tkinter as tk


class BatchProcessProgressWindow(object):
    def __init__(self, master):
        self.master = master
        self.process_parameters = master.process_parameters
        self.output_parameters = master.output_parameters

    def create_window(self):
        top = tk.Toplevel()
        top.maxsize(width=250, height=0)
        top.protocol('WM_DELETE_WINDOW', self.close)
        top.title('Batch Process Progress Window')

        if self.process_parameters.calibration_file or \
                self.process.quantitation_file:
            reading_label = tk.Label(top, text='Data Reading')
            reading_label.pack(fill='both', expand=True)
            reading_progress_bar = ProgressBar(top)
            reading_progress_bar.reset_bar()
            reading_progress_bar.bar.pack(fill='both', expand=True)
            self.reading_progress_bar = reading_progress_bar

        if self.process_parameters.calibration_file:
            calibration_label = tk.Label(top, text='Calibration')
            calibration_label.pack(fill='both', expand=True)
            calibration_progress_bar = ProgressBar(top)
            calibration_progress_bar.reset_bar()
            calibration_progress_bar.bar.pack(fill='both', expand=True)
            self.calibration_progress_bar = calibration_progress_bar

        if self.process_parameters.quantitation_file:
            quantitation_label = tk.Label(top, text='Quantitation')
            quantitation_label.pack(fill='both', expand=True)
            quantitation_progress_bar = ProgressBar(top)
            quantitation_progress_bar.reset_bar()
            quantitation_progress_bar.bar.pack(fill='both', expand=True)
            self.quantitation_progress_bar = quantitation_progress_bar

        if self.output_parameters.pdf_report.get() == True:
            report_label = tk.Label(top, text='PDF Report Generation')
            report_label.pack(fill='both', expand=True)
            report_progress_bar = ProgressBar(top)
            report_progress_bar.reset_bar()
            report_progress_bar.bar.pack(fill='both', expand=True)
            self.report_progress_bar = report_progress_bar

        self.top = top

    def close(self):
        self.top.destroy()
