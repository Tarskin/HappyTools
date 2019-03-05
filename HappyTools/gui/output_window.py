import tkinter as tk


class OutputWindow(object):
    def __init__(self, master):
        """Create a pop-up enabling output selection.

        This function creates a pop up box that allows the user to specify
        what output should be shown in the final summary. The default value
        for all variables is off (0) and by ticking a box it is set to on
        (1).

        Keyword arguments:
        none
        """
        top = tk.Toplevel()
        top.protocol('WM_DELETE_WINDOW', self.close_output_settings)
        top.title('Output Options')

        selAll = tk.Button(top, text='Select All',
                           command=self.select_all_output_settings)
        selAll.grid(row=0, column=0, sticky=tk.W)
        none = tk.Button(top, text='Select None',
                         command=self.select_no_output_settings)
        none.grid(row=0, column=1, sticky=tk.E)

        text1 = tk.Label(top, text='Base Outputs', font='bold')
        text1.grid(row=1, column=0, sticky=tk.W)
        text2 = tk.Label(top, text='Output Modifiers', font='bold')
        text2.grid(row=1, column=1, sticky=tk.W)

        ai = tk.Checkbutton(
                top, text=u'Analyte Intensity\u00B9',
                variable=master.output_parameters.absolute_intensity,
                onvalue=1, offvalue=0)
        ai.grid(row=2, column=0, sticky=tk.W)
        bck = tk.Checkbutton(
                top, text=u'\u00B9Background subtracted Intensities',
                variable=master.output_parameters.background_subtraction,
                onvalue=1, offvalue=0)
        bck.grid(row=2, column=1, sticky=tk.W)

        ri = tk.Checkbutton(
                top, text=u'Relative Intensity\u00B9',
                variable=master.output_parameters.relative_intensity,
                onvalue=1, offvalue=0)
        ri.grid(row=3, column=0, sticky=tk.W)
        exp = tk.Checkbutton(
                top, text=u'\u00B9Gaussian Intensities',
                variable=master.output_parameters.gaussian_intensity,
                onvalue=1, offvalue=0)
        exp.grid(row=3, column=1, sticky=tk.W)

        pq = tk.Checkbutton(
                top, text='Peak Quality Criteria',
                variable=master.output_parameters.analyte_quality_criteria,
                onvalue=1, offvalue=0)
        pq.grid(row=4, column=0, sticky=tk.W)

        bn = tk.Checkbutton(
                top, text='Background and Noise',
                variable=master.output_parameters.background_and_noise,
                onvalue=1, offvalue=0)
        bn.grid(row=5, column=0, sticky=tk.W)

        pdf = tk.Checkbutton(
                top, text='Generate PDF Reports',
                variable=master.output_parameters.pdf_report,
                onvalue=1, offvalue=0)
        pdf.grid(row=6, column=0, sticky=tk.W)

        ok = tk.Button(top,text = 'Ok', command=
                       self.close_output_settings)
        ok.grid(row=7, column=0, sticky = tk.W)
        save = tk.Button(top, text='Save', command=
                         self.save_output_settings)
        save.grid(row=7, column=1, sticky=tk.E)

        top.lift()

        self.master = master
        self.top = top

    def close_output_settings(self):
        self.top.destroy()

    def save_output_settings(self):
        self.master.output_parameters.save_to_disk()

    def select_all_output_settings(self):
        """Set all variables to on (1).
        """
        self.master.output_parameters.absolute_intensity.set(1)
        self.master.output_parameters.relative_intensity.set(1)
        self.master.output_parameters.gaussian_intensity.set(1)
        self.master.output_parameters.background_subtraction.set(1)
        self.master.output_parameters.background_and_noise.set(1)
        self.master.output_parameters.analyte_quality_criteria.set(1)
        self.master.output_parameters.pdf_report.set(1)

    def select_no_output_settings(self):
        """Set all variables to off (0).
        """
        self.master.output_parameters.absolute_intensity.set(0)
        self.master.output_parameters.relative_intensity.set(0)
        self.master.output_parameters.gaussian_intensity.set(0)
        self.master.output_parameters.background_subtraction.set(0)
        self.master.output_parameters.background_and_noise.set(0)
        self.master.output_parameters.analyte_quality_criteria.set(0)
        self.master.output_parameters.pdf_report.set(0)
