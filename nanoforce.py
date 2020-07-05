#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File      : afm.py
# Author    : Chris Jones
# Email     : crj341@student.bham.ac.uk
# Date      : 26/06/2020
# Version   : 0.1

import  numpy                   as      np
import  plotly.graph_objects    as      go
import  easygui                 as      gui
from    scipy                   import  interpolate as  interp
import  re
import  os

class AFM:

    '''
    The nanofoce.AFM() class is used for importing, storing and analysing AFM
    force curves produced by Nanoscope 6. Details of the parameters, functions
    and an example running script are given below.

    Prerequisits
    ------------
    pip install numpy
    pip install plotly
    pip install scipy
    pip install easygui

    Parameters
    ----------  
    approach_raw: [samples, curve_number]
        Raw AFM force curve approach data (nN) imported from nanoscope
        files without baseline or contact corrections applied, produced
        by input_files function
    retract_raw: [samples, curve_number]
        Raw AFM force curve retract data (nN) imported from nanoscope
        files without baseline or contact corrections applied, produced
        by input_files function
    z_position_raw: [samples, curve_number]
        Raw AFM force curve tip position data (nm) imported from nanoscope
        files without baseline or contact corrections applied, produced
        by input_files function
    approach: [samples, curve_number]
        AFM force curve approach data (nN) produced by the baseline function
        and updated by teh contact function
    retract: [samples, curve_number]
        AFM force curve retract data (nN) produced by the baseline function
        and updated by teh contact function
    z_position:[samples, curve_number]
        AFM force curve tip position data (nm) produced by the baseline
        function and updated by teh contact function
    adhesion: [curve_number]
        Adhesion (nN) value calculated for each force curve by calc_adhesion
    modulus: [curve_number]
        Modulus (MPa) value calculated for each force curve by calc_modulus
    def_sens:
        Deflection sensetivity (nm/V), imported manually from first force
        curve file and manually overriden with set_def_sense function
    spr_const:
        Spring constant (N/m), imported manually from first force
        curve file and manually overriden with set_def_sense function
    file_name:
        First file name (usually with .000 extension)
    filelist: [curve_number]
        Full list of imported file names (.000 to .n extension)
    run_name:
        Name of experiment, specified by set_run_name function

    Functions
    ---------
    set_run_name:
        Manually set 'run_name' variable <.set_run_name(Sample name here)>
    set_def_sens:
        Manually set 'def_sens' variable <.set_def_sens(Value)>
    set_spr_const:
        Manually set 'spr_const' variable <.set_spr_const(Value)>
    run:
        Basic run script calling:
            .input_files()
            .nanoscope_params()
            .nanoscope_read()
            .baseline()
            .contact()
            .plot_curves()
            .calc_adhesion()
            .calc_modulus()
    input_files:
        Function to find input file names. If called with no inputs a GUI
        window will appear prompting the user to select the first Nanoscope
        file (usually with '.000' extension). Alternativley the file path can
        be given as an input to the function <.input_files(\path\to\file\.000)>.
        The function will then search for all subsequent files with the same
        name and consecutive file extensions (.000, .001, .002, ..., .nnn)
    nanoscope_params:
        Function to extract relevant parameters from first Nanoscope file
    nanoscope_read:
        Function to import raw force curve data
    plot_raw:
        Function to plot individual, unprocessed force curves. Curve number to
        plot must be given as an input <.plot_raw(n)>
    baseline:
        Function to adjust force curve baseline to zero by taking the mean value
        of a specified region of the approach curve. The default is from 0.45*z
        to 0.8*z, starting at the contact side of the curve. The start and end
        points may be specified using the 'start_pos' and 'end_pos' inputs to
        account for a longer/shorter baseline or regions of noise. Noisy curves
        may be removed by specifying a maximum standard deviation for the baseline
        region separatley for the approach and retract curves using the inputs
        'max_approach_noise' and 'max_retract_noise' respectivley. The defaults for
        each are 1 nN.
    contact:
        Function to align the contact point of the approach curve with zero by
        finding the first point to cross the axis and interpolating between this
        and the previous value. Note this method may be unstable for noisy curves.
    plot_adjusted:
        Function to plot individual, processed force curves. Curve number to
        plot must be given as an input <.plot_adjusted(n)>
    plot_curves:
        Function to plot all processed force curves
    delete_curve:
        Function to manually remove erroneous or highly noisy curves. Specify the
        curve(s) to remove in the input <.delete_curve(i,j,k)>
    save_data:
        Function to save current variables in '.npz' format. The file name and path
        for saving can be given as an input. If no input is given a GUI will prompt
        the user to select a location and file name. In both case no file extension
        is needed.
    load_data:
        Function to import previously saved variables in '.npz' format. As with saving,
        the file name and path can be given as an input or left blank to use a GUI.
    calc_adhesion:
        Funciton to calculate the adhesion value for each curve by finding the minimum
        point on the retact curve. A histogram of calculated values can be plotted by
        setting the 'plot_hist' input to 'True'
    overlay_adhesion_hist:
        Funciton to plot a histogram comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    overlay_adhesion_box:
        Funciton to produce a box plot comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    overlay_adhesion_bar:
        Funciton to produce a bar chart comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    calc_modulus:
        Funciton to calculate the modulus value for each curve using the Hertz model.
        A histogram of calculated values can be plotted by setting the 'plot_hist'
        input to 'True'
    overlay_modulus_hist:
        Funciton to plot a histogram comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    overlay_modulus_box:
        Funciton to produce a box plot comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    overlay_modulus_bar:
        Funciton to produce a bar chart comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
    plot_adhesion_modulus:
        Function to produce a scatter plot of adhesion vs. modulus values for a single
        experiment

    Example
    -------
    An example script and full tutorial for using the AFM class is available from:
    https://github.com/crj341/nanoforce

    '''

    def __init__(self):

        self.approach = []
        self.retract = []
        self.z_position = []

        self.approach_raw = []
        self.retract_raw = []
        self.z_position_raw = []

        self.adhesion = []
        self.modulus = []

        self.data_offset = []
        self.data_length = []
        self.ramp_size = []
        self.samps_line = []
        self.zsens = []
        self.zscale = []
        self.def_sens = []
        self.spr_const = []

        self.file_name = []
        self.filelist = []
        self.run_name = []


    def set_run_name(self,_run_name):
        '''
        Manually set 'run_name' variable <.set_run_name(Sample name here)>
        '''
        self.run_name = str(_run_name)

    def set_def_sens(self,_def_sens):
        '''
        Manually set 'def_sens' variable <.set_def_sens(Value)>
        '''
        self.def_sens = _def_sens

    def set_spr_const(self,_spr_const):
        '''
        Manually set 'spr_const' variable <.set_spr_const(Value)>
        '''
        self.spr_const = _spr_const

    def run(self):
        '''
        Basic run script calling:
            .input_files()
            .nanoscope_params()
            .nanoscope_read()
            .baseline()
            .contact()
            .plot_curves()
            .calc_adhesion()
            .calc_modulus()
        '''

        self.input_files()
        self.nanoscope_params()
        self.nanoscope_read()
        self.baseline()
        self.contact()
        self.plot_curves()
        self.calc_adhesion()
        self.calc_modulus()

    def input_files(
        self,
        gui_on = True
    ):
        '''
        Function to find input file names. If called with no inputs a GUI
        window will appear prompting the user to select the first Nanoscope
        file (usually with '.000' extension). Alternativley the file path can
        be given as an input to the function <.input_files(\path\to\file\.000)>.
        The function will then search for all subsequent files with the same
        name and consecutive file extensions (.000, .001, .002, ..., .nnn)
        '''

        if gui_on:
            file = gui.fileopenbox()
        else:
            file = self.file_name

        file_no_ext = file[:-3]

        ext = '000'
        i = 0
        j = 0
        k = 0

        _filelist = list('')

        while os.path.isfile(file_no_ext + ext):
            
            _filelist.append(file_no_ext + ext)

            k += 1
            if k == 10:
                k = 0
                j += 1
            
            if j == 10:
                j = 0
                i += 1

            ext = str(i) + str(j) + str(k)

        self.filelist = _filelist

        

    def nanoscope_params(self):
        '''
        Function to extract relevant parameters from first Nanoscope file
        '''

        filename = self.filelist[0]
        data = open(filename,"r")

        nanoscope_chk = data.readline().strip()
        if not nanoscope_chk.startswith('\*Force file list'):
            Warning('Input file not nanoscope')
            return

        while data:

            line = data.readline().strip()

            if line.startswith('\*File list end'):
                break

            if len(line) < 1:
                continue
            
            if re.search('\Data offset:', line):
                _data_offset = [int(s) for s in line.split() if s.isdigit()]

            elif re.search('\Data length:', line):
                _data_length = [int(s) for s in line.split() if s.isdigit()]

            elif re.search('\@4:Ramp size:', line):
                _ramp_size = re.findall(r'\d+',line)[1] + '.' + re.findall(r'\d+',line)[2]

            elif re.search('\Samps/line', line):
                _samps_line = [int(s) for s in line.split() if s.isdigit()]

            elif re.search('@Sens. Zsens:', line):
                _zsens = re.findall(r'\d+',line)[0] + '.' + re.findall(r'\d+',line)[1]

            elif re.search('\@4:Z scale: V ', line):
                _zscale = re.findall(r'\d+',line)[1] + '.' + re.findall(r'\d+',line)[2]

            elif re.search('\@Sens. DeflSens:', line):
                _def_sens = re.findall(r'\d+',line)[0] + '.' + re.findall(r'\d+',line)[1]
                
            elif re.search('\Spring Constant:', line):
                _spr_const = re.findall(r'\d+',line)[0] + '.' + re.findall(r'\d+',line)[1]
        data.close()

        self.data_offset = int(_data_offset[0])
        self.data_length = int(_data_length[0])
        self.ramp_size = float(_ramp_size)
        self.samps_line = int(_samps_line[0])
        self.zsens = float(_zsens)
        self.zscale = float(_zscale)
        self.def_sens = float(_def_sens)
        self.spr_const = float(_spr_const)

    def nanoscope_read(self):
        '''
        Function to import raw force curve data
        '''

        i = 0

        self.approach_raw = np.zeros([self.samps_line, len(self.filelist)])
        self.retract_raw = np.zeros([self.samps_line, len(self.filelist)])
        self.z_position_raw = np.zeros([self.samps_line, len(self.filelist)])

        for filename in self.filelist:

            data = open(filename,"rb")

            data.seek(self.data_offset)
            _approach_curve = np.fromfile(file = data, dtype = 'int16', count = self.samps_line) * self.zscale * self.spr_const * self.def_sens
            _retract_curve = np.fromfile(file = data, dtype = 'int16', count = self.samps_line) * self.zscale * self.spr_const * self.def_sens

            hscale = (self.zsens * self.ramp_size * 2 ** 16 / self.samps_line ) / 2
            _z_position = hscale * np.arange(1, self.samps_line + 1)

            self.approach_raw[:,i] = np.array(_approach_curve)
            self.retract_raw[:,i] = np.array(_retract_curve)
            self.z_position_raw[:,i] = np.array(_z_position)

            i += 1

    def plot_raw(
        self,
        curve_number = 1
    ):
        '''
        Function to plot individual, unprocessed force curves. Curve number to
        plot must be given as an input <.plot_raw(n)>
        '''
        fig = go.Figure()
        fig.add_trace(go.Scatter(y = self.approach_raw[:,curve_number], x = self.z_position_raw[:,curve_number], mode='lines', name = 'Approach curve (raw)'))
        fig.add_trace(go.Scatter(y = self.retract_raw[:,curve_number], x = self.z_position_raw[:,curve_number], mode='lines', name = 'Retract curve (raw)'))
        fig.update_layout(
            yaxis_title_text='Force (nN)',
            xaxis_title_text='Tip separation (nm)',
        )
        fig.show()

    def baseline(
        self,
        start_pos = 0.45,
        end_pos = 0.8,
        max_approach_noise = 1,
        max_retract_noise = 1
    ):
        '''
        Function to adjust force curve baseline to zero by taking the mean value
        of a specified region of the approach curve. The default is from 0.45*z
        to 0.8*z, starting at the contact side of the curve. The start and end
        points may be specified using the 'start_pos' and 'end_pos' inputs to
        account for a longer/shorter baseline or regions of noise. Noisy curves
        may be removed by specifying a maximum standard deviation for the baseline
        region separatley for the approach and retract curves using the inputs
        'max_approach_noise' and 'max_retract_noise' respectivley. The defaults for
        each are 1 nN.
        '''
        self.approach = np.zeros([self.samps_line, len(self.filelist)])
        self.retract = np.zeros([self.samps_line, len(self.filelist)])
        self.z_position = np.zeros([self.samps_line, len(self.filelist)])

        i = 0
        
        for curve_number in range(0, self.z_position_raw.shape[1]):

            _approach = self.approach_raw[:,curve_number]
            _retract = self.retract_raw[:,curve_number]
            _z = self.z_position_raw[:,curve_number]

            baseline_start = int(len(_z) * start_pos)
            baseline_end = int(len(_z) * end_pos)
            baseline_region = np.arange(baseline_start,baseline_end)

            _approach_region = _approach[baseline_region]
            _retract_region = _retract[baseline_region]
            _z_region = _z[baseline_region]

            _approach_mean = np.mean(_approach_region)

            adjusted_approach = _approach - _approach_mean
            adjusted_retract =  _retract - _approach_mean

            approach_dev = np.std(_approach_region)
            retract_dev = np.std(_retract_region)

            if approach_dev < max_approach_noise and retract_dev < max_retract_noise:
                self.approach[:,i] = adjusted_approach
                self.retract[:,i] = adjusted_retract
                self.z_position[:,i] = _z

                i += 1

        self.approach = self.approach[:,:i]
        self.retract = self.retract[:,:i]
        self.z_position = self.z_position[:,:i]

    def contact(self):
        '''
        Function to align the contact point of the approach curve with zero by
        finding the first point to cross the axis and interpolating between this
        and the previous value. Note this method may be unstable for noisy curves.
        '''
        
        for i in range(0, self.z_position.shape[1]):
            _approach = self.approach[:,i]
            _retract = self.retract[:,i]
            _z = self.z_position[:,i]

            j = 0

            while _approach[j] > 0:
                j += 1

            y = _z[j-1], _z[j]
            x = _approach[j-1], _approach[j]

            f = interp.interp1d(x,y)
            intercept = f(0)

            self.z_position[:,i] = self.z_position[:,i] - intercept


    def plot_adjusted(
        self,
        curve_number = 1
    ):
        '''
        Function to plot individual, processed force curves. Curve number to
        plot must be given as an input <.plot_adjusted(n)>
        '''
        fig = go.Figure()
        fig.add_trace(go.Scatter(y = self.approach[:,curve_number], x = self.z_position[:,curve_number], mode='lines', name = 'Approach curve'))
        fig.add_trace(go.Scatter(y = self.retract[:,curve_number], x = self.z_position[:,curve_number], mode='lines', name = 'Retract curve'))
        fig.update_layout(
            yaxis_title_text='Force (nN)',
            xaxis_title_text='Tip separation (nm)',
        )
        fig.show()

    def plot_curves(self):
        '''
        Function to plot all processed force curves
        '''

        trace_1 = []
        trace_2 = []

        for step in range(0, self.z_position.shape[1]):

            trace_1.append(go.Scatter(
                        visible=False,
                        name= 'Approach ' + str(step),
                        y = self.approach[:, step],
                        x = self.z_position[:, step]
                    ))

            trace_2.append(go.Scatter(
                        visible=False,
                        name= 'Retract ' + str(step),
                        y = self.retract[:, step],
                        x = self.z_position[:, step]
                    ))

        fig = go.Figure(data = trace_1 + trace_2)

        n_trace = int(len(fig.data) / 2)

        # Make 1st trace visible
        fig.data[0].visible = True
        fig.data[n_trace].visible = True

        # Create and add slider
        steps = []
        for i in range(n_trace):
            step = dict(
                method="update",
                args=[{"visible": [False] * len(fig.data)},
                    {"title": "Slider switched to experiment: " + str(i)}],
            )
            j = int(n_trace+i)
            step["args"][0]["visible"][i] = True
            step["args"][0]["visible"][j] = True
            steps.append(step)

        sliders = [dict(
            pad={"t": 50},
            steps=steps
        )]

        fig.update_layout(
            sliders=sliders,
            yaxis_title_text='Force (nN)',
            xaxis_title_text='Tip separation (nm)',
        )

        fig.show()

    def delete_curve(self,*args):
        '''
        Function to manually remove erroneous or highly noisy curves. Specify the
        curve(s) to remove in the input <.delete_curve(i,j,k)>
        '''

        for curve_number in args:

            self.approach = np.delete(self.approach, curve_number, axis = 1)
            self.retract = np.delete(self.retract, curve_number, axis = 1)
            self.z_position = np.delete(self.z_position, curve_number, axis = 1)

            try:
                self.adhesion = np.delete(self.adhesion, curve_number, axis = 1)
            except:
                pass

            try:
                self.modulus = np.delete(self.modulus, curve_number, axis = 1)
            except:
                pass
    
    def save_data(
        self,
        filename = []
    ):
        '''
        Function to save current variables in '.npz' format. The file name and path
        for saving can be given as an input. If no input is given a GUI will prompt
        the user to select a location and file name. In both case no file extension
        '''
        if len(filename) == 0:
            output_path = gui.filesavebox()
        else:
            output_path = filename

        np.savez(
            output_path,
            approach = self.approach,
            retract = self.retract,
            z_position = self.z_position,
            adhesion = self.adhesion,
            modulus = self.modulus,
            data_offset = self.data_offset,
            data_length =self.data_length,
            ramp_size = self.ramp_size,
            samps_line = self.samps_line,
            zsens = self.zsens,
            zscale = self.zscale,
            def_sens = self.def_sens,
            spr_const = self.spr_const,
            file_name = self.file_name,
            run_name = self.run_name
        )

    def load_data(
        self,
        filename = []
    ):
        '''
        Function to import previously saved variables in '.npz' format. As with saving,
        the file name and path can be given as an input or left blank to use a GUI.
        '''

        if len(filename) == 0:
            input_path = gui.fileopenbox()
        else:
            input_path = filename
        
        data = np.load(input_path)

        self.approach = data['approach']
        self.retract = data['retract']
        self.z_position = data['z_position']
        self.adhesion = data['adhesion']
        self.modulus = data['modulus']
        self.data_offset = data['data_offset']
        self.data_length = data['data_length']
        self.ramp_size = data['ramp_size']
        self.samps_line = data['samps_line']
        self.zsens = data['zsens']
        self.zscale = data['zscale']
        self.def_sens = data['def_sens']
        self.spr_const = data['spr_const']
        self.file_name = data['file_name']
        self.run_name = data['run_name']

    def calc_adhesion(
        self,
        plot_hist = False
    ):
        '''
        Funciton to calculate the adhesion value for each curve by finding the minimum
        point on the retact curve. A histogram of calculated values can be plotted by
        setting the 'plot_hist' input to 'True'
        '''

        self.adhesion = np.zeros(self.z_position.shape[1])

        for curve_number in range(0, self.z_position.shape[1]):
            self.adhesion[curve_number] = np.amin(self.retract[:, curve_number]) * -1

        if plot_hist is True:
            fig = go.Figure(data = [go.Histogram(x = self.adhesion)])
            fig.update_layout(
                xaxis_title_text='Adhesion (nN)',
                yaxis_title_text='Count',
            )
            fig.show()

        print('')
        try:
            print(self.run_name)
        except:
            pass

        print('Mean adhesion = %(mean).2f nN' % {'mean': np.mean(self.adhesion)})
        print('Standard deviation = %(std).2f nN' % {'std': np.std(self.adhesion)})

    def overlay_adhesion_hist(
        self,
        *all_adhesion
    ):
        '''
        Funciton to plot a histogram comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()

        for single in all_adhesion:
            fig.add_trace(go.Histogram(x = single.adhesion,
            name = single.run_name
        ))

        fig.update_layout(
                xaxis_title_text='Adhesion (nN)',
                yaxis_title_text='Count',
                barmode='overlay'
            )
        fig.update_traces(opacity=0.75)
        fig.show()

    def overlay_adhesion_box(
        self,
        *all_adhesion
    ):
        '''
        Funciton to produce a box plot comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()

        for single in all_adhesion:
            fig.add_trace(go.Box(x = single.adhesion,
            name = single.run_name
        ))

        fig.update_layout(
                xaxis_title_text='Adhesion (nN)'
            )
            
        fig.show()

    def overlay_adhesion_bar(
        self,
        *all_adhesion
    ):
        '''
        Funciton to produce a bar chart comparing adhesion values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()

        mean_adhesion = []
        std_adhesion = []
        trace = []

        for single in all_adhesion:
            mean_adhesion = np.append(mean_adhesion, np.mean(single.adhesion))
            std_adhesion = np.append(std_adhesion, np.std(single.adhesion))
            trace = np.append(trace, single.run_name)
        
        fig.add_trace(go.Bar(
            y = mean_adhesion,
            x = trace,
            error_y = dict(type = 'data', array = std_adhesion)
        ))

        fig.update_layout(
                yaxis_title_text='Adhesion (nN)'
            )
            
        fig.show()

    def calc_modulus(
        self,
        poisson_ratio = 0.5,
        indenter_radius = 10 * 10e-9,
        plot_hist = False
    ):
        '''
        Funciton to calculate the modulus value for each curve using the Hertz model.
        A histogram of calculated values can be plotted by setting the 'plot_hist'
        input to 'True'
        '''
        self.modulus = np.zeros(self.z_position.shape[1])

        for curve_number in range(0, self.z_position.shape[1]):
            _approach = self.approach[:,curve_number]
            _z = self.z_position[:,curve_number]

            force = _approach[0] * 10e-9
            depth = _z[0] * 10e-9
            ratio = depth / indenter_radius

            term_1 = (4/3)
            term_2 = 1 / (1 - poisson_ratio ** 2)
            term_3 = 1 / force
            term_4 = np.sqrt(indenter_radius)
            term_5 = np.absolute(depth) ** (3/2)
            # term_6 = (1 - (1/10) * ratio - (1/840) * ratio ** 2 + (11/15120) * ratio ** 3 + (1357/6652800) * ratio ** 4)

            terms = term_1 * term_2 * term_3 * term_4 * term_5 #* term_6

            self.modulus[curve_number] = (1 / terms) * 10e-6
        
        if plot_hist is True:
            fig = go.Figure(data=[go.Histogram(x = self.modulus)])
            fig.update_layout(
                xaxis_title_text='Modulus (MPa)',
                yaxis_title_text='Count',
            )
            fig.show()

        print('')
        try:
            print(self.run_name)
        except:
            pass

        print('Mean modulus = %(mean).2f MPa' % {'mean': np.mean(self.modulus)})
        print('Standard deviation = %(std).2f MPa' % {'std': np.std(self.modulus)})

    def overlay_modulus_hist(
        self,
        *all_modulus
    ):
        '''
        Funciton to plot a histogram comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()

        for single in all_modulus:
            fig.add_trace(go.Histogram(x = single.modulus,
            name = single.run_name
        ))

        fig.update_layout(
                xaxis_title_text='Modulus (MPa)',
                yaxis_title_text='Count',
                barmode='overlay'
            )
        fig.update_traces(opacity=0.75)
        fig.show()
    
    def overlay_modulus_box(
        self,
        *all_modulus
    ):
        '''
        Funciton to produce a box plot comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()

        for single in all_modulus:
            fig.add_trace(go.Box(x = single.modulus,
            name = single.run_name,
        ))

        fig.update_layout(
                xaxis_title_text='Modulus (MPa)'
            )
            
        fig.show()

    def overlay_modulus_bar(
        self,
        *all_modulus
    ):
        '''
        Funciton to produce a bar chart comparing modulus values for multiple experiments.
        A separate instance of the AFM class containing the data for each should be
        given as the inputs
        '''
        fig = go.Figure()
        mean_modulus = []
        std_modulus= []
        trace = []

        for single in all_modulus:
            mean_modulus = np.append(mean_modulus, np.mean(single.modulus))
            std_modulus = np.append(std_modulus, np.std(single.modulus))
            trace = np.append(trace, single.run_name)
        
        fig.add_trace(go.Bar(
            y = mean_modulus,
            x = trace,
            error_y = dict(type = 'data', array = std_modulus)
        ))

        fig.update_layout(
                yaxis_title_text='Modulus (MPa)',
            )

        fig.show()

    def plot_adhesion_modulus(self):
        '''
        Function to produce a scatter plot of adhesion vs. modulus values for a single
        experiment
        '''
        fig = go.Figure(go.Scatter(x = self.adhesion, y = self.modulus, mode='markers'))
        fig.update_layout(
                xaxis_title_text='Adhesion (nN)',
                yaxis_title_text='Modulus (MPa)'
            )
        fig.show()