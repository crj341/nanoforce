#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File      : afm_analysis.py
# Author    : Chris Jones
# Email     : crj341@student.bham.ac.uk
# Date      : 26/06/2020

# Import AFM class
from  nanoforce   import  AFM

# Name varialbles for each experiment
expt_1 = AFM()
expt_2 = AFM()
expt_3 = AFM()

# Name each experiment
expt_1.set_run_name('Sample 1')
expt_2.set_run_name('Sample 2')
expt_3.set_run_name('Sample 3')

# Select files to import
expt_1.input_files()
expt_2.input_files()
expt_3.input_files()

# Extract parameters from input file
expt_1.nanoscope_params()
expt_2.nanoscope_params()
expt_3.nanoscope_params()

# Manually set deflection sensetivity
expt_1.set_def_sens(100)
expt_2.set_def_sens(100)
expt_3.set_def_sens(100)

# Manually set spring constant
expt_1.set_spr_const(0.32)
expt_2.set_spr_const(0.32)
expt_3.set_spr_const(0.32)

# Read raw data
expt_1.nanoscope_read()
expt_2.nanoscope_read()
expt_3.nanoscope_read()

# Adjust the baseline
expt_1.baseline()
expt_2.baseline()
expt_3.baseline()

# Set contact point
expt_1.contact()
expt_2.contact()
expt_3.contact()

# Plot all curves
expt_1.plot_curves()
expt_2.plot_curves()
expt_3.plot_curves()

# Calculate adhesion
expt_1.calc_adhesion()
expt_2.calc_adhesion()
expt_3.calc_adhesion()

# Calculate modulus
expt_1.calc_modulus()
expt_2.calc_modulus()
expt_3.calc_modulus()

# Save files
expt_1.save_data('sample_1')
expt_2.save_data('sample_2')
expt_3.save_data('sample_3')

# Plot histogram to compare adhesion
AFM.overlay_adhesion_hist(expt_1,expt_2,expt_3)

# PLot bar chart to compare modulus
AFM.overlay_modulus_bar(expt_1,expt_2,expt_3)