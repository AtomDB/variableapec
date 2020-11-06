"""This module contains methods for looking at emission and line ratios
as a function of varying atomic data from the AtomDB files. Requires
PyAtomdB, Python 3, LaTeX, (and astroML if you wish to use LaTeX font with plots)."""

__version__ = "1.0.2" #Nov 3, 2020
__author__ = "Keri Heuer"
__email__ = "kh3286@drexel.edu"
__status__ = "Development"

import matplotlib.pyplot as plt, matplotlib.ticker as mtick, scipy.stats as stats, \
pyatomdb, numpy, pickle, pathlib, csv, os, errno, hashlib, requests, urllib.request, \
urllib.parse, urllib.error, time, subprocess, shutil, wget, glob, datetime, ftplib, \
pathlib, collections, operator, requests, matplotlib.pylab as pylab, glob, math
from io import StringIO
from matplotlib.ticker import FormatStrFormatter, LogFormatterSciNotation
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.table import Table, Column, vstack
from matplotlib.offsetbox import AnchoredText
from decimal import Decimal
from matplotlib import cm
try:
  import astropy.io.fits as pyfits
except:
  import pyfits

""" User can/should update plot settings before calling
any plotting routine. Below are the defaults.
cmap : str
    matplotlib color map string name
fontsize : int
    size of x and y axis labels
ticksize : int
    size of x and y ticks
legendsize : str
    string of legend size 
legendloc : int or str 
    legend location (default is 'best')
tickstep : float
    size of y ticks for error (bottom) panels of plots
important_ticks : tuple or list
    list of important y ticks for error (bottom) panels of plots, 
    i.e. 0.1 to draw a tick at 0.1 fractional error if not already in default ticks
grid_kwargs : dict
    dictionary of optional arguments for turning on axis grid
    for error (bottom) panels, i.e. keys include 'color', 'alpha',
    'linestyle', 'which', 'axis', etc.
use_latex_font : bool 
    True to use LaTeX font on plots (default)
    will be False if cannot find astroML or LaTeX installed"""

global plot_settings
plot_settings = {'cmap': 'hsv', 'fontsize': 12, 'ticksize': 12,
                 'legendsize': 'xx-small', 'legendloc': 'best', 'tickstep': 0.1,
                 'important_ticks': {}, 'grid_kwargs': {}, 'use_latex_font': False}
if plot_settings['use_latex_font'] == True:
    try:
        from astroML.plotting import setup_text_plots
    except:
        use_latex_font = False
        print("Using regular matplotlib font for plotting. Install astroML to setup plots with LaTeX font.")
#if setup_text_plots in globals():
    #setup_text_plots(uselatex = plot_settings['use_latex_font'])

#set up datacache anytime variableapec is imported
global d
d = {}

""" CSD errors: to be added to
These are typical systematic uncertainties from experiments
measuring ionization and recombination rate coefficients. 
To add to this, {key: val} should be
{Z: {'ionize': {z1: fractional error}, 'recomb':{z1: fractional error}}}

*note: the max_ionize and max_recomb dictionaries do not need to have
errors for every z1, variableapec will fill in the dictionary using 
the average of adjacent ions. Adding even a single z1 error will be
useful for users to calculate CSD error estimates."""

systematics = {8: {'ionize': {8: 0.1, 7: 0.2, 6: 0.05, 5: 0.05, 4: 0.05, 3: 0.05, 2: 0.05, 1: 0.05},
                   'recomb': {8: 0.18, 7: 0.2, 6: 0.16, 5: 0.3, 4: 0.3, 3: 0.3, 2: 0.3, 1: 0.3}},
               26: {'ionize': {26: 0.15, 25: 0.14, 24: 0.13, 23: 0.265, 22: 0.4, 21: 0.4, 20: 0.4,
                                   19: 0.4, 18: 0.15, 17: 0.16, 16: 0.2, 15: 0.26, 14: 0.06, 13: 0.06,
                                   12: 0.14, 11: 0.09, 10: 0.09, 9: 0.16, 8: 0.12, 7: 0.12, 6: 0.12,
                                   5: 0.12, 4: 0.13, 3: 0.12, 2: 0.12, 1: 0.12},
                    'recomb':{26: 0.2, 25: 0.2, 24: 0.2, 23: 0.2, 22: 0.2, 21: 0.2, 20: 0.2,
                                  19: 0.2, 18: 0.2, 17: 0.2, 16: 0.2, 15: 0.26, 14: 0.18, 13: 0.18,
                                  12: 0.13, 11: 0.25, 10: 0.25, 9: 0.29, 8: 0.25, 7: 0.25, 6: 0.25,
                                  5: 0.25, 4: 0.25, 3: 0.25, 2: 0.25, 1: 0.25}}}

def ionize(Z, z1, Te, dens, in_range, pop_fraction, datacache):
    z1_drv = z1-1

    #get rates and set up CR matrix
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1_drv, Te, dens, do_la= True, \
                            do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
    
    lvdat = pyatomdb.atomdb.get_data(Z, z1_drv, 'LV', datacache=datacache)
    lvdat = lvdat[1].data
    nlev = len(lvdat)
    drv_matrix = numpy.zeros((nlev,nlev))
    drv_B = numpy.zeros(nlev)
    
    #populate full CR matrix by summing rates for all processes
    for x,y,z in zip(final, init, rates):
        drv_matrix[x][y] += z
        
    #set up and solve CR matrix for level populations
    drv_matrix[0][:], drv_B[0] = 1.0, 1.0
    drv_lev_pop = numpy.linalg.solve(drv_matrix,drv_B)

    #get level populations and make linelist
    ion_levpop = pyatomdb.apec.calc_ioniz_popn(drv_lev_pop*pop_fraction[z1_drv-1], Z, z1, z1_drv, Te, dens, datacache=datacache, do_xi=True)
    ion_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    ion_linelist['upperlev'], ion_linelist['lowerlev'] = in_range['UPPER_LEV'], in_range['LOWER_LEV']

    #use observational wavelength if exists, otherwise theoretical
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == True)
    ion_linelist[i]['lambda'] = in_range[i]['WAVELEN']
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == False)
    ion_linelist[i]['lambda'] = in_range[i]['WAVE_OBS']

    #calculate epsilon
    ion_linelist['epsilon'] = in_range['EINSTEIN_A'] * ion_levpop[in_range['UPPER_LEV'] - 1]
    return ion_linelist
            
def recombine(Z, z1, Te, dens, in_range, pop_fraction, datacache):
    z1_drv = z1+1
    if z1 < Z:
        init, final, rates = pyatomdb.apec.gather_rates(Z, z1_drv, Te, dens, do_la=True, \
                                                        do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)

        lvdat = pyatomdb.atomdb.get_data(Z, z1_drv, 'LV', datacache=datacache)
        lvdat = lvdat[1].data
        nlev = len(lvdat)
        drv_matrix = numpy.zeros((nlev, nlev))
        drv_B = numpy.zeros(nlev)

        # populate full CR matrix by summing rates for all processes
        for x, y, z in zip(final, init, rates):
            drv_matrix[x][y] += z

        # set up and solve CR matrix for level populations
        drv_matrix[0][:], drv_B[0] = 1.0, 1.0
        drv_lev_pop = numpy.linalg.solve(drv_matrix, drv_B)

    else:
        # declare 1 level, fully populated
        drv_lev_pop = numpy.ones(1, dtype=float)

    # get level populations and make linelist
    recomb_levpop = pyatomdb.apec.calc_recomb_popn(drv_lev_pop * pop_fraction[z1_drv - 1], Z, z1, z1_drv, Te, dens,
                                                       drlevrates=0, rrlevrates=0,datacache=datacache)
    recomb_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    recomb_linelist['upperlev'], recomb_linelist['lowerlev'] = in_range['UPPER_LEV'], in_range['LOWER_LEV']

    #use observational wavelength if exists, otherwise theoretical
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == True)
    recomb_linelist[i]['lambda'] = in_range[i]['WAVELEN']
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == False)
    recomb_linelist[i]['lambda'] = in_range[i]['WAVE_OBS']

    #calculate epsilon
    recomb_linelist['epsilon'] = in_range['EINSTEIN_A'] * recomb_levpop[in_range['UPPER_LEV'] - 1]
    return recomb_linelist

def set_up(Z, z1, Te, dens, extras={}):
    """
    Z : int
        element
    z1 : int
        ion charge +1
    Te : int
        temperature in K
    dens :
        density in cm-3
    extras : dict
        dictionary of extra inputs if varying rate,
        dictionary is empty if running set_up() on its own
        (i.e. to get original epsilons at Te, dens)
        should be dict of {}

    Sets up linelist from excitation, recombination and ionization
    and creates sensitivity tables. Returns a dictionary "inputs"
    and the transition varied if extras is not empty."""

    print("\n***********************************************************************************")
    print("Setting up linelists at Te={:.3e} and dens={:.3e}".format(Te, dens))

    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)

    lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV')
    lvdat = lvdat[1].data
    nlev = len(lvdat)

    matrix = numpy.zeros((nlev, nlev))
    B = numpy.zeros(nlev)
    # populate full CR matrix by summing rates for all processes
    for x,y,z in zip(final, init, rates):
        matrix[x][y] += z

    # set up and solve CR matrix for level populations
    matrix[0][:] = 1.0
    B[0] = 1.0
    lev_pop = numpy.linalg.solve(matrix, B)

    # convert level populations into line lists & intensities for excitation only
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
    in_range = ladat[1].data

    # find fraction of each ion in plasma
    pop_frac = pyatomdb.apec._solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

    #get epsilon from excitation only
    linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    linelist['epsilon'] = in_range['EINSTEIN_A'] * lev_pop[in_range['UPPER_LEV'] - 1] * pop_frac[z1-1]
    linelist['upperlev'], linelist['lowerlev'] = in_range['UPPER_LEV'], in_range['LOWER_LEV']

    # get wavelengths - use observational lambda if exists, otherwise theoretical
    wavelens = numpy.zeros(len(in_range))
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == True)[0]
    wavelens[i] = in_range[i]['WAVELEN']
    i = numpy.where(numpy.isnan(in_range['WAVE_OBS']) == False)[0]
    wavelens[i] = in_range[i]['WAVE_OBS']
    linelist['lambda'] = wavelens

    # set up complete line list (emiss only due to excitation at this point)
    full_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    full_linelist['lambda'] = linelist['lambda']
    full_linelist['epsilon'] = linelist['epsilon']

    # now add emissivity from ionization and recombination to excitation linelist (depending on z1)
    if z1 == 1:  # skip ionization
        recomb_emiss = recombine(Z, z1, Te, dens,in_range, pop_frac, d)
        full_linelist['epsilon'] += recomb_emiss['epsilon']
    elif z1 == Z + 1:  # skip recombination
        ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_frac, d)
        full_linelist['epsilon'] += ion_emiss['epsilon']
    else:  # do both
        recomb_emiss = recombine(Z, z1, Te, dens,in_range, pop_frac, d)
        ion_emiss = ionize(Z, z1, Te, dens,in_range, pop_frac, d)
        full_linelist['epsilon'] += recomb_emiss['epsilon']
        full_linelist['epsilon'] += ion_emiss['epsilon']

    # set up sensitivity & partial derivatives tables
    table = Table([full_linelist['lambda'], in_range['UPPER_LEV'].astype(int), in_range['LOWER_LEV'].astype(int), full_linelist['epsilon']], \
                  names=('Lambda', 'Upper', 'Lower', 'Epsilon_orig'))
    new_table = Table([full_linelist['lambda'], in_range['UPPER_LEV'].astype(int), in_range['LOWER_LEV'].astype(int), full_linelist['epsilon']],
                      names=('Lambda', 'Upper', 'Lower', 'Epsilon_orig'))

    # save variables
    inputs = {'Z': Z, 'z1': z1, 'Te': Te, 'dens': dens}
    values = {'matrix': matrix, 'B': B, 'in_range': in_range, 'linelist': linelist, 'table': table,
              'new_table': new_table}

    if extras != {}:
        process, delta_r, transition, transition_2, wavelen, \
        Te_range, dens_range, corrthresh, e_signif, pop_fraction = [extras.get(k) for k in extras]

        #print out emissivity from excitation, ionization, recombination
        if process == 'exc':
            i = numpy.where([(linelist['upperlev'] == transition[1]) & (linelist['lowerlev'] == transition[0])])[1]
            #check for multiple wavelengths with same (up, lo), prompt for theoretical wavelen if so
            if len(i) > 1:
                wave = input("Found multiple lines w/ (up,lo) = " + str(
                    transition) + ", input theoretical wavelength in A of your transition to continue: ")
                i = numpy.where([math.isclose(in_range['WAVELEN'], float(wave), rel_tol=1e-5) == True])[1][0]
                transition = (in_range[i]['UPPER_LEV'], in_range[i]['LOWER_LEV'])
            else:
                i = i[0]
            
            print("Excitation emiss =", linelist[i]['epsilon'])
            i = numpy.where([(ion_emiss['upperlev'] == transition[1]) & (ion_emiss['lowerlev'] == transition[0])])[1][0]
            print("Ionization emiss =", ion_emiss[i]['epsilon'])
            i = numpy.where([(recomb_emiss['upperlev'] == transition[1]) & (recomb_emiss['lowerlev'] == transition[0])])[1][0]
            print("Recombination emiss =", recomb_emiss[i]['epsilon'])
        elif process == 'A':
            i = numpy.where([(linelist['upperlev'] == transition[0]) & (linelist['lowerlev'] == transition[1])])[1]
            print(i)
            #check for multiple wavelengths with same (up, lo), prompt for theoretical wavelen if so
            if len(i) > 1:
                wave = input("Found multiple lines w/ (up,lo) = " + str(
                    transition) + ", input theoretical wavelength in A of your transition to continue: ")
                i = numpy.where([math.isclose(in_range['WAVELEN'], float(wave), rel_tol=1e-5) == True])
                transition = (in_range[i]['UPPER_LEV'], in_range[i]['LOWER_LEV'])
            else:
                i = i[0]
            print("Excitation emiss =", linelist[i]['epsilon'])
            i = numpy.where([(ion_emiss['upperlev'] == transition[0]) & (ion_emiss['lowerlev'] == transition[1])])[1][0]
            print("Ionization emiss =", ion_emiss[i]['epsilon'])
            i = numpy.where([(recomb_emiss['upperlev'] == transition[0]) & (recomb_emiss['lowerlev'] == transition[1])])[1][0]
            print("Recombination emiss =", recomb_emiss[i]['epsilon'])

        inputs.update({'process': process, 'delta_r': delta_r, 'transition': transition, 'transition_2': transition_2,
                       'wavelen': wavelen, 'Te_range': Te_range, 'dens_range': dens_range,
                       'corrthresh': corrthresh, 'e_signif': e_signif, 'pop_fraction': pop_frac})
        #use pop_fraction specified by user if not {}
        if isinstance(pop_fraction, numpy.ndarray):
            inputs['pop_fraction'] = pop_fraction
        values.update({'index': i})
        return inputs, values, transition
    else:
        return inputs, values

def vary_a(inputs, values, which_transition):
    """
    inputs : dict
        returned by set_up()
    values : dict
        returned by set_up()
    which_transition : tuple
        returned by set_up()

    Recalculates the emissivities from varying
    the A value of which_transition and updates
    sensitivity tables and values dictionary.
    Returns the dictionaries input and values."""

    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    wavelen, Te_range, dens_range, corrthresh, e_signif, pop_fraction = [inputs.get(k) for k in inputs]

    matrix, B, in_range, linelist, table, new_table, idx = [values.get(k) for k in values]

    print("Varying A-value for", str(which_transition), "by", str(delta_r*100)+"%")

    initial_lev, final_lev = which_transition[0], which_transition[1]

    a_index = numpy.where([(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)])[1][0]
    old_a = in_range[a_index]['EINSTEIN_A']
    table['Epsilon_orig'].unit = 'A'

    if old_a == 0: old_a = 1e-40

    # vary A values
    new_a = [old_a * (1 - delta_r), old_a * (1 + delta_r)]
    q_min, q_max = new_a[0], new_a[1]

    # store the original ARAD_TOT, and invert
    lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV', datacache=d)
    orig_arad_tot = lvdat[1].data['ARAD_TOT'][initial_lev - 1] * 1.0

    excitation, ionization, recombination = [], [], []
    for index, x in enumerate(new_a):

        # update LA data for specified transition
        in_range['EINSTEIN_A'][a_index] = x
        lvdat[1].data['ARAD_TOT'][initial_lev - 1] = orig_arad_tot + x - old_a

        # get new CR matrix and resolve level pops with new A
        frac = str(round(x / old_a, 2))
        new_matrix = matrix.copy()
        new_matrix[final_lev - 1, initial_lev - 1] += (x - old_a)  # off diagonal term
        new_matrix[initial_lev - 1, initial_lev - 1] -= (x - old_a)  # diagonal term
        new_matrix[0][:] = 1.0

        # find new level populations and get new epsilon values from excitation
        new_lev_pop = numpy.linalg.solve(new_matrix, B) * pop_fraction[z1-1]

        #fill in new_linelist
        new_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
        new_linelist['epsilon'] = in_range['EINSTEIN_A'] * new_lev_pop[in_range['UPPER_LEV'] - 1]
        new_linelist['upperlev'], new_linelist['lowerlev'] = in_range['UPPER_LEV'], in_range['LOWER_LEV']

        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += ion_emiss['epsilon']
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, d)
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
            new_linelist['epsilon'] += ion_emiss['epsilon']

        #update sensitivity table
        if index == 0:
            print("Adding the min epsilon", new_linelist[idx]['epsilon'])
            excitation.append(new_linelist[idx]['epsilon'])
            if 'ion_emiss' in locals(): ionization.append(ion_emiss[idx]['epsilon'])
            if 'recomb_emiss' in locals(): recombination.append(recomb_emiss[idx]['epsilon'])
            new_col = Column(name='Epsilon_min', data=new_linelist['epsilon'], unit=frac + ' A')
        elif index == 1:
            print("Adding the max epsilon", new_linelist[idx]['epsilon'])
            excitation.append(new_linelist[idx]['epsilon'])
            if 'ion_emiss' in locals(): ionization.append(ion_emiss[idx]['epsilon'])
            if 'recomb_emiss' in locals(): recombination.append(recomb_emiss[idx]['epsilon'])
            new_col = Column(name='Epsilon_max', data=new_linelist['epsilon'], unit=frac + ' A')

        table.add_columns([new_col])

    #reset LA data to old_a
    in_range['EINSTEIN_A'][a_index] = old_a
    lvdat[1].data['ARAD_TOT'][initial_lev - 1] = orig_arad_tot

    for label, array in zip(['excitation', 'ionization', 'recombination'], [excitation, ionization, recombination]):
        print("Min and max", str(label), "emiss =", array)
    print("***********************************************************************************")

    values = {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, 'q_max': q_max, 'q_min': q_min}
    return inputs, values

def vary_exc(inputs, values, which_transition):
    """
    inputs : dict
        returned by set_up()
    values : dict
        returned by set_up()
    which_transition : tuple
        returned by set_up()

    Recalculates the emissivities from varying
    the direct excitation rate of which_transition
    and updates sensitivity tables and values dictionary.
    Returns the dictionaries input and values."""

    Z, z1, Te, dens, process, delta_r, transition, transition_2, wavelen, \
    Te_range, dens_range, corrthresh, e_signif, pop_fraction = [inputs.get(k) for k in inputs]

    matrix, B, in_range, linelist, table, new_table, idx = [values.get(k) for k in values]

    print("Varying excitation rate for", str(which_transition), "by", str(delta_r*100)+"%")

    exc_init, exc_final, exc_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=False, \
                                                                do_ec=True, do_ir=False, do_pc=False, do_ai=False,
                                                                datacache=d)

    initial_lev, final_lev = which_transition[0], which_transition[1]

    #get original excitation rate of which_transition
    try:
        i = numpy.where([(exc_init == which_transition[0]-1) & (exc_final == which_transition[1]-1)])[1][0]
        old_rate = exc_rates[i]
        if old_rate == 0: old_rate = 1e-40
    except:
        print("Could not find transition", which_transition, " - please check input transition levels")

    table['Epsilon_orig'].unit = 'orig rate'

    # vary rate
    new_rate = [(1-delta_r) * old_rate, (1+delta_r) * old_rate]
    q_max, q_min = new_rate[-1], new_rate[0]

    excitation, ionization, recombination = [], [], []
    for index, x in enumerate(new_rate):
        # loop through varied rates, get new matrix and resolve level pops
        frac = str(round(x / old_rate, 2))
        new_matrix = matrix.copy()
        new_matrix[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        new_matrix[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term
        new_matrix[0][:] = 1.0

        # find new level populations and get new epsilon values from excitation
        new_lev_pop = numpy.linalg.solve(new_matrix, B) * pop_fraction[z1 - 1]

        new_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
        new_linelist['epsilon'] = in_range['EINSTEIN_A'] * new_lev_pop[in_range['UPPER_LEV'] - 1]
        new_linelist['upperlev'], new_linelist['lowerlev'] = in_range['UPPER_LEV'], in_range['LOWER_LEV']

        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += ion_emiss['epsilon']
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, d)
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, d)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
            new_linelist['epsilon'] += ion_emiss['epsilon']

        # update sensitivity table and print out excitation, ionization, and recombination emiss
        if index == 0:
            excitation.append(new_linelist[idx]['epsilon'])
            if 'ion_emiss' in locals(): ionization.append(ion_emiss[idx]['epsilon'])
            if 'recomb_emiss' in locals(): recombination.append(recomb_emiss[idx]['epsilon'])
            new_col = Column(name='Epsilon_min', data=new_linelist['epsilon'], unit=frac + ' rate')
        elif index == 1:
            excitation.append(new_linelist[idx]['epsilon'])
            if 'ion_emiss' in locals(): ionization.append(ion_emiss[idx]['epsilon'])
            if 'recomb_emiss' in locals(): recombination.append(recomb_emiss[idx]['epsilon'])
            new_col = Column(name='Epsilon_max', data=new_linelist['epsilon'], unit=frac + ' rate')

        table.add_columns([new_col])

    for label, array in zip(['excitation', 'ionization', 'recombination'], [excitation, ionization, recombination]):
        print("Min and max", label, "emiss =", array)
    print("***********************************************************************************")
    values = {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, 'q_max': q_max, 'q_min': q_min}
    return inputs, values

def get_tables(inputs, values):
    """
    inputs: dict outputted by vary_exc() or vary_a()
    values: dict outputted by vary_exc() or vary_a()

    Updates and returns the following:
    table (epsilon values and dE/dR)
    new_table (epsilon changes dE/E, sorted by greatest to smallest)
    inputs (dict)
    results (dict)

    Tables will be filtered if user specifies corrthresh, e_signif,
    i.e. table will only have values of dE/E > corrthresh
    and lines with original epsilon > e_signif

    Note:
    dE/dR = difference in epsilon calculated from
        rate-error and rate+error divided by the change in the rate.
    dE/E = averaged fractional change in emissivity
        (i.e. half the change in epsilon from rate-error and rate+error
        divided by the original emissivity)."""

    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    wavelen, Te_range, dens_range, corrthresh, e_signif, pop_fraction = [inputs.get(k) for k in inputs]
    table, new_table, new_linelist, q_max, q_min = [values.get(k) for k in values]
    min, max = table['Epsilon_min'], table['Epsilon_max']

    #add partial derivative dE/dR (change in emissivity dE normalized by varied rate dR) and epsilon orig
    new_table['dE/dR'] = abs((max-min)/(q_max-q_min)).round(5)
    new_table['dE/dR'][new_table['dE/dR'] < 0.0001] = 0
    new_table['dE/dR'].unit = None
    
    #add "correlation factor" dE/E -> **redefined dE/E to be the +/- error from orig epsilon**
    new_table['|dE/E|'] = ((numpy.abs(max-min)/table['Epsilon_orig'])/2).round(5)
    new_table['|dE/E|'][new_table['|dE/E|'] < 0.0001] = 0
    new_table['|dE/E|'].unit = None
    try:
      new_table.sort('|dE/E|', reverse=True)
    except TypeError:
      new_table.sort('|dE/E|')
      new_table = new_table[::-1]
    
    #apply filters
    if corrthresh != 0.0:     #only show lines whose "epsilon correlation" >= than specified value
        new_table = new_table[new_table['|dE/E|'] >= corrthresh]
    elif e_signif != 0.0:    #only show lines with partial epsilon/partial rate derivative is >= specified value
        new_table = new_table[new_table['Epsilon_orig'] >= e_signif]

    results = {'inputs': inputs, 'wavelength': table['Lambda'], 'upper': new_table['Upper'],
               'lower': new_table['Lower'], 'dE/dR': new_table['dE/dR'], 'epsilon_orig': new_table['Epsilon_orig'],\
                '|dE/E|': new_table['|dE/E|'], 'min_eps': min,'max_eps': max}
    return table, new_table, inputs, results

def is_integer(value: str, *, base: int=10) -> bool:
    try:
        int(value, base=base)
        return True
    except ValueError:
        return False

def get_file_num(filename, ext):
    a = os.listdir('./')
    b = [f.strip(ext) for f in a if (ext and filename) in f]
    numbers = [int(x[-1]) for x in b if is_integer(x[-1])]

    if len(numbers) == 0:
        return 1
    else:
        return max(numbers)+1

def plot_ratio(fnames, ratio={}, opacity=0.4, show=True, labels=[]):
    """ Plots line ratio and error from specified files.
    fnames : str or list
        contains fits file names of runs to plot
    ratio : str
        {} = regular 2 line ratio (default)
        'r' = R ratio
        'g' = G ratio
        'b' = 3 line blended ratio
    opacity : float or list
        float less than 1 supplied as alpha arg to plt.plot()
        if float, uses same opacity for all fnames
        supply list of same length of fnames for individual alphas
    show : boolean
        Show to screen if True
        (plot will be saved as pdf either way)
    labels : list
        list of legend labels in same order as files in fname
    important_ticks : list
        list of floats specifying any important ticks
        for fractional error, i.e. 0.1 for a tick to
        be drawn at 10% error
    """
    print("Plotting files:", fnames)

    #check inputs for colormap, alpha, and how many files to plot
    if isinstance(fnames, (list, tuple)):
        clist = get_cmap(len(fnames)+1, name=plot_settings['cmap'])
        if isinstance(opacity, float):
            alphas = [opacity] * len(fnames)
        else: alphas = opacity
    else:
        fnames = [fnames]
        if isinstance(opacity, float): alphas = [opacity]
        clist = get_cmap(1, name=plot_settings['cmap'])
    if labels == []:
        labels = ['']*len(fnames)

    name = fnames[0].split('_')[0]

    gs_kw = dict(width_ratios=[3], height_ratios=[2, 1])
    fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True, gridspec_kw=gs_kw)
    ax2.set_ylabel('Fractional Error', fontsize=plot_settings['labelsize'])
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for axis in [ax, ax2]:
        for which in ['x', 'y']: axis.tick_params(axis=which, labelsize=plot_settings['labelsize'])

    if ratio == {}:  # plot regular 2 line ratio
        print("Plotting", name, "line ratio")
        ax.set_ylabel('Line Ratio', fontsize=plot_settings['fontsize'])
        max_y_tick = 0
        for i, file in enumerate(fnames):
            with fits.open(file) as hdul:
                data = hdul[1].data
                if 'temp' in file:
                    type = 'temp'
                    temps, Te_orig, Te_min, Te_max = data['temps'], data['Te_orig'], data['Te_min'], data['Te_max']
                    if 'density' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                        string = "{:.0e}".format(data['density'][0]).split("e+")
                        if string[1] == '00': density = "$N_e = " + "{}".format(string[0]) + " cm^{-3}$"
                        else:
                            density = "$N_e = " + "{}\\times 10^{}".format(string[0], round(int(string[1]))) + " cm^{-3}$"

                elif 'dens' in file:
                    ax2.set_xlabel('Density in cm$^{-3}$', fontsize=plot_settings['fontsize'])
                    type = 'dens'
                    dens, dens_orig, dens_min, dens_max = data['dens'], data['dens_orig'], data['dens_min'], data['dens_max']
                    if 'temperatures' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                        string = "{:.0e}".format(data['temperatures'][0]).split("e+")
                        temperature = "$T_e$ = " + "{}\\times 10^{}$K".format(string[0], round(int(string[1])))

            if type == 'temp':  # plot emissivity versus temperature
                ax.semilogx(temps, Te_orig, color=clist(i))
                if labels[i] == '':
                    ax.fill_between(temps, Te_min, Te_max, alpha=alphas[i], color=clist(i), label=density)
                else:
                    ax.fill_between(temps, Te_min, Te_max, alpha=alphas[i], color=clist(i), label=labels[i])
                error = (abs(Te_max - Te_min) / Te_orig)/2
                if max(error) > max_y_tick: max_y_tick = max(error)
                ax2.semilogx(temps, error, color=clist(i))
                ax2.set_xlabel('Temperature in K', fontsize=plot_settings['fontsize'])
                if Te_max[-1] > Te_max[0]: left = True
                else: left = False

            elif type == 'dens':
                ax.semilogx(dens, dens_orig, color=clist(i))
                if labels[i] == '':
                    ax.fill_between(dens, dens_min, dens_max, alpha=alphas[i], color=clist(i), label=temperature)
                else:
                    ax.fill_between(dens, dens_min, dens_max, alpha=alphas[i], color=clist(i), label=labels[i])
                error = (abs(dens_max - dens_min) / dens_orig)/2
                if max(error) > max_y_tick: max_y_tick = max(error)
                ax2.semilogx(dens, error, color=clist(i))
                if dens_max[-1] > dens_max[0]: left = True
                else: left = False

        # set y ticks for error subplot if user didn't specify
        if tickstep == {}:
            if max_y_tick < 0.05: tickstep = 0.01
            elif 0.05 < max_y_tick < 0.1: tickstep = 0.02
            elif 0.1 < max_y_tick < 0.2: tickstep = 0.05
            elif 0.2 < max_y_tick < 0.4: tickstep = 0.1
            elif 0.4 < max_y_tick < 0.8: tickstep = 0.2
            elif 0.8 < max_y_tick < 1.5: tickstep = 0.25
            elif 1.5 < max_y_tick: tickstep = 0.5
        ticks = numpy.arange(0, max_y_tick + tickstep, tickstep)
        if plot_settings['important_ticks'] != {}:
            ticks += numpy.array(plot_settings['important_ticks'])
        ax2.set_yticks(ticks)
        
        if plot_settings['grid_kwargs'] != {}:
            ax2.grid(plot_settings['grid_kwargs'])

        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.95)
        if plot_settings['legendloc'] != 'best': ax.legend(fontsize=plot_settings['legendsize'], loc=plot_settings['legendloc'])
        else:
            if left == True: ax.legend(fontsize=plot_settings['legendsize'], loc='upper left')
            else: ax.legend(fontsize=plot_settings['legendsize'], loc='upper right')
        fig.align_ylabels()
        plt.savefig(name + " line ratio.pdf")

    elif ratio.lower() == 'r': # density dependent R ratio
        print("Plotting", name, "R ratio")
        ax.set_ylabel(name + ' R Ratio', fontsize=plot_settings['fontsize'])
        ax2.set_xlabel('Density in cm$^{-3}$', fontsize=plot_settings['fontsize'])

        # retrieve diagnostics
        max_y_tick = 0
        for i, file in enumerate(fnames):
            with fits.open(file) as hdul:
                data = hdul[1].data
                dens_bins, r_orig, r_min, r_max = data['dens'], data['r orig'], data['r min'], data['r max']
                f_orig, f_min, f_max = data['f orig'], data['f min'], data['f max']
                i_orig, i_min, i_max = data['i orig'], data['i min'], data['i max']
                i2_orig, i2_min, i2_max = data['i2 orig'], data['i2 min'], data['i2 max']
                if 'temperature' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                    string = "{:.0e}".format(data['temperature'][0]).split("e+")
                    temperature = "$T_e = " + "{}\\times 10^{}$K".format(string[0], round(int(string[1])))
            # do math
            R_orig = f_orig/(i_orig+i2_orig)
            numpy.nan_to_num(R_orig)
            # error propagation for positive dg
            dr, df, di, di2 = (r_max - r_orig), (f_max - f_orig), (i_max - i_orig), (i2_max - i2_orig)
            f_term = df ** 2 / (i_orig + i2_orig) ** 2
            i_term = di ** 2 * (-f_orig / (i_orig + i2_orig) ** 2) ** 2
            i2_term = di2 ** 2 * (-f_orig / (i_orig + i2_orig) ** 2) ** 2
            dR = numpy.sqrt(f_term + i_term + i2_term)
            R_max = R_orig + dR

            #error propagation for negative dg
            dr, df, di, di2 = (r_orig - r_min), (f_orig - f_min), (i_orig - i_min), (i2_orig - i2_min)
            f_term = df ** 2 / (i_orig + i2_orig) ** 2
            i_term = di ** 2 * (-f_orig / (i_orig + i2_orig) ** 2) ** 2
            i2_term = di2 ** 2 * (-f_orig / (i_orig + i2_orig) ** 2) ** 2
            dR = numpy.sqrt(f_term + i_term + i2_term)
            R_min = R_orig - dR

            if labels[i] != '':
                ax.semilogx(dens_bins, R_orig, color=clist(i), label=labels[i])
            else:
                ax.semilogx(dens_bins, R_orig, color=clist(i), label=temperature)
            ax.fill_between(dens_bins, R_min, R_max, alpha=alphas[i], color=clist(i))
            error = (abs(R_max-R_min)/R_orig)/2
            if max(error) > max_y_tick: max_y_tick = max(error)
            ax2.semilogx(dens_bins, error, color=clist(i))
            if R_orig[-1] > R_orig[0]: left = True
            else: left = False

        # set y ticks for error subplot if user didn't specify
        if tickstep == {}:
            if max_y_tick < 0.05: tickstep = 0.01
            elif 0.05 < max_y_tick < 0.1: tickstep = 0.02
            elif 0.1 < max_y_tick < 0.2: tickstep = 0.05
            elif 0.2 < max_y_tick < 0.4: tickstep = 0.1
            elif 0.4 < max_y_tick < 0.8: tickstep = 0.2
            elif 0.8 < max_y_tick < 1.5: tickstep = 0.25
            elif 1.5 < max_y_tick: tickstep = 0.5
        ticks = numpy.arange(0, max_y_tick + tickstep, tickstep)
        if plot_settings['important_ticks'] != {}:
            ticks += numpy.array(plot_settings['important_ticks'])
        ax2.set_yticks(ticks)

        if plot_settings['grid_kwargs'] != {}:
            ax2.grid(plot_settings['grid_kwargs'])

        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.95)
        if plot_settings['legendloc'] != 'best':
            ax.legend(fontsize=plot_settings['legendsize'], loc=plot_settings['legendloc'])
        else:
            if left == True:
                ax.legend(fontsize=plot_settings['legendsize'], loc='upper left')
            else:
                ax.legend(fontsize=plot_settings['legendsize'], loc='upper right')
        fig.align_ylabels()
        fig.savefig(name + ' R ratio.pdf')

    elif ratio.lower() == 'g':    # temp dependent g ratio
        print("Plotting", name, "G ratio")
        ax2.set_xlabel('Temperature in K', fontsize=plot_settings['fontsize'])
        ax.set_ylabel(name + ' G Ratio', fontsize=plot_settings['fontsize'])

        max_y_tick = 0
        for i, file in enumerate(fnames):
            with fits.open(file) as hdul:
                data = hdul[1].data
                temp_bins, r_orig, r_min, r_max = data['temp'], data['r orig'], data['r min'], data['r max']
                f_orig, f_min, f_max = data['f orig'], data['f min'], data['f max']
                i_orig, i_min, i_max = data['i orig'], data['i min'], data['i max']
                i2_orig, i2_min, i2_max = data['i2 orig'], data['i2 min'], data['i2 max']
                if 'density' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                    string = "{:.0e}".format(data['density'][0]).split("e+")
                    if string[1] == '00':
                        density = "$N_e = " + "{}".format(string[0]) + " cm^{-3}$"
                    else:
                        density = "$N_e = " + "{}\\times 10^{}".format(string[0], round(int(string[1]))) + " cm^{-3}$"
            # do math
            g_orig = (f_orig + i_orig + i2_orig) / r_orig
            numpy.nan_to_num(g_orig)
            # error propagation for positive dg
            dr, df, di, di2 = (r_max - r_orig), (f_max - f_orig), (i_max - i_orig), (i2_max - i2_orig)
            r_term = (((-f_orig - i_orig - i2_orig) / r_orig ** 2) ** 2) * (dr ** 2)
            f_term = df ** 2 / r_orig ** 2
            i_term = di ** 2 / r_orig ** 2
            i2_term = di2 ** 2 / r_orig ** 2
            dg = numpy.sqrt(r_term + f_term + i_term + i2_term)
            g_max = g_orig + dg

            # error propagation for negative dg
            dr, df, di, di2 = (r_orig - r_min), (f_orig - f_min), (i_orig - i_min), (i2_orig - i2_min)
            r_term = (((-f_orig - i_orig - i2_orig) / r_orig ** 2) ** 2) * (dr ** 2)
            f_term = df ** 2 / r_orig ** 2
            i_term = di ** 2 / r_orig ** 2
            i2_term = di2 ** 2 / r_orig ** 2
            dg = numpy.sqrt(r_term + f_term + i_term + i2_term)
            g_min = g_orig - dg

            ax.semilogx(temp_bins, g_orig, color=clist(i))
            if labels[i] == '':
                ax.fill_between(temp_bins, g_min, g_max, alpha=alphas[i], color=clist(i), label=density)
            else:
                ax.fill_between(temp_bins, g_min, g_max, alpha=alphas[i], color=clist(i), label=labels[i])

            error = (abs(g_max-g_min)/g_orig)/2
            error = error[3:]
            if max(error) > max_y_tick: max_y_tick = max(error)
            ax2.semilogx(temp_bins[3:], error, color=clist(i))
            if g_orig[-1] > g_orig[0]: left = True
            else: left = False

        # set y ticks for error subplot if user didn't specify
        if tickstep == {}:
            if max_y_tick < 0.05: tickstep = 0.01
            elif 0.05 < max_y_tick < 0.1: tickstep = 0.02
            elif 0.1 < max_y_tick < 0.2: tickstep = 0.05
            elif 0.2 < max_y_tick < 0.4: tickstep = 0.1
            elif 0.4 < max_y_tick < 0.8: tickstep = 0.2
            elif 0.8 < max_y_tick < 1.5: tickstep = 0.25
            elif 1.5 < max_y_tick: tickstep = 0.5
        ticks = numpy.arange(0, max_y_tick + tickstep, tickstep)
        if plot_settings['important_ticks'] != {}:
            ticks += numpy.array(plot_settings['important_ticks'])
        ax2.set_yticks(ticks)

        if plot_settings['grid_kwargs'] != {}:
            ax2.grid(plot_settings['grid_kwargs'])

        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.95)
        if plot_settings['legendloc'] != 'best':
            ax.legend(fontsize=plot_settings['legendsize'], loc=plot_settings['legendloc'])
        else:
            if left == True: ax.legend(fontsize=plot_settings['legendsize'], loc='upper left')
            else: ax.legend(fontsize=plot_settings['legendsize'], loc='upper right')
        fig.align_ylabels()
        fig.savefig(name + ' G ratio.pdf')

    elif ratio.lower() == 'b':    #plot blended ratio
        print("Plotting", name, "blended ratio")
        ax.set_ylabel('Blended Line Ratio', fontsize=plot_settings['fontsize'])

        max_y_tick = 0
        for i, file in enumerate(fnames):
            if 'dens' in file:
                ax2.set_xlabel('Density in cm$^{-3}$', fontsize=plot_settings['fontsize'])
                with fits.open(file) as hdul:
                    data = hdul[1].data
                    dens_bins, dens_total_min, dens_total_max, dens_total_orig = data['dens'], data['dens_min'], data['dens_max'], data['dens_orig']
                    if 'temperatures' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                        string = "{:.0e}".format(data['temperatures'][0]).split("e+")
                        temperature = "$T_e$ = " + "{}\\times 10^{}$K".format(string[0], round(int(string[1])))
                ax.semilogx(dens_bins, dens_total_orig, color=clist(current))
                if labels[i] == '':
                    ax.fill_between(dens_bins, dens_total_min, dens_total_max, label=temperature, color=clist(i), alpha=alphas[i])
                else:
                    ax.fill_between(dens_bins, dens_total_min, dens_total_max, label=label, color=clist(i), alpha=alphas[i])

                error = (abs(dens_total_max-dens_total_min)/dens_total_orig)/2
                if max(error) > max_y_tick: max_y_tick = max(error)
                ax2.semilogx(dens_bins, error, color=clist(i))
                if dens_total_orig[-1] > dens_total_orig[0]: left = True
                else: left = False

            elif 'Te' in file:
                ax2.set_xlabel('Temperature in K', fontsize=plot_settings['fontsize'])
                with fits.open(file) as hdul:
                    data = hdul[1].data
                    temp_bins, Te_total_min, Te_total_max, Te_total_orig = data['temp'], data['Te_min'], data['Te_max'], data['Te_orig']
                    if 'density' in [hdul[1].header[x] for x in hdul[1].header.keys()]:
                        string = "{:.0e}".format(data['density'][0]).split("e+")
                        density = "$N_e = " + "{}\\times 10^{}".format(string[0], round(int(string[1]))) + " cm^{-3}$"
                ax.semilogx(temp_bins, Te_total_orig)
                if label == '':
                    ax.fill_between(temp_bins, Te_total_min, Te_total_max, label=density, color=clist(i), alpha=alphas[i])
                else:
                    ax.fill_between(temp_bins, Te_total_min, Te_total_max, label=label, color=clist(i), alpha=alphas[i])
                error = (abs(Te_total_max-Te_total_min)/Te_total_orig)/2
                if max(error) > max_y_tick: max_y_tick = max(error)
                ax2.semilogx(temp_bins, error, color=clist(i))
                if Te_total_orig[-1] > Te_total_orig[0]: left = True
                else: left = False

        # set y ticks for error subplot if user didn't specify
        if tickstep == {}:
            if max_y_tick < 0.05: tickstep = 0.01
            elif 0.05 < max_y_tick < 0.1: tickstep = 0.02
            elif 0.1 < max_y_tick < 0.2: tickstep = 0.05
            else: tickstep = 0.1
        ticks = numpy.arange(0, max_y_tick + tickstep, tickstep)
        if plot_settings['important_ticks'] != {}:
            ticks += numpy.array(plot_settings['important_ticks'])
        ax2.set_yticks(ticks)

        if plot_settings['grid_kwargs'] != {}:
            ax2.grid(plot_settings['grid_kwargs'])

        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.95)
        if plot_settings['legendloc'] != 'best':
            ax.legend(fontsize=plot_settings['legendsize'], loc=plot_settings['legendloc'])
        else:
            if left == True: ax.legend(fontsize=plot_settings['legendsize'], loc='upper left')
            else: ax.legend(fontsize=plot_settings['legendsize'], loc='upper right')
        fig.align_ylabels()
        plt.savefig(name + ' blended ratio.pdf')

    if show == True:
        plt.show()

##update for better plotting
def run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range={}, dens_range={}, num={}, pop_fraction={}, show=True, makefiles=False):
    """ Run line diagnostics for specified Z, z1 where up, lo are transition levels at
    specified Te and dens. Vary is 'exc' rate or 'A' value, delta_r is fractional error.
    Can specify range of Te and dens values with tuple of (min, max) values and number
    of values to calculate at. Default Te_range is 20 values over (Te/10, Te*10)
    and default dens_range is 10 values over (1, 1e16). Will plot if set to True.
    Default is both temp and dens diagnostics if left blank. Can pick one type and use
    default ranges by setting either Te_range or dens_range = -1."""
    from timeit import default_timer as timer

    #set default values
    if num == {}: temp_num, dens_num = 20, 10
    else: temp_num, dens_num = num, num

    #set default ranges blank
    if Te_range == -1: Te_range = (Te/10, Te*10)
    if dens_range == -1: dens_range = (1, 1e16)

    #check what type of diagnostics
    if ((Te_range == {}) & (dens_range == {})) or ((Te_range == -1) & (dens_range == -1)): #run both
        type, temperature, density = 'both', True, True
    elif (Te_range == {}) & (dens_range != {}): #run dens diagnostics
        type, temperature, density = 'dens', False, True
    elif (Te_range != {}) & (dens_range == {}): #run temp diagnostics
        type, temperature, density = 'temp', True, False

    start = timer()
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)

    extras = {'process': vary, 'delta_r': delta_r, 'transition': (up, lo), 'transition_2': [],
              'wavelen': (10, 20), 'Te_range': Te_range, 'dens_range': dens_range, 'corrthresh': 10e-5, 'e_signif': 0.0,
              'pop_fraction': {}}

    if temperature == True:
        temp_bins = list(numpy.geomspace(Te_range[0], Te_range[1], num=temp_num))  # 20
        Te_eps_orig, Te_eps_min, Te_eps_max =[], [], []
        for counter, temp_Te in enumerate(temp_bins):
            if vary == 'A':
                Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, extras=extras)
                Te_new_inputs, Te_new_values = vary_a(Te_inputs, Te_values, (up, lo))
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, extras=extras)
                Te_new_inputs, Te_new_values = vary_exc(Te_inputs, Te_values, (lo, up))
            Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
            i = numpy.where((Te_table['Upper'] == up) & (Te_table['Lower'] == lo))[0]
            for col, array in zip(['Epsilon_orig', 'Epsilon_min', 'Epsilon_max'], [Te_eps_orig, Te_eps_min, Te_eps_max]):
                if Te_table[i][col] == 0: Te_table[i][col] = 1e-40
                array.append(Te_table[i][col])
            print(str(temp_num - counter), 'temperatures left\n')

    if density == True:
        dens_bins = list(numpy.geomspace(dens_range[0], dens_range[1], num=dens_num))
        dens_eps_orig, dens_eps_min, dens_eps_max =[],[],[]

        for counter, temp_dens in enumerate(dens_bins):
            print("density is:", temp_dens)
            if vary == 'A':
                dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, extras=extras)
                dens_new_inputs, dens_new_values = vary_a(dens_inputs, dens_values, (up, lo))
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, extras=extras)
                dens_new_inputs, dens_new_values = vary_exc(dens_inputs, dens_values, (lo, up))
            dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)

            i = numpy.where((dens_table['Upper'] == up) & (dens_table['Lower'] == lo))[0]
            for col, array in zip(['Epsilon_orig', 'Epsilon_min', 'Epsilon_max'],
                                  [dens_eps_orig, dens_eps_min, dens_eps_max]):
                if dens_table[i][col] == 0: dens_table[i][col] = 1e-40
                array.append(dens_table[i][col])
            print(str(dens_num-counter), 'densities left\n')

    if type == 'temp':
        line_diagnostics = {'type': 'temp', 'temps': list(temp_bins),'orig': list(Te_eps_orig),
                                 'min': list(Te_eps_min),'max': list(Te_eps_max), 'density': [dens]*len(temp_bins)}
        table = Table([temp_bins, Te_eps_orig, Te_eps_min, Te_eps_max, [dens] * len(temp_bins)],
                      names=('temps', 'orig', 'min', 'max', 'density'))
    elif type == 'dens':
        line_diagnostics = {'type': 'dens', 'dens': list(dens_bins),'orig': list(dens_eps_orig),
                                 'min': list(dens_eps_min),'max': list(dens_eps_max), 'temperature': [Te]*len(dens_bins)}
        table = Table([dens_bins, dens_eps_orig, dens_eps_min, dens_eps_max, [Te]*len(dens_bins)], names=('dens', 'orig', 'min', 'max', 'temperature'))
    elif type == 'both':
        line_diagnostics = {'type': 'both', 'temps': list(temp_bins), 'dens': list(dens_bins),
                                 'Te_orig': list(Te_eps_orig), 'Te_min': list(Te_eps_min),
                                 'Te_max': list(Te_eps_max), 'dens_orig': list(dens_eps_orig),
                                 'dens_min': list(dens_eps_min), 'dens_max': list(dens_eps_max)}
        table = Table([dens_bins, dens_eps_orig, dens_eps_min, dens_eps_max, [Te] * len(dens_bins).
                       temp_bins, Te_eps_orig, Te_eps_min, Te_eps_max, [dens]*len(temp_bins)],
                      names=('dens', 'dens_orig', 'dens_min', 'dens_max', 'temperature',
                             'temps', 'Te_orig', 'Te_min', 'Te_max', 'density'))

    settings = {'Z': Z, 'z1': z1, 'transition': (up, lo), 'vary': vary, 'delta_r': delta_r}

    if makefiles == True:
        fname = element + ' ' + ion + '_' + str(up) + '-' + str(lo) + ' ' + type + ' '
        number = get_file_num(fname, '.fits')
        table.write(fname+str(number) + '.fits')
        print("Wrote to file:", fname+str(number) + '.fits')

    if show == True:
        plot_line_diagnostics(line_diagnostics, settings)
    return line_diagnostics

def plot_line_diagnostics(line_diagnostics, settings):  ##need to update
    """ Line_diagnostics is a dictionary outputted by run_line_diagnostics. It holds arrays of
    the temperature and density bins, the original, minimum, and maximum emissivities from varying
    temperature and then density, the name of the element and ion, label of which process varied
    and by how much, and the transition.

    Plots emissivity as a function of temperature and density for the specified single transition."""

    print("Plotting line diagnostics now.")
    type = line_diagnostics.get('type')
    Z, z1, transition, vary, delta_r = [settings.get(k) for k in settings]
    text = '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1])

    if type == 'temp':  # plot emissivity versus temperature
        type, temps, Te_orig, Te_min, Te_max = [line_diagnostics.get(k) for k in line_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r * 100) + '%')
        ax.set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        ax.semilogx(temps, Te_orig, label='Original', color='b')
        ax.fill_between(temps, Te_min, Te_max, alpha=0.5, color='b', label="Range")
        error = [abs(x - y) / z for x, y, z in zip(Te_max, Te_min, Te_orig)]
        ax2.semilogx(temps, error, color='b')
        ax2.set_xlabel('Temperature in K')
        ax2.set_ylabel('Fractional Error')
        ax.legend(fontsize='xx-small')
    elif type == 'dens':
        type, dens, dens_orig, dens_min, dens_max = [line_diagnostics.get(k) for k in line_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r * 100) + '%')
        ax.set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        ax.semilogx(dens, dens_orig, label='Original', color='b')
        ax.fill_between(dens, dens_min, dens_max, alpha=0.5, color='b', label="Range")
        error = [abs(x-y)/z for x,y,z in zip(dens_max, dens_min, dens_orig)]
        ax2.semilogx(dens, error, color='b')
        ax2.set_xlabel('Density in cm$^{-3}$')
        ax2.set_ylabel('Fractional Error')
        ax.legend(fontsize='xx-small')
    elif type == 'both':
        type, temps, dens, Te_orig, Te_min, Te_max, dens_orig, dens_min, dens_max = [line_diagnostics.get(k) for k in
                                                                                     line_diagnostics]
        fig, ax = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col')
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r * 100) + '%')
        ax[0, 0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        ax[1, 0].set_xlabel('Temperature in K')
        ax[1, 1].set_xlabel('Density in cm$^{-3}$')
        ax[1, 0].set_ylabel('Fractional Error')
        ax[0, 1].semilogx(dens, dens_orig, label='Original', color='b')
        ax[0, 0].semilogx(temps, Te_orig, label='Original', color='b')
        ax[0, 1].fill_between(dens, dens_min, dens_max, alpha=0.5, color='b', label="Range")
        ax[0, 0].fill_between(temps, Te_min, Te_max, color='b', alpha=0.5, label="Range")
        dens_error = [abs(x - y) / z for x, y, z in zip(dens_max, dens_min, dens_orig)]
        Te_error = [abs(x - y) / z for x, y, z in zip(Te_max, Te_min, Te_orig)]
        ax[1,1].semilogx(dens, dens_error, color='b')
        ax[1,0].semilogx(temps, Te_error, color='b')
        for ax in [ax[0, 0], ax[0, 1], ax[1, 1], ax[1, 0]]:
            ax.legend(fontsize='xx-small')

    plt.tight_layout()
    plt.subplots_adjust(top=0.86)
    plt.align_ylabels()
    plt.savefig('line diagnostics.pdf')
    plt.show()
    plt.close('all')

def run_line_ratio_diagnostics(line_diagnostics_1, line_diagnostics_2, settings):
    """ Table1 and table2 are tables from individually run get_tables() on the two transitions
    specified by the user. Inputs1, values1, inputs2, values2, are  dictionaries holding the
    inputs and the sensitivity table/emissivity values for each transition respectively.

    Varies temperature and density separately for each transition and recalculates emissivity.

    Returns dictionary line_ratio_diagnostics containing arrays of the temperature and density bins
    for each transition, as well as the original, minimum, and maximum line ratios calculated
    from varying temperature and density independently."""
    type = line_diagnostics_1.get('type')
    Z, z1, transition1, transition2, vary, delta_r = [settings.get(k) for k in settings]
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)

    if type == 'temp':
        type1, temp_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1, density1 = [line_diagnostics_1.get(k) for k in
                                                                     line_diagnostics_1]
        type2, temp_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2, density2 = [line_diagnostics_2.get(k) for k in
                                                                     line_diagnostics_2]
        Te_line_ratios = [x / y for x, y in zip(Te_eps_orig1, Te_eps_orig2)]
        Te_line_ratios_min = [x / y for x, y in zip(Te_eps_min1, Te_eps_max2)]
        Te_line_ratios_max = [x / y for x, y in zip(Te_eps_max1, Te_eps_min2)]

        line_ratio_diagnostics = {'type': 'temp', 'temps': numpy.asarray(temp_bins1),
                                  'orig': numpy.asarray(Te_line_ratios),
                                  'min': numpy.asarray(Te_line_ratios_min),
                                  'max': numpy.asarray(Te_line_ratios_max), 'density': density1}
        table = Table([temp_bins1, Te_line_ratios, Te_line_ratios_min, Te_line_ratios_max, density1],
                      names=('temps', 'Te_orig', 'Te_min', 'Te_max', 'density'))

    elif type == 'dens':
        type1, dens_bins1, dens_eps_orig1, dens_eps_min1, dens_eps_max1, temp1 = [line_diagnostics_1.get(k) for k in
                                                                           line_diagnostics_1]
        type2, dens_bins2, dens_eps_orig2, dens_eps_min2, dens_eps_max2, temp2 = [line_diagnostics_2.get(k) for k in
                                                                           line_diagnostics_2]
        dens_line_ratios = [x / y for x, y in zip(dens_eps_orig1, dens_eps_orig2)]
        dens_line_ratios_min = [x / y for x, y in zip(dens_eps_min1, dens_eps_max2)]
        dens_line_ratios_max = [x / y for x, y in zip(dens_eps_max1, dens_eps_min2)]

        line_ratio_diagnostics = {'type': 'dens', 'dens': numpy.asarray(dens_bins1),
                                  'orig': numpy.asarray(dens_line_ratios),
                                  'min': numpy.asarray(dens_line_ratios_min),
                                  'max': numpy.asarray(dens_line_ratios_max), 'temperature': temp1}

        table = Table([dens_bins1, dens_line_ratios, dens_line_ratios_min, dens_line_ratios_max, temp1],
                      names=('dens', 'dens_orig', 'dens_min', 'dens_max', 'temperature'))
    elif type == 'both':
        #### add density a nd temperature
        type1, temp_bins1, dens_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1, dens_eps_orig1, \
        dens_eps_min1, dens_eps_max1 = [line_diagnostics_1.get(k) for k in line_diagnostics_1]
        type2, temp_bins2, dens_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2, dens_eps_orig2, \
        dens_eps_min2, dens_eps_max2 = [line_diagnostics_2.get(k) for k in line_diagnostics_2]

        Te_line_ratios = [x / y for x, y in zip(Te_eps_orig1, Te_eps_orig2)]
        Te_line_ratios_min = [x / y for x, y in zip(Te_eps_min1, Te_eps_max2)]
        Te_line_ratios_max = [x / y for x, y in zip(Te_eps_max1, Te_eps_min2)]
        dens_line_ratios = [x / y for x, y in zip(dens_eps_orig1, dens_eps_orig2)]
        dens_line_ratios_min = [x / y for x, y in zip(dens_eps_min1, dens_eps_max2)]
        dens_line_ratios_max = [x / y for x, y in zip(dens_eps_max1, dens_eps_min2)]

        line_ratio_diagnostics = {'type': 'both', 'temps': numpy.asarray(temp_bins1), 'dens': numpy.asarray(dens_bins1),
                                  'Te_orig': numpy.asarray(Te_line_ratios), 'Te_min': numpy.asarray(Te_line_ratios_min),
                                  'Te_max': numpy.asarray(Te_line_ratios_max),
                                  'dens_orig': numpy.asarray(dens_line_ratios),
                                  'dens_min': numpy.asarray(dens_line_ratios_min),
                                  'dens_max': numpy.asarray(dens_line_ratios_max)}
        table = Table([temp_bins1, dens_bins1, Te_line_ratios, Te_line_ratios_min, Te_line_ratios_max,
                       dens_line_ratios, dens_line_ratios_min, dens_line_ratios_max], names=('temps', 'dens',
                    'Te_orig', 'Te_max', 'dens_orig', 'dens_min', 'dens_max'))

    fname = element+' '+ion+'_'+str(transition1[0])+'-'+str(transition1[1])+'-'+str(transition2[0])+'-' + str(transition2[1])+' ratio '+type+' '
    number = get_file_num(fname, '.fits')
    table.write(fname + str(number) + '.fits')
    print("Wrote line diagnostics to file:", fname + str(number) + '.fits')
    fname = fname + str(number) + '.fits'
    return fname, line_ratio_diagnostics

def line_ratio_diagnostics(Z, z1, up, lo, up2, lo2, Te, dens, vary, delta_r, Te_range={}, dens_range={}, pop_fraction={}, num=10, show=True):
    """ Varies 'exc' or 'A' value of each line individually and recalculates
    line ratio with the various errors."""
    # set defaults
    if Te_range == -1: Te_range = (Te / 10, Te * 10)
    if dens_range == -1: dens_range = (10e0, 10e16)
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)

    if dens_range == {}:
        print("Running line ratio diagnostics for", element, ion, "and transition:", up, "->", lo, "/", up2, "->", lo2,
          "over", num, "temperatures in range", Te_range, "with delta_r =", delta_r)
    elif Te_range == {}:
        print("Running line ratio diagnostics for", element, ion, "and transition:", up, "->", lo, "/", up2, "->", lo2,
              "over", num, "densities in range", dens_range, "with delta_r =", delta_r)
    else:
        print("Running line ratio diagnostics for", element, ion, "and transition:", up, "->", lo, "/", up2, "->", lo2,
              "over", num, "temperatures in range", Te_range, "and densities in range", dens_range, "with delta_r =", delta_r)

    if isinstance(delta_r, float): errors = [delta_r]
    else:
        errors = delta_r

    filenames, line_ratio_diags = [], []
    for delta_r in errors:
        settings = {'Z': Z, 'z1': z1, 'transition1': (up, lo), 'transition2': (up2, lo2), 'vary': vary, 'delta_r': delta_r}
        line1 = run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range=Te_range,dens_range=dens_range, num=num, makefiles=False)
        line2 = run_line_diagnostics(Z, z1, up2, lo2, Te, dens, vary, delta_r, Te_range=Te_range,dens_range=dens_range, num=num, makefiles=False)
        fname, line_ratio = run_line_ratio_diagnostics(line1, line2, settings)
        
        if len(errors) == 1:
            plot_ratio(fname, label=str(delta_r*100)+'%', show=show)
            return fname, line_ratio
        else: 
            filenames.append(fname)
            line_ratio_diags.append(line_ratio)
    
    if len(errors) > 1:
        plot_ratio(filenames[::-1], labels=[str(x*100)+'%' for x in errors][::-1], cmap=plot_settings['cmap'], show=show)
        return filenames, line_ratio_diags

def plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, vary, lines, wavelen={}, corrthresh=1e-4, e_signif=1e-20, max_ionize={}, max_recomb={},
                              plot_settings={}, show=True):
    """ Plots sensitive epsilons for multiple transitions all on one plot.
    Lines is dict of {'name': (up, lo)} or just list of [(up, lo), (up2, lo2)]
    Can specify wavelength range as well as correlation threshold and minimum
    emissivity (e_signif) for data plotted. Corrthresh refers to fractional change
    in emissivity, |dE/E|. Vary is either 'A' value or 'exc' rate."""
    if vary == 'A': processes = ['A']
    elif vary == 'exc': processes = ['exc']
    elif vary == 'both': processes = ['A', 'exc']
    if isinstance(lines,list):
        set = {}
        for (up, lo) in lines:
            set.update({str(up)+'->'+str(lo): (up,lo)})
    elif isinstance(lines, dict):
        set = lines

    temp = '%.E' % Decimal(Te) + 'K'
    percentage = '{0:g}%'.format(delta_r * 100)
    density = "{0:.1e}".format(dens) + " cm$^{-3}$"
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)
    clist = get_cmap(len(set) + 2, name=plot_settings['cmap'])

    # try bar plot
    width = 0.25
    x2 = [x + width / 2 for x in numpy.arange(1, len(set)+1)]
    x1 = [x - width / 2 for x in numpy.arange(1, len(set) + 1)]
    fig2, ax2 = plt.subplots(figsize=(10,8))
    ax2.set_xticks(x2)
    ax2.set_xticklabels([x for x in set])
    fig2.suptitle(element+' '+ion+ " lines affected in range \n" + str(wavelen[0]) + "-" + str(wavelen[1])
                  + " $\AA$ at Te =" + temp + ' and ' + density)
    ax2.set_xlabel('Modified Transition by $\pm$ ' + percentage)
    ax2.set_ylabel('Additional Line Emissivities Changed >= ' + '{0:g}%'.format(corrthresh * 100))

    for process, i in zip(processes, range(0,len(processes))):
        if process == 'A': text = ('$\Delta$A value $\pm$ ' + percentage)
        elif process == 'exc': text = ('$\Delta$exc rate $\pm$ ' + percentage)
        fig, ax = plt.subplots(figsize=(10,8))
        ax.set_xlabel('Modified Transition ' + text, fontsize=plot_settings['fontsize'])
        ax.set_ylabel('Fractional Emissivity Change', fontsize=plot_settings['fontsize'])
        for which in ['x', 'y']:
            ax.tick_params(axis=which, labelsize=plot_settings['ticksize'])
        fig.suptitle("Te="+temp + ' dens=' + density)

        t = Table(names=('Upper', 'Lower', 'Lambda', 'Epsilon_orig'))
        plot = {}
        for k, idx in zip(set, numpy.arange(1, len(set)+1)):
            transition, de = set.get(k), []
            if process == 'A':
                extras = {'process': process, 'delta_r': delta_r, 'transition': transition, 'transition_2': None,
                           'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': 0.0, 'e_signif': 0.0,
                          'pop_fraction': {}}
                inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_a(inputs, values, transition)
            elif process == 'exc':
                extras = {'process': process, 'delta_r': delta_r, 'transition': transition[::-1], 'transition_2': None,
                          'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': 0.0, 'e_signif': 0.0}
                inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_exc(inputs, values, transition)

            table, new_table, inputs, results = get_tables(new_inputs, new_values)
            cutoff_data = new_table[new_table['|dE/E|'] >= corrthresh]
            if wavelen != {}:
                cutoff_data = cutoff_data[wavelen[0] < cutoff_data['Lambda']]
                cutoff_data = cutoff_data[cutoff_data['Lambda'] < wavelen[1]]
            if e_signif != {}:
                cutoff_data = cutoff_data[cutoff_data['Epsilon_orig'] >= e_signif]

            #add new affected lines to table
            for line in cutoff_data:
                if (line['Upper'], line['Lower']) not in zip(t['Upper'], t['Lower']):
                    t.add_row([int(line['Upper']), int(line['Lower']), line['Lambda'], line['Epsilon_orig']]+['0']*(len(t.colnames)-4))

            #now compare t to all line data, only add dE for lines in t
            count = 0
            for x in t:
                for line in new_table:
                    if (line['Upper'], line['Lower']) == (x['Upper'], x['Lower']):
                        de.append(line['|dE/E|'])
                        if line['|dE/E|'] != 0 and (line['Upper'], line['Lower']) != transition:
                            count += 1
            new = Column(data=de, name='dE/E ' + process + ' ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1]))
            t.add_columns([new])
            plot.update({k: count})
        print("Number of lines affected in wavelength range and by more than emiss cutoff for", process, "is:", plot)
        #plot only emissivity change of strongest line affected (largest epsilon origs
        for col, idx, k in zip(t.colnames[4:], numpy.arange(1, len(set)+1), set):
            print("\nPlotting lines affected by:", col, " i.e.", k)
            for line in t:
                if line[col] >= corrthresh:
                    if str(int(line['Upper'])) + '->' + str(int(line['Lower'])) == col.split("%")[1][1:]:
                        plt.semilogy(idx, line[col], linestyle='', marker='o', color='k')
                        print("Plotted own transition dE/E=", line[col], "(", line['Upper'], "->", line['Lower'], ")")
                    else:
                        plt.semilogy(idx, line[col], linestyle='', marker='o', color=clist(idx-1))
                        print("Plotted", line['Upper'], "->", line['Lower'], "lambda (", line['Lambda'], "A) with dE/E =", line[col])
        ax.set_xticks(numpy.arange(1, len(set)+1))
        ax.set_xticklabels([x for x in set])
        fig.savefig(element + " " + ion + " sensitive lines.pdf")

        try:
            t.sort('Epsilon_orig', reverse=True)
        except:
            t.sort('Epsilon_orig')
            t = t[::-1]
        print(t)

        fname = 'sensitive lines ' + process + ' ' + element + str(z1) + ' '
        number = get_file_num(fname, '.csv')
        t.write(fname + str(number) + '.csv', format='csv')
        print("Wrote line diagnostics to file:", fname + str(number) + '.csv')

        #plot to bar chart for each process
        if process == 'A':
            ax2.bar(x1, [plot.get(x) for x in plot], width, color=clist(0), label='A value')
        elif process == 'exc':
            ax2.bar(x2, [plot.get(x) for x in plot], width, color=clist(1), label='Direct exc rate')

    ax2.legend(fontsize=plot_settings['fontsize'])
    fig.subplots_adjust(left=0.17, right=0.94, bottom=0.15)
    fig.savefig(fname+'.pdf')
    if show == True:
        plt.show()

def blended_line_ratio(Z, z1, Te, dens, vary, delta_r, transition_list, denom, type={}, num=10, Te_range={}, dens_range={}, max_ionize={}, max_recomb={}, show=True):
    #specify equation of blended line ratio
    #i.e. denom=1, then blended line ratio = [line 1 + line 2] / line 3
    #denm=2, then blended line ratio = line 1 / [line 2 + line 3]

    if type == 'temp':
        print("Running blended line ratio diagnostics for", transition_list, "over", num, "temperatures from Te = ", Te_range)
    elif type == 'dens':
        print("Running blended line ratio diagnostics for", transition_list, "over", num, "densitiess at dens = ", dens_range)

    #input checking
    if type == {}: type = 'both'

    emiss1, emiss2, emiss3 = {}, {}, {}
    emiss_list = (1,2,3)
    extras = {'process': vary, 'delta_r': delta_r, 'transition': [], 'transition_2': [],
              'wavelen': (10, 20), 'Te_range': {}, 'dens_range': {}, 'corrthresh': 10e-5, 'e_signif': 0.0,
              'pop_fraction': {}}
    element = pyatomdb.atomic.Ztoelsymb(Z)

    if type == 'temp':
        if Te_range == {} or Te_range == -1: Te_range = (Te / 10, Te * 10)
        extras.update({'Te_range': Te_range})
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating temperature diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A': extras.update({'transition': transition})
            elif vary == 'exc': extras.update({'transition': transition[::-1]})
            inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
            table = values.get('table')

            lines=[]
            i = numpy.where((table['Upper'] == up) & (table['Lower'] == lo))[1][0]
            lines.append(table[i]['Lambda'])

            diagnostics = run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r,
                            Te_range=Te_range, num=num, plot=False)
            type, temp_bins, Te_eps_orig, Te_eps_min, Te_eps_max, density = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'temp_bins': temp_bins, 'Te_min': numpy.array(Te_eps_min), 'Te_max': numpy.array(Te_eps_max), 
                          'Te_orig': numpy.array(Te_eps_orig), 'density': density}
            if emiss == 2:
                emiss2 = {'temp_bins': temp_bins, 'Te_min': numpy.array(Te_eps_min), 'Te_max': numpy.array(Te_eps_max), 
                          'Te_orig': numpy.array(Te_eps_orig), 'density': density}
            if emiss == 3:
                emiss3 = {'temp_bins': temp_bins, 'Te_min': numpy.array(Te_eps_min), 'Te_max': numpy.array(Te_eps_max), 
                          'Te_orig': numpy.array(Te_eps_orig), 'density': density}

        temp_bins1, Te_eps_min, Te_eps_max, Te_eps_orig, density = [emiss1.get(k) for k in emiss1]
        temp_bins2, Te_eps_min2, Te_eps_max2, Te_eps_orig2, density2 = [emiss2.get(k) for k in emiss2]
        temp_bins3, Te_eps_min3, Te_eps_max3, Te_eps_orig3, density3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            Te_total_min = Te_eps_min/(Te_eps_max2+Te_eps_max3)
            Te_total_max = Te_eps_max/(Te_eps_min2+Te_eps_min3)
            Te_total_orig = Te_eps_orig/(Te_eps_orig2+Te_eps_orig3)
        elif denom==1:
            Te_total_min = (Te_eps_min + Te_eps_min2) / Te_eps_max3
            Te_total_max = (Te_eps_max + Te_eps_min2) / Te_eps_min3
            Te_total_orig = (Te_eps_orig + Te_eps_orig2) / Te_eps_orig3

        t1 = Table([temp_bins1, Te_total_min, Te_total_max, Te_total_orig, density],names=('temp', 'Te_min', 'Te_max', 'Te_orig', 'density'))
        fname = element+ '_'+str(z1) + '_blended_ratio_Te_'
        number = get_file_num(fname, '.fits')
        t1.write(fname + str(number) + '.fits', format='fits')
        print("Wrote blended line diagnostics to file:", fname + str(number) + '.fits')

        plot_ratio(fname + str(number) + '.fits', ratio='b', show=show)

    elif type == 'dens':
        if dens_range == {} or dens_range == -1: dens_range = (1, 1e16)
        extras.update({'dens_range': dens_range})
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating density diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A': extras.update({'transition': transition})
            elif vary == 'exc': extras.update({'transition': transition[::-1]})
            inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
            table = values.get('table')

            lines = []
            i = numpy.where((table['Upper']==up) & (table['Lower']==lo))
            lines.append(table[i]['Lambda'])

            diagnostics = run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, dens_range=dens_range, num=num, plot=False)

            type, dens_bins, dens_eps_orig, dens_eps_min, dens_eps_max, temperature = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'dens_bins': dens_bins, 'dens_min': numpy.array(dens_eps_min), 'dens_max': numpy.array(dens_eps_max),
                          'dens_orig': numpy.array(dens_eps_orig), 'temperature': temperature}
            if emiss == 2:
                emiss2 = {'dens_bins': dens_bins, 'dens_min': numpy.array(dens_eps_min), 'dens_max': numpy.array(dens_eps_max),
                          'dens_orig': numpy.array(dens_eps_orig), 'temperature': temperature}
            if emiss == 3:
                emiss3 = {'dens_bins': dens_bins, 'dens_min': numpy.array(dens_eps_min), 'dens_max': numpy.array(dens_eps_max),
                          'dens_orig': numpy.array(dens_eps_orig), 'temperature': temperature}

        dens_bins, dens_eps_min, dens_eps_max, dens_eps_orig, temperature = [emiss1.get(k) for k in emiss1]
        dens_bins2, dens_eps_min2, dens_eps_max2, dens_eps_orig2, temperature2 = [emiss2.get(k) for k in emiss2]
        dens_bins3, dens_eps_min3, dens_eps_max3, dens_eps_orig3, temperature3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            dens_total_min = dens_eps_min / (dens_eps_max2 + dens_eps_max3)
            dens_total_max = dens_eps_max / (dens_eps_min2 + dens_eps_min3)
            dens_total_orig = dens_eps_orig / (dens_eps_orig2 + dens_eps_orig3)

        elif denom==1:
            dens_total_min = (dens_eps_min + dens_eps_min2) / dens_eps_max3
            dens_total_max = (dens_eps_max + dens_eps_min2) / dens_eps_min3
            dens_total_orig = (dens_eps_orig + dens_eps_orig2) / dens_eps_orig3

        t2 = Table([dens_bins, dens_total_min, dens_total_max, dens_total_orig, temperature], names=('dens', 'dens_min', 'dens_max', 'dens_orig', 'temperature'))
        fname = element+ '_'+str(z1) + '_blended_ratio_dens_'
        number = get_file_num(fname, '.fits')
        t2.write(fname + str(number) + '.fits', format='fits')
        print("Wrote line diagnostics to file:", fname + str(number) + '.fits')

        plot_ratio(fname+str(number)+'.fits', ratio='b', show=show)

    elif type == 'both':    #prob get rid of this option, need to add clist
        if dens_range == {} or dens_range == -1: dens_range = (1, 1e16)
        if Te_range == {} or Te_range == -1: Te_range = (Te/10, Te*10)
        extras.update({'dens_range': dens_range, 'Te_range': Te_range})
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating temperature and density diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A': extras.update({'transition': transition})
            elif vary == 'exc': extras.update({'transition': transition[::-1]})
            inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
            table = values.get('table')

            lines = []
            i = numpy.where((table['Upper'] == up) & (table['Lower'] == lo))[1][0]
            lines.append(table[i]['Lambda'])

            diagnostics = run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r,
                                               Te_range=Te_range, dens_range=dens_range, num=num, plot=False)

            type, temps, dens, Te_orig, Te_min, Te_max, dens_orig, dens_min, dens_max, density, temperature =\
                [diagnostics.get(k) for k in diagnostics]

            if emiss == emiss1:
                emiss1 = {'temp_bins': temps, 'dens_bins': dens, 'Te_min': numpy.array(Te_min), 'Te_max': numpy.array(Te_max),
                          'Te_orig': numpy.array(Te_orig), 'dens_min': numpy.array(dens_min), 'dens_max': numpy.array(dens_max),
                          'dens_orig': numpy.array(dens_orig), 'density': density, 'temperature': temperature}
            if emiss == emiss2:
                emiss2 = {'temp_bins': temps, 'dens_bins': dens, 'Te_min': numpy.array(Te_min), 'Te_max': numpy.array(Te_max),
                          'Te_orig': numpy.array(Te_orig), 'dens_min': numpy.array(dens_min), 'dens_max': numpy.array(dens_max),
                          'dens_orig': numpy.array(dens_orig), 'density': density, 'temperature': temperature}
            if emiss == emiss3:
                emiss3 = {'temp_bins': temps, 'dens_bins': dens, 'Te_min': numpy.array(Te_min), 'Te_max': numpy.array(Te_max),
                          'Te_orig': numpy.array(Te_orig), 'dens_min': numpy.array(dens_min), 'dens_max': numpy.array(dens_max),
                          'dens_orig': numpy.array(dens_orig), 'density': density, 'temperature': temperature}

        temps1, dens1, Te_min, Te_max, Te_orig, dens_min, dens_max, dens_orig, density1, temperature1 = [emiss1.get(k) for k in emiss1]
        temps2, dens2, Te_min2, Te_max2, Te_orig2, dens_min2, dens_max2, dens_orig2, density2, temperature2 = [emiss2.get(k) for k in emiss2]
        temps3, dens3, Te_min3, Te_max3, Te_orig3, dens_min3, dens_max3, dens_orig3, density3, temperature3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            Te_total_min = Te_min/(Te_max2+Te_max3)
            Te_total_max = Te_max/(Te_min2+Te_min3)
            Te_total_orig = Te_orig/(Te_orig2+Te_orig3)
            dens_total_min = dens_min/(dens_max2+dens_max3)
            dens_total_max = dens_max/(dens_min2+dens_min3)
            dens_total_orig = dens_orig/(dens_orig2+dens_orig3)
        elif denom==1:
            Te_total_min = (Te_min + Te_max2) / Te_max3
            Te_total_max = (Te_max + Te_min2) / Te_min3
            Te_total_orig = (Te_orig + Te_orig2) / Te_orig3
            dens_total_min = (dens_min + dens_max2) / dens_max3
            dens_total_max = (dens_max + dens_min2) / dens_min3
            dens_total_orig = (dens_orig + dens_orig2) / dens_orig3

        t1 = Table([temps1, Te_total_min, Te_total_max, Te_total_orig, density1], names=('temp', 'Te_min', 'Te_max', 'Te_orig', 'density'))
        t2 = Table([dens1, dens_total_min, dens_total_max, dens_total_orig, temperature1], names=('dens', 'dens_min', 'dens_max', 'dens_orig', 'temperature'))
        fname1 = element + '_' + str(z1) + '_blended_ratio_Te_'
        number = get_file_num(fname1, '.fits')
        t1.write(fname1+str(number)+'.fits', format='fits')
        fname2 = element+ '_' + str(z1) + '_blended_ratio_dens'
        t2.write(fname2 + str(number) + '.fits', format='fits')
        print("Wrote temperature diagnostics to file:", fname1 + str(number) + '.fits')
        print("Wrote density diagnostics to file:", fname2 + str(number)+ '.fits')

        if show == True:
            blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
            fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig.suptitle(blended_ratio)
            ax1.set_ylabel('Blended ratio', fontsize=plot_settings['fontsize'])
            ax2.set_xlabel('Temperature in K', fontsize=plot_settings['fontsize'])
            ax2.set_ylabel('Fractional Error', fontsize=plot_settings['fontsize'])
            for axis in [ax1, ax2]:
                for which in ['x', 'y']:
                    axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])
            ax1.semilogx(temp_bins, Te_total_orig, label='Original', color='b')
            ax1.fill_between(temp_bins, Te_total_min, Te_total_max, label='Range', color='b', alpha=0.5)
            error = (abs(Te_total_max-Te_total_min)/Te_total_orig)/2
            ax2.semilogx(temp_bins, error, color='b')
            plt.tight_layout()
            fig.subplots_adjust(top=0.86)
            fig.align_ylabels()
            ax1.legend(fontsize=plot_settings['legendsize'])

            fig2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig2.suptitle(blended_ratio)
            ax2.set_xlabel('Density in cm$^{-3}$', fontsize=12)
            ax2.set_ylabel('Fractional Error', fontsize=12)
            ax1.set_ylabel('Blended Ratio', fontsize=12)
            for axis in [ax1, ax2]:
                for which in ['x', 'y']:
                    axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])
            ax1.semilogx(dens_bins, dens_total_orig, label='Original', color='b')
            ax1.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='b', alpha=0.5)
            error = abs(dens_total_max-dens_total_min)/dens_total_orig
            ax2.semilogx(dens_bins, error, color='b')
            plt.tight_layout()
            fig2.subplots_adjust(top=0.86)
            fig2.align_ylabels()
            ax1.legend(fontsize=plot_settings['legendsize'])

            fname1 = element + ' ' + str(z1) + 'blended ratio Te '
            number = get_file_num(fname, '.pdf')
            fig1.savefig(fname1 + str(number) + '.pdf')
            print("Saved temp diagnostics plot to file:", fname1 + str(number) + '.pdf')
            fname2 = element + ' ' + str(z1) + 'blended ratio dens '
            fig2.savefig(fname2 + str(number) + '.pdf')
            print("Saved dens diagnostics plot to file:", fname2 + str(number) + '.pdf')
            plt.show()
            plt.close('all')

def four_line_diagnostic(Z, z1, Te, dens, vary, delta_r, Te_range={}, dens_range={}, type={}, num=10):
    """ Calculates new line emissivities for forbidden,
    resonance, and intercombination lines depending on Z, z1.

    Z: int
        element
    z1: int
        ion charge +1
    Te: int or list/tuple
        single or list of temperatures in K
        will write a different fits file for each Te
    dens: int or list/tuple
        single or list of densities
        will write a different fits file for each dens
    vary: str
        'exc' or 'A'
    delta_r: float
        fractional error to vary rate by, i.e. 0.1 for 10%
    Te_range: list or tuple
        min and max temperature in K
        leave empty or set to -1 for default (Te/10, Te*10)
    dens_range: list or tuple
        min and max densities
        leave empty or set to -1 for default (1, 1e16)
    type: str
        'dens', 'temp', or 'both'
    num: int
        number of temperature or density points
    """
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)

    if type == {}:
        type = 'both'
        if Te_range == -1 or Te_range == {}: Te_range = (Te / 10, Te * 10)
        if (dens_range == -1) or dens_range == {}: dens_range = (1, 1e16)
    elif type == 'temp':
        if Te_range == -1 or Te_range == {}: Te_range = (Te/10, Te*10)
    elif type == 'dens':
        if (dens_range == -1) or dens_range == {}: dens_range = (1, 1e16)

    if (Z % 2) == 0: lines = {'r': (7, 1), 'f': (2, 1), 'i': (6, 1), 'i2': (5, 1)}
    else: lines = {'r': (7, 1), 'f': (2, 1), 'i': (4, 1), 'i2': (5, 1)}

    for line, transition in lines.items():
        up, lo = transition[0], transition[1]
        diagnostics = run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range=Te_range,
                                                              dens_range=dens_range, num=num, makefiles=False)
        if line == 'r': r_diagnostics = diagnostics
        elif line == 'f': f_diagnostics = diagnostics
        elif line == 'i': i_diagnostics = diagnostics
        elif line == 'i2': i2_diagnostics = diagnostics

    # write diagnostics
    if type == 'temp':
        type, temp_bins1, Te_r_orig, Te_r_min, Te_r_max, density = [r_diagnostics.get(k) for k in r_diagnostics]
        type2, temp_bins2, Te_f_orig, Te_f_min, Te_f_max, density2 = [f_diagnostics.get(k) for k in f_diagnostics]
        type3, temp_bins3, Te_i_orig, Te_i_min, Te_i_max, density3 = [i_diagnostics.get(k) for k in i_diagnostics]
        type4, temp_bins4, Te_i_orig2, Te_i_min2, Te_i_max2, density4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temp_bins1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2, density], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max', 'density'))
        fname = element + ' ' + ion + '_four_line_Te_'
        number = get_file_num(fname, '.fits')
        table.write(fname + str(number) + '.fits', format='fits')
        print("Wrote line diagnostics to file:", fname + str(number) + '.fits')
        fname = element + '_' + ion + '_four_line_Te_' + str(number) + '.fits'

    elif type == 'dens':
        type, dens_bins1, dens_r_orig, dens_r_min, dens_r_max, temperature = [r_diagnostics.get(k) for k in r_diagnostics]
        type2, dens_bins2, dens_f_orig, dens_f_min, dens_f_max, temperature2 = [f_diagnostics.get(k) for k in f_diagnostics]
        type3, dens_bins3, dens_i_orig, dens_i_min, dens_i_max, temperature3 = [i_diagnostics.get(k) for k in i_diagnostics]
        type4, dens_bins4, dens_i_orig2, dens_i_min2, dens_i_max2, temperature4 = [i2_diagnostics.get(k) for k in i2_diagnostics]
        
        table = Table([dens_bins1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2, temperature], names=
                       ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max', 'temperature'))
        fname = element + ' ' + ion + '_four_line_dens_'
        number = get_file_num(fname, '.fits')
        table.write(fname + str(number) + '.fits', format='fits')
        print("Wrote line diagnostics to file:", fname + str(number) + '.fits')
        fname = element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits'

    elif type == 'both':
        ##still need to add temperature and density
        t1, temps1, dens1, Te_r_orig, Te_r_min, Te_r_max, dens_r_orig, dens_r_min, dens_r_max, density, temperature = [r_diagnostics.get(k) for k in r_diagnostics]
        t2, temps2, dens2, Te_f_orig, Te_f_min, Te_f_max, dens_f_orig, dens_f_min, dens_f_max, density2, temperature2 = [f_diagnostics.get(k) for k in f_diagnostics]
        t3, temps3, dens3, Te_i_orig, Te_i_min, Te_i_max, dens_i_orig, dens_i_min, dens_i_max, density3, temperature3 = [i_diagnostics.get(k) for k in i_diagnostics]
        t4, temps4, dens4, Te_i_orig2, Te_i_min2, Te_i_max2, dens_i_orig2, dens_i_min2, dens_i_max2, name4, density4, temperature4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temps1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2, density], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max', 'density'))
        table2 = Table([dens1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2, temperature], names=
            ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
             'i2 orig', 'i2 min', 'i2 max', 'temperature'))
        fname1 = element + '_' + str(z1) + '_four_line_Te_'
        fname2 = element + '_' + str(z1) + '_four_line_dens_'
        number = get_file_num(fname, '.fits')
        table2.write(fname1 + str(number) + '.fits', format='fits')
        table.write(fname2 + str(number) + '.fits', format='fits')
        fname = [fname1+str(number)+'.fits', fname2+str(number)+'.fits']
        print("Wrote", fname1, "and", fname2)

    return fname

def g_ratio(Z, z1, dens, vary, delta_r, Te_range={}, num=10, need_data=True, show=True):
    """ Default Te_range is (Te/10, Te*10) where Te = peak ion temperature.
    dens can be a single temp or list densities, different file written each time."""

    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)
    print("Calculating G ratio for", element, ion, "with", str(delta_r*100), "% on", vary)

    if Te_range == {}: Te_range = (Te/10, Te*10)
    if isinstance(dens, (tuple, list)) == False: dens = [dens]

    fnames = []
    if need_data == True:
        for tmp_dens in dens:
            fname = four_line_diagnostic(Z, z1, find_peak_Te(Z, z1), tmp_dens, vary, delta_r, Te_range=Te_range, dens_range={}, type='temp', num=num)
            fnames.append(fname)

    plot_ratio(fnames, ratio='g', show=show)

def r_ratio(Z, z1, Te, vary, delta_r, dens_range={}, num=10, need_data=True, show=True):
    """ Default dens_range is (1, 1e16) with 10 densities.
    Te can be a single temperature in K or a list/tuple of several,
    will write different file for each temp."""

    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)
    print("Calculating R ratio for", element, ion, "with", str(delta_r * 100), "% on", vary)

    if dens_range == {}: dens_range = (1,1e16)
    if isinstance(Te, (tuple, list)) == False: Te = [Te]
    fnames = []

    if need_data == True:
        for tmp_Te in Te:
            fname = four_line_diagnostic(Z, z1, tmp_Te, 1, vary, delta_r, Te_range={}, dens_range=dens_range, type='dens', num=num)
            fnames.append(fname)

    plot_ratio(fnames, ratio='r', show=show)

def solve_ionrec(Telist, ionlist, reclist, Z):
    popn = numpy.zeros([len(Telist), Z + 1])
    for ite in range(len(Telist)):
        Te, ion, rec = Telist[ite], ionlist[ite, :], reclist[ite, :]
        b, a = numpy.zeros(Z + 1, dtype=float), numpy.zeros([Z + 1, Z + 1], dtype=float)

        for iZ in range(0, Z):
            a[iZ, iZ] -= (ion[iZ])
            a[iZ + 1, iZ] += (ion[iZ])

            a[iZ, iZ + 1] += (rec[iZ])
            a[iZ + 1, iZ + 1] -= (rec[iZ])

            # conservation of population
        for iZ in range(0, Z + 1):
            a[0, iZ] = 1.0
        b[0] = 1.0

        c = numpy.linalg.solve(a, b)
        arr = numpy.ndarray.tolist(c)
        for i in range(len(arr)):
            if abs(arr[i]) <= 1e-10:
                arr[i] = 0
        popn[ite, :] = arr
    return popn

def monte_carlo_csd(Z, max_ionize, max_recomb, runs=100, Te_range=[1e4,1e9], distribution='f', Telist={}, makefiles=False, show=True):
    """ Varies CSD with Monte Carlo calculations with a max error on
    ionization and recombination rates. max_ionize and max_recomb can
    either be a float of a single fractional error for all ions,
    or a dict of individual errors with ion z1 being the key,
    i.e. max_ionize = {z1: fractional error}.
    Distribution is either 'f' (default) for flat distribution or 'g' for Gaussian.
    Default is 100 Monte Carlo runs, Te range of 1e4 to 1e9 K, and plot=False.
    If need varied eigenfiles, set makefiles=True."""

    if Telist == {}: Telist = numpy.logspace(numpy.log10(Te_range[0]), numpy.log10(Te_range[1]), 1251)

    element = pyatomdb.atomic.Ztoelsymb(Z)
    clist = get_cmap(Z+2, name=plot_settings['cmap'])

    ionlist, reclist = numpy.zeros([len(Telist), Z]), numpy.zeros([len(Telist), Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True, datacache=d)
        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    eqpopn = get_orig_popn(Z, Telist)

    mc_popn = numpy.zeros([runs, Z+1, len(Telist)])
    random_ion, random_rec = numpy.zeros([len(Telist), Z]), numpy.zeros([len(Telist), Z])

    if isinstance(max_ionize, dict):
        print("Filling in any empty z1 errors as avg of adjacent ions with supplied ionization errors")
        ionize_errors = get_errors(Z, max_ionize)
    elif isinstance(max_ionize, float):
        print("Using fractional error =", max_ionize, "for all ionization rates.")
        ionize_errors = [max_ionize]*(Z)
    elif isinstance(max_ionize, numpy.ndarray): #means another routine already ran get_errors()
        ionize_errors = max_ionize
    if isinstance(max_recomb, dict):
        print("Filling in any empty z1 errors as avg of adjacent ions with supplied recombination errors")
        recomb_errors = get_errors(Z, max_recomb)
    elif isinstance(max_recomb, float):
        print("Using fractional error =", max_recomb, "for all recombination rates.")
        recomb_errors = [max_recomb]*(Z)
    elif isinstance(max_recomb, numpy.ndarray): #means another routine already ran get_errors()
        recomb_errors = max_recomb

    for run in range(runs):
        random1, random2 = [], []
        #get random ionization errors for each ion
        for max_error in ionize_errors:
            if distribution.lower() == 'g':
                lower, upper = -2 * max_error, 2 * max_error
                mu, sigma = 0, max_error
                random1.append(stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=1))
            elif distribution.lower() == 'f':
                random1.append(numpy.random.uniform(low=-max_error, high=max_error, size=1))
        #get random recombination errors
        for max_error in recomb_errors:
            # use average if error not supplied for that ion
            if distribution.lower() == 'g':
                lower, upper = -2 * max_error, 2 * max_error
                mu, sigma = 0, max_error
                random2.append(stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=1))
            elif distribution.lower() == 'f':
                random2.append(numpy.random.uniform(low=-max_error, high=max_error, size=1))
        for col in range(Z):
            random_ion[:,col] = random1[col]
            random_rec[:,col] = random2[col]

        # vary rates at all temps by a different random number for each ion
        iontmp = ionlist * 1.0
        rectmp = reclist * 1.0
        for z1, rand_ion, rand_rec in zip(range(1, Z + 1), random_ion, random_rec):
            rectmp[:, z1 - 1] *= (1 + random_rec[:,z1-1])
            iontmp[:, z1 - 1] *= (1 + random_ion[:,z1-1])
        newpopn = solve_ionrec(Telist, iontmp, rectmp, Z)

        if makefiles==True:
            generate_varied_xspec_ionbal_files(Z, run, iontmp, rectmp)

        for z1 in range(1, Z + 2):
            mc_popn[run, z1 - 1, :] = newpopn[:, z1 - 1]

    # find mean and standard dev from all runs for each ion at each temp
    min = numpy.zeros([len(Telist), Z + 1])
    max = numpy.zeros([len(Telist), Z + 1])
    for z1 in range(1, Z + 2):
        for i in range(len(Telist)):
            pop = mc_popn[:, z1 - 1, i]
            pop_list = numpy.sort(pop)
            min_index = int(0.16 * runs - 1)
            max_index = int(0.84 * runs - 1)
            min[i, z1 - 1] = pop_list[min_index]
            max[i, z1 - 1] = pop_list[max_index]

    if show == True:
        gs_kw = dict(width_ratios=[3], height_ratios=[2, 1])
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True, gridspec_kw=gs_kw)
        ax2.set_xlabel('Temperature (K)', fontsize=plot_settings['fontsize'])
        ax.set_ylabel('Ionic Concentration', fontsize=plot_settings['fontsize'])
        ax.set_ylim(0,1.02)
        ax2.set_ylabel('Fractional error', fontsize=plot_settings['fontsize'])

        for axis in [ax, ax2]:
            for which in ['x', 'y']:
                axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])

        max_y = 0
        for z1 in range(1, Z + 2):
            ax.semilogx(Telist, eqpopn[:, z1 - 1], color=clist(z1 - 1), linestyle='-')
            ax.fill_between(Telist, min[:, z1 - 1], max[:, z1 - 1], color=clist(z1 - 1), alpha=0.4)
            error, temp_bins = [], []
            for x, y, z, Te in zip(max[:, z1 - 1], min[:, z1 - 1], eqpopn[:, z1 - 1], Telist):
                if z >= 10e-2:
                    error.append(((x - y) / z)/2)
                    temp_bins.append(Te)
            ax2.semilogx(temp_bins, error, color=clist(z1 - 1), linewidth=1)
            if numpy.max(error) > max_y: max_y = numpy.max(error)

        ticks = numpy.arange(0, max_y+plot_settings['tickstep'], plot_settings['tickstep'])
        ax2.set_yticks(ticks)
        ax2.set_ylim(0, max_y+0.03)
        for y in ticks[1:-1]:
            ax2.plot(Telist, [y] * len(Telist), color='grey', linewidth=1, alpha=0.3)
        fig.tight_layout()
        fig.subplots_adjust(top=0.86)
        fig.align_ylabels()
        fig.savefig(element + ' Monte Carlo CSD.pdf')
        plt.show()
        plt.close('all')

    return eqpopn, min, max

def temp_change(Z, frac, max_ionize, max_recomb, distribution='f', precision=0.1):
    """
    Z : int
        element
    frac : float
        fraction of peak abundance to calculate temp change
        (i.e. half peak abundance, frac = 0.5)
    max_ionize : float or dict
        maximum ionization error
        if float, max_ionize will be same for all z1 ions
        dict should look like {5: 0.2, 6: 0.3} of form Z : delta_r
    max_recomb : float or dict
        maximum recombination error
        if float, max_recomb will be same for all z1 ions
        dict should look like {5: 0.2, 6: 0.3} of form Z : delta_r
    distribution : str
        type of distribution for error calculation
        'f' for flat, 'g' for Gaussian
    precision : float
        precision in keV of temp shifts given, default is 0.1

    Calculates the absolute (keV) and relative (%)
    temperature shift for a given Z at the specified
    fractional ion concentration after Monte Carlo
    calculations to perturb the CSD. Writes to csv file."""
    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    cols = ['Ion', 'i/r error', 'Te_rise ' + str(frac) + '*peak (keV)', 'Rising Error (keV)', 'Rising Error (%)',
            'Te_fall ' + str(frac) + '*peak (keV)', 'Falling Error (keV)', 'Falling Error (%)']
    num_decimals = len(str(precision/1000).split('.')[1])

    if isinstance(max_ionize, dict):
        ionize_errors = get_errors(Z, max_ionize)
    elif isinstance(max_ionize, float):
        ionize_errors = [max_ionize] * (Z + 1)

    if isinstance(max_recomb, dict):
        recomb_errors = get_errors(Z, max_recomb)
    elif isinstance(max_recomb, float):
        recomb_errors = [max_recomb] * (Z + 1)

    eqpopn, negpopn, pospopn = monte_carlo_csd(Z, ionize_errors, recomb_errors, distribution=distribution)
    fname = element + ' temp change '
    num = get_file_num(fname, '.csv')

    with open(fname+str(num)+'.csv', mode='w') as write_obj:
        dict_writer = csv.DictWriter(write_obj, fieldnames=cols)
        dict_writer.writeheader()
        for z1 in range(1, Z+2):
            if z1 == 1:
                recomb, ionize = 'n/a', str(round(ionize_errors[0]*100))+'%'
                ion = element
            elif z1 == Z+1:
                ionize, recomb = 'n/a', str(round(recomb_errors[-1]*100))+'%'
                ion = element + ' ' + str(z1 - 1) + '+'
            else:
                ionize, recomb = str(round(ionize_errors[z1-1]*100))+'%', str(round(recomb_errors[z1-1]*100))+'%'
                ion = element + ' ' + str(z1 - 1) + '+'

            peak, index = numpy.max(eqpopn[:, z1 - 1]), numpy.argmax(eqpopn[:, z1-1])
            # calculate rising temps first:
            try:
                rising_min_temp = numpy.interp(frac * peak, negpopn[:index, z1 - 1], Telist[:index])/11604525.0061657
                rising_max_temp = numpy.interp(frac * peak, pospopn[:index, z1 - 1], Telist[:index])/11604525.0061657
                rising_orig_temp = round(numpy.interp(frac * peak, eqpopn[:index, z1 - 1], Telist[:index])/11604525.0061657, num_decimals)
                rising_keV_change = round(abs(rising_max_temp - rising_min_temp)/2, num_decimals)
                rising_percent_change = round(((abs(rising_max_temp - rising_min_temp)/rising_orig_temp)/2)*100,1)
            except ValueError:
                print("No rising curve for z1=", z1, "at", frac, "*peak")
                rising_orig_temp, rising_keV_change, rising_percent_change = '--', '--', '--'

            #calculate falling temps
            try:
                falling_min_temp = numpy.interp(frac * peak, negpopn[index:, z1 - 1][::-1], Telist[index:][::-1])/11604525.0061657
                falling_max_temp = numpy.interp(frac * peak, pospopn[index:, z1 - 1][::-1], Telist[index:][::-1])/11604525.0061657
                falling_orig_temp = round(numpy.interp(frac * peak, eqpopn[index:, z1-1][::-1], Telist[index:][::-1])/11604525.0061657, num_decimals)
                falling_keV_change = round(abs(falling_max_temp - falling_min_temp) / 2, num_decimals)
                falling_percent_change = round(((abs(falling_max_temp - falling_min_temp) / falling_orig_temp)/2) * 100, 1)
            except ValueError:
                print("No falling curve on z1 =", z1, "at", frac, "*peak")
                falling_orig_temp, falling_keV_change, falling_percent_change = '--', '--', '--'

            temp_changes = {'Ion': ion, 'i/r error': ionize + ', ' + recomb,
                            'Te_rise ' + str(frac) + '*peak (keV)': rising_orig_temp,
                            'Rising Error (keV)': rising_keV_change, 'Rising Error (%)': rising_percent_change,
                            'Te_fall ' + str(frac) + '*peak (keV)': falling_orig_temp,
                            'Falling Error (keV)': falling_keV_change, 'Falling Error (%)': falling_percent_change}
            dict_writer.writerow(temp_changes)

def plot_temp_shift(Z, z1, max_ionize, max_recomb, frac, show=True):
    """Plots the unperturbed and perturbed CSD from a given
    ionization and recombination error with lines denoting
    the horizontal temperature shift at the specified
    fraction of peak ionic concentration on the rising
    and falling CSD curve.

    Z : int
        element
    z1 : int
        ion charge +1
    max_ionize : float or dict
        float for single ionization error for all ions or
        dict of {z1: error}
    max_recomb : float or dict
        float for single recombination error for all ions or
         dict of {z1: error}
    frac : float
        fraction of peak ionic concentration
        i.e. 0.5 for half peak ion fraction
    show : bool
        True to show plot on screen (default)
        figure saves either way"""

    eqpopn, min, max = monte_carlo_csd(Z, max_ionize, max_recomb, show=False)
    clist, Telist = get_cmap(Z+2, cmap=plot_settings['cmap']), numpy.logspace(4,9,1251)
    idx = numpy.where(eqpopn[:, z1 - 1] >= 10e-2)[0]
    fig, ax = plt.subplots()
    ax.semilogx(Telist[idx], eqpopn[idx, z1 - 1], color=clist(z1 - 1), linestyle='-')
    ax.fill_between(Telist[idx], min[idx, z1 - 1], max[idx, z1 - 1], color=clist(z1 - 1), alpha=0.4)
    ax.tick_params(axis='y', which='minor')
    ax.xaxis.set_minor_formatter(LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(2, 0.4)))
    locs = ax.get_yticks()
    ax.set_ylim(0, locs[-1])

    #interpolate
    frac, peak, index = frac, numpy.max(eqpopn[:, z1 - 1]), numpy.argmax(eqpopn[:, z1 - 1])
    rising_min_temp = numpy.interp(frac * peak, min[:index, z1 - 1], Telist[:index])
    rising_max_temp = numpy.interp(frac * peak, max[:index, z1 - 1], Telist[:index])
    rising_orig_temp = numpy.interp(frac * peak, eqpopn[:index, z1 - 1], Telist[:index])

    falling_min_temp = numpy.interp(frac * peak, min[index:, z1 - 1][::-1], Telist[index:][::-1])
    falling_max_temp = numpy.interp(frac * peak, max[index:, z1 - 1][::-1], Telist[index:][::-1])
    falling_orig_temp = numpy.interp(frac * peak, eqpopn[index:, z1 - 1][::-1], Telist[index:][::-1])

    ax.axhline(frac * peak, color='k', alpha=0.1, linestyle='--')
    ax.errorbar(rising_orig_temp, frac * peak,
                 xerr=[[rising_orig_temp - rising_min_temp], [rising_max_temp - rising_orig_temp]], color='k',
                 marker='.', barsabove=True, capsize=5, linewidth=1.5, markeredgewidth=1.5, zorder=3)
    ax.errorbar(falling_orig_temp, frac * peak,
                 xerr=[[falling_orig_temp - falling_min_temp], [falling_max_temp - falling_orig_temp]], color='k',
                 marker='.', barsabove=True, capsize=5, linewidth=1.5, markeredgewidth=1.5, zorder=3)

    element = pyatomdb.atomic.Ztoelsymb(Z)
    fig.savefig(element+' '+str(z1-1)+'+ temp shift ' + str(frac)+' peak.pdf')
    if show == True:
        plt.show()

def vary_csd(Z, z1, varyir, delta_r, Te_range=[1e4,1e9], errors={}):
    """ Varies CSD by changing z1 'i' or 'r' rate only.
    Delta_r is fractional change, i.e. 0.10 for 10%.
    Plots new CSD for Z.
    Default temperature range is 1e4 to 1e9 K.
    Supply list of positive errors to check rate of change in
    ion concentration as a function of rate error.
    Returns dictionary of orig min, max populations."""

    if plot_settings['latex_font'] == True:
        try:
            setup_text_plots(usetex=True)
        except:
            setup_text_plots(usetex=False)
            print("Using regular matplotlib font for plotting, could not find an installation of LaTeX")

    z1_test, clist = z1, []
    Telist = numpy.logspace(numpy.log10(Te_range[0]), numpy.log10(Te_range[1]), 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)

    eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1_test, varyir, delta_r)
    ret = {'orig': eqpopn, 'max': pospopn, 'min': negpopn}

    fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
    for z1 in range(1, Z + 2):
        if z1 > 1:
            label = element + ' ' + str(z1 - 1) + '+'
        elif z1 == 1:
            label = element
        line, = ax.semilogx(Telist, eqpopn[:, z1 - 1], label=label)
        clist.append(line.get_color())

    for z1 in range(1, Z + 2):
        ax.fill_between(Telist, negpopn[:, z1 - 1], pospopn[:, z1 - 1], color=clist[z1 - 1], alpha=0.5)

        i = numpy.where(eqpopn[:, z1 - 1] > 1e-5)[0]
        if z1 > 1: label = element + ' ' + str(z1 - 1) + '+'
        elif z1 == 1: label = element
        ax2.plot(Telist[i], eqpopn[i, z1 - 1] / negpopn[i, z1 - 1], color=clist[z1 - 1], label=label)
        ax2.plot(Telist[i], eqpopn[i, z1 - 1] / pospopn[i, z1 - 1], color=clist[z1 - 1], label=label)

    #add grid lines to error panel if grid_kwargs were supplied
    if plot_settings['grid_kwargs'] != {}:
        ax2.grid(plot_settings['grid_kwargs'])

    fig.suptitle('CSD from new ' + varyir + ' rate on ' + str(z1_test - 1) + '+')
    ax.set_ylabel("Ionic Fraction", fontsize=plot_settings['fontsize'])
    ax2.set_ylabel("Original/New Fraction", fontsize=plot_settings['fontsize'])
    ax2.set_xlabel("Temperature (K)", fontsize=plot_settings['fontsize'])
    for axis in [ax, ax2]:
        for which in ['x', 'y']:
            axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])
    fig.align_ylabels()
    fig.savefig(element+' '+str(z1_test-1)+'+ ' + varyir + ' ' +str(delta_r*100)+'% CSD.pdf')
    plt.show()
    return ret

def csd_slope(Z, z1, varyir, errors={}, show=True):
    """Plots the rate of change in CSD of z1 and z1+1
    due to perturbing the 'i' or 'r' rate of z1 as a
    function of various errors.
    Z: int
        element
    z1 : int
        ion charge +1
    varyir : str
        rate coefficient to change
        'i' for ionization, 'r' for recombination
    errors : list or tuple
        list of positive errors to vary rate by
    """

    if plot_settings['use_latex_font'] == True:
        try:
            setup_text_plots(usetex=True)
        except:
            setup_text_plots(usetex=False)
            print("Using regular matplotlib font for plotting, could not find an installation of LaTeX")

    z1_test, clist = z1, []
    element = pyatomdb.atomic.Ztoelsymb(Z)

    #do percent changes and slope plot
    fig2, (ax, ax2) = plt.subplots(ncols=2, figsize=(8, 3))
    fig2.subplots_adjust(wspace=0.28)
    ax.set_xlabel('Fractional Error', fontsize=plot_settings['fontsize'])
    ax.set_ylabel('New ' + element + "$^{" + str(z1_test - 1) + "+}$" + '/Original', fontsize=plot_settings['fontsize'])
    ax2.set_xlabel('Fractional Error', fontsize=plot_settings['fontsize'])
    ax2.set_ylabel('New ' + element + "$^{" + str(z1_test) + "+}$" + '/Original', fontsize=plot_settings['fontsize'])
    for axis in [ax, ax2]:
        for which in ['x', 'y']:
            axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])

    percent_names = [' (1%)', ' (10%)', ' (peak)', ' (10%)', ' (1%)']
    # find 5 temperatures of interest
    temp_bins = temps_of_interest(Z, z1_test)
    # find 5 temperatures of interest for neighbor ion
    next_temp_bins = temps_of_interest(Z, z1_test+1)

    if errors == {}: new_errors = [-0.20, -0.15, -0.10, +0.10, +0.15, +0.20]
    else:
        new_errors = errors+[-x for x in errors]
        new_errors.sort()

    frac, next_frac = numpy.zeros([len(temp_bins), len(new_errors)]), numpy.zeros([len(temp_bins), len(new_errors)])
    ionlist, reclist = numpy.zeros([len(temp_bins), Z]), numpy.zeros([len(temp_bins), Z])
    next_ionlist, next_reclist = numpy.zeros([len(next_temp_bins), Z]), numpy.zeros([len(next_temp_bins), Z])

    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(temp_bins, False, Z=Z, z1=z1, extrap=True, datacache=d)
        ionlist[:, z1 - 1], reclist[:, z1-1] = iontmp, rectmp

        next_iontmp, next_rectmp = pyatomdb.atomdb.get_ionrec_rate(next_temp_bins, False, Z=Z, z1=z1, extrap=True, datacache=d)
        next_ionlist[:, z1 - 1], next_reclist[:, z1-1] = next_iontmp, next_rectmp

    eqpopn = solve_ionrec(temp_bins, ionlist, reclist, Z)
    next_eqpopn = solve_ionrec(next_temp_bins, ionlist, reclist, Z)

    for i in range(len(new_errors)):
        delta_r = new_errors[i]
        # copy rates
        iontmp, rectmp = ionlist * 1.0, reclist * 1.0

        # multiply rates by + factor
        if varyir.lower() == 'r': rectmp[:, z1_test - 1] *= (1 + delta_r)
        elif varyir.lower() == 'i': iontmp[:, z1_test - 1] *= (1 + delta_r)
        newpopn = solve_ionrec(temp_bins, iontmp, rectmp, Z)
        next_newpopn = solve_ionrec(next_temp_bins, iontmp, rectmp, Z)

        frac[:, i] = (newpopn[:, z1_test - 1]) / (eqpopn[:, z1_test - 1])
        next_frac[:, i] = (next_newpopn[:, z1_test]) / (next_eqpopn[:, z1_test])

    #plot z1
    for i, percent in zip(range(len(temp_bins)), percent_names):
        Te, y = temp_bins[i], frac[i, :]
        label = numpy.format_float_scientific(Te, exp_digits=1, precision=1) + ' K' + percent
        ax.plot(new_errors, y, label=label, marker='.')

    #plot z1+1
    for i, percent in zip(range(len(next_temp_bins)), percent_names):
        Te, y = next_temp_bins[i], next_frac[i, :]
        label = numpy.format_float_scientific(Te, exp_digits=1, precision=1) + ' K' + percent
        ax2.plot(new_errors, y, label=label, marker='.')

    for axis in (ax, ax2):
        axis.legend(fontsize=plot_settings['legendsize'])
        for which in ['x', 'y']:
            axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])

    plt.tight_layout()
    fig2.savefig(element+' '+str(z1_test-1)+ '+ ' + varyir + ' CSD % change.pdf')

    if show == True:
        plt.show()
        plt.close()

def get_errors(Z, error_dict):
    """
    Z : int
        element
    error_dict : dict
        dictionary of {z1: fractional_error}

    Helper function to fill in any z1 in error dict
    that user didn't supply max ionization or recombination
    error with the average of errors for adjacent ions with
    error given.

    i.e. for Z=8, if error_dict given is {2:0.12, 4:0.2, 5:0.3, 8:0.4}
    routine returns: [0.12 0.12 0.16 0.2  0.3  0.35 0.35 0.4]
    """
    errors = numpy.zeros([Z])
    # first populate error array with existing errors supplied
    for z1, error in error_dict.items():
        errors[z1 - 1] = error

    # now go through and check for z1's with no error supplied
    for z1 in range(1, Z + 1):
        if errors[z1 - 1] == 0:
            # if z1=1 empty, set error equal to z1=2 error
            if z1 == 1:
                errors[z1 - 1] = errors[z1]
            # if z1=Z empty, set error equal to Z-1 ion error
            elif z1 == Z:
                errors[z1 - 1] = errors[z1 - 2]
            # for any other ion, average the two adjacent ions with errors
            else:
                for left_ion in range(z1, 1, -1):
                    if errors[left_ion - 1] != 0:
                        left_error = errors[left_ion - 1]
                        break
                for right_ion in range(z1, Z + 1):
                    if errors[right_ion - 1] != 0:
                        right_error = errors[right_ion - 1]
                        break
                avg = (left_error + right_error) / 2

                # now set any z1 ions with 0 error in between left and right ion to avg
                for tmp_z1 in range(left_ion, right_ion + 1):
                    if errors[tmp_z1 - 1] == 0: errors[tmp_z1 - 1] = avg

    return errors.round(2)

def wrapper_monte_carlo_csd(list, ionize_errors, recomb_errors, runs, makefiles=False, plot=False):
    """List is [], errors is []"""
    for Z in list:
        for max_ionize, max_recomb in zip(ionize_errors, recomb_errors):
            monte_carlo_csd(Z, max_ionize, max_recomb, runs, Telist, makefiles, plot)

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def extremize_csd(Z, max_error, runs, makefiles=False, plot=False):
    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    clist = get_cmap(Z + 2, name=plot_settings['cmap'])

    # set up plot
    plt.figure(figsize=(10,8))
    ax = plt.subplot(111)
    plt.xlabel('Log T(K)', fontsize=12)
    plt.ylabel('Ion Fraction', fontsize=12)

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True, datacache=d)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

    #case 1: ionize + error, recomb - error
    mc_popn = numpy.zeros([runs, Z + 1, len(Telist)])
    random_ion = numpy.zeros([len(Telist), Z])
    random_rec = numpy.zeros([len(Telist), Z])
    for run in range(runs):
        lower, upper = -2 * max_error, 2 * max_error
        mu, sigma = 0, max_error
        random = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
        for col in range(Z):
            random_ion[:, col] = random[col]
            random_rec[:, col] = -random[col]

            # vary rates at all temps by a different random number for each ion
            iontmp = ionlist * 1.0
            rectmp = reclist * 1.0
            for z1, rand_ion, rand_rec in zip(range(1, Z + 1), random_ion, random_rec):
                rectmp[:, z1 - 1] *= (1 + random_rec[:, z1 - 1])
                iontmp[:, z1 - 1] *= (1 + random_ion[:, z1 - 1])
            newpopn = solve_ionrec(Telist, iontmp, rectmp, Z)
            if makefiles == True:
                generate_varied_xspec_ionbal_files(Z, 'pm_'+str(run), iontmp, rectmp)

            for z1 in range(1, Z + 2):
                mc_popn[run, z1 - 1, :] = newpopn[:, z1 - 1]

    if plot == True:
        # find mean and standard dev of each column
        for z1 in range(1, Z + 2):
            median = []
            middle, min, max = [], [], []
            for i in range(len(Telist)):
                pop = mc_popn[:, z1 - 1, i]
                pop_list = numpy.sort(pop)
                median.append(numpy.median(pop_list))
                min_index = int(0.16 * runs - 1)
                max_index = int(0.84 * runs - 1)
                min.append(pop_list[min_index])
                max.append(pop_list[max_index])

            if z1 == 1:
                label = element
            else:
                label = element + ' ' + str(z1 - 1) + '+'
            plt.semilogx(Telist, median, label=label, color=clist(z1 - 1), linestyle='-')
            plt.fill_between(Telist, min, max, color=clist(z1 - 1), alpha=0.4)

        plt.title("Monte Carlo CSD, ionize + " + str(max_error * 100) + "%, recomb - " + str(max_error*100) + "%")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol = 9, fontsize='xx-small')
        plt.savefig(element + ' i+ r- ' + str(max_error*100) +'% Monte Carlo CSD.pdf')
        plt.close('all')

    # case 2: ionize - error, recomb + error
    plt.figure(figsize=(10,8))
    ax = plt.subplot(111)
    plt.xlabel('Log T(K)', fontsize=12)
    plt.ylabel('Ion Fraction', fontsize=12)

    mc_popn = numpy.zeros([runs, Z + 1, len(Telist)])
    random_ion = numpy.zeros([len(Telist), Z])
    random_rec = numpy.zeros([len(Telist), Z])
    for run in range(runs):
        lower, upper = -2 * max_error, 2 * max_error
        mu, sigma = 0, max_error
        random = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
        print("Random errors should be all pos fractions", random)
        for col in range(Z):
            random_ion[:, col] = -random[col]
            random_rec[:, col] = +random[col]

            # vary rates at all temps by a different random number for each ion
            iontmp = ionlist * 1.0
            rectmp = reclist * 1.0
            for z1, rand_ion, rand_rec in zip(range(1, Z + 1), random_ion, random_rec):
                rectmp[:, z1 - 1] *= (1 + random_rec[:, z1 - 1])
                iontmp[:, z1 - 1] *= (1 + random_ion[:, z1 - 1])
            newpopn = solve_ionrec(Telist, iontmp, rectmp, Z)
            if makefiles == True:
                generate_varied_xspec_ionbal_files(Z, 'mp_'+str(run), iontmp, rectmp)

            for z1 in range(1, Z + 2):
                mc_popn[run, z1 - 1, :] = newpopn[:, z1 - 1]

    if plot==True:
        # find mean and standard dev of each column
        for z1 in range(1, Z + 2):
            median = []
            middle, min, max = [], [], []
            for i in range(len(Telist)):
                pop = mc_popn[:, z1 - 1, i]
                pop_list = numpy.sort(pop)
                median.append(numpy.median(pop_list))
                min_index = int(0.16 * runs - 1)
                max_index = int(0.84 * runs - 1)
                min.append(pop_list[min_index])
                max.append(pop_list[max_index])
            if z1 == 1:
                label = element
            else:
                label = element + ' ' + str(z1 - 1) + '+'
            plt.semilogx(Telist, median, label=label, color=clist(z1 - 1), linestyle='-')
            plt.fill_between(Telist, min, max, color=clist(z1 - 1), alpha=0.4)

        plt.title("Monte Carlo CSD, ionize - " + str(max_error * 100) + "%, recomb + " + str(max_error * 100) + "%")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol = 9, fontsize='xx-small')
        plt.savefig(element + ' i- r+ ' + str(max_error*100) +'% Monte Carlo CSD.pdf')
        plt.close('all')

def generate_varied_xspec_ionbal_files(Z, filesuffix, ionlist, reclist):
    import scipy.linalg

    Telist = numpy.logspace(4, 9, 1251)

    # outputs:
    feqb = numpy.zeros([len(Telist), Z + 1])
    vl_out = numpy.zeros([len(Telist), Z ** 2])
    vr_out = numpy.zeros([len(Telist), Z ** 2])
    eig_out = numpy.zeros([len(Telist), Z])

    for ite in range(len(Telist)):
        Te = Telist[ite]
        ion = ionlist[ite, :]
        rec = reclist[ite, :]

        b = numpy.zeros(Z + 1, dtype=float)
        a = numpy.zeros([Z + 1, Z + 1], dtype=float)

        for iZ in range(0, Z):
            a[iZ, iZ] -= (ion[iZ])
            a[iZ + 1, iZ] += (ion[iZ])

            a[iZ, iZ + 1] += (rec[iZ])
            a[iZ + 1, iZ + 1] -= (rec[iZ])

        # conservation of population
        for iZ in range(0, Z + 1):
            a[0, iZ] = 1.0
        b[0] = 1.0

        c = numpy.linalg.solve(a, b)
        c[0] = 1 - sum(c[1:])
        c[c < 1e-10] = 0.0
        feqb[ite] = c

        ZZ = len(ion) + 1
        ndim = ZZ
        AA = numpy.zeros((ndim - 1, ndim - 1), dtype=float)
        # populate with stuff

        for iCol in range(ndim - 1):
            for iRow in range(ndim - 1):

                if (iRow == 0):
                    if (iCol == 0):
                        if (Z >= 2):
                            AA[0, iCol] = -(ion[0] + ion[1] + rec[0])
                        else:
                            AA[0, iCol] = -(ion[0] + rec[0])

                    if (iCol == 1): AA[0, iCol] = rec[1] - ion[0]
                    if (iCol > 1):
                        AA[0, iCol] = -ion[0]
                else:
                    if (iRow == iCol + 1):  AA[iRow, iCol] = ion[iRow]
                    if (iRow == iCol):
                        if (iRow + 2 < ndim):

                            AA[iRow, iCol] = -(ion[iRow + 1] + rec[iRow])
                        else:
                            AA[iRow, iCol] = -rec[iRow]

                    if (iRow == iCol - 1):
                        AA[iRow, iCol] = rec[iRow + 1]

        w, vr = numpy.linalg.eig(AA)
        if (w.dtype != 'float64'):
            print("nooooooooooO", w.dtype)
        leftevec = numpy.zeros(Z ** 2)
        rightevec = numpy.zeros(Z ** 2)

        # The value VL in which is stored is not actually the left eigenvecotr,
        # but is instead the inverse of vr.

        vl = numpy.matrix(vr) ** -1

        for i in range(Z):
            for j in range(Z):
                leftevec[i * Z + j] = vl[i, j]
                rightevec[i * Z + j] = vr[j, i]

        vr_out[ite] = rightevec
        vl_out[ite] = leftevec
        eig_out[ite] = w

    hdu0 = pyfits.PrimaryHDU()
    now = datetime.datetime.utcnow()

    hdu0.header['DATE'] = now.strftime('%d/%m/%y')
    hdu0.header['FILENAME'] = "Python routine"
    hdu0.header['ORIGIN'] = ("ATOMDB", os.environ['USER'] + ", AtomDB project")

    # secondary HDU, hdu1:
    hdu1 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(
        [pyfits.Column(name='FEQB',
                       format='%iD' % (Z + 1),
                       array=feqb),
         pyfits.Column(name='EIG',
                       format='%iD' % (Z),
                       array=eig_out),
         pyfits.Column(name='VR',
                       format='%iD' % (Z ** 2),
                       array=vr_out),
         pyfits.Column(name='VL',
                       format='%iD' % (Z ** 2),
                       array=vl_out)]))

    hdulist = pyfits.HDUList([hdu0, hdu1])
    hdulist[1].header['EXTNAME'] = 'EIGEN'

    element = pyatomdb.atomic.Ztoelsymb(Z).lower()
    fname = 'varied_eigen_%s_%s.fits' % (element, filesuffix)
    hdulist.writeto(fname, checksum=True, overwrite=True)

def update_eigen(Z, filename, session):
    """
    Updates eigenvector file data for element Z from file filename

    """
    d = pyatomdb.pyfits.open(filename)

    session.spectra.datacache['data']['misc']['EIGEN'][Z] = d
    session.spectra.datacache['datasums']['misc']['EIGEN'][Z] = d[1].header['DATASUM']

def extremize_ionreclist(Z, Telist, ion_deltar={}, rec_deltar={}):
    """Will shift all rates by specified delta_r. Can provide
    delta_r for ionization and/or recombination rate, will shift
    accordingly and return new ion and rec lists."""

    #vary all ions of element Z, vary either ion rate, rec rate or both

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    #get original rates
    for z1 in range(1, Z + 1):
        iontmp, rectmp = atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True,
                                                settings=False, datacache=d)
        if ion_deltar != {}:
            #change ion rate
            iontmp = [x + (x * (ion_deltar / 100)) for x in iontmp]

        if rec_deltar != {}:
            #change rec rate
            rectmp = [x + (x * (rec_deltar / 100)) for x in rectmp]

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    return ionlist, reclist

def gaussian_varied_rates(Z, max_error, Telist):
    """Applies random error from Gaussian distribution
    of sigma = max_error for all rate coefficients
    in an element over Telist range.
    Returns arrays of varied ion and rec rates."""

    #get original rates
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True, datacache=d)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    #get random errors from truncated Gaussian
    random_ion = numpy.zeros([len(Telist), Z])
    random_rec = numpy.zeros([len(Telist), Z])
    lower, upper = -2 * max_error, 2 * max_error
    mu, sigma = 0, max_error

    random1 = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
    random2 = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
    for col in range(Z):
        random_ion[:, col] = random1[col]
        random_rec[:, col] = random2[col]

    # vary rates at all temps by a different random number for each ion
    iontmp = ionlist * 1.0
    rectmp = reclist * 1.0
    for z1, rand_ion, rand_rec in zip(range(1, Z + 1), random_ion, random_rec):
        rectmp[:, z1 - 1] *= (1 + random_rec[:, z1 - 1])
        iontmp[:, z1 - 1] *= (1 + random_ion[:, z1 - 1])

    return iontmp, rectmp

#still need to fix
def get_partial_deriv(Z, z1, vary, errors, Te={}, dens=1, lines={}, wavelen={}, corrthresh=10e-4, e_cutoff={}):
    """Creates .csv file with partial log derivatives
    dE/dR and dE/E for specified error(s) for key lines
    and other affected lines.

    dens: int or str 'critical'
        If no dens supplied, uses dens=1
        If 'critical', will use critical dens
    errors: list or float
        float for single error
        list of multiple errors to run for
    vary: str of either 'exc', 'A', or 'both'
    Te: float
        If no Te supplied, uses Te where ion peaks in plasma
    lines: list []
        If key lines not specified,
        will vary H-like, He-like, or Fe line complexes.
        Can either list lines as (up, lo) or int of wavelength in A
    wavelen: tuple or list
        min and max wavelength in A for range
    corrthresh : int
        Minimum change in emissivity from perturbing line
        Default is 10e-4
    e_cutoff: int
        Minimum emissivity cutoff after factoring in elemental abundance.
        i.e. 1e-17."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    print("Calculating partial derivatives for", element, str(z1))

    if Te == {}:
        #find temperature where ion fraction peaks
        Te = find_peak_Te(Z, z1, unit='K')
        print("Using peak Te =", Te)

    if lines == {}:
        if (Z, z1) == (26,17):
            #trans_list = [(2,1), (3,1), (5,1), (17, 1), (23,1), (27,1)]
            trans_list = [(17, 1), (23,1), (27,1)]
        if (Z, z1) == (26, 19):
            trans_list = [(53,1), (68,1), (71,1), (74, 1), (76, 4)]
        if Z-z1 == 1:
            trans_list = [(2,1), (5,1), (6,1), (7,1), (13,1)]
            print("Varying He-like lines.")
        if Z==z1:
            print("Varying H-like lines.")
            trans_list=[(3,1), (4,1), (6,1), (7,1)]
    elif lines != {}: trans_list = lines

    if dens == 'critical':
        # define critical densities for ions
        critical_dens = {}
        elsymb = pyatomdb.atomic.Ztoelsymb(Z)
        if (elsymb, z1) == ('O', 7): critical_dens = 10e14
        if (elsymb, z1) == ('N', 6): critical_dens = 10e11
        if (elsymb, z1) == ('Ne', 9): critical_dens = 10e13
        if (elsymb, z1) == ('Mg', 11): critical_dens = 10e14
        if (elsymb, z1) == ('Si', 13): critical_dens = 10e15
        if (elsymb, z1) == ('S', 15): critical_dens = 10e15
        if (elsymb, z1) == ('Al', 12): critical_dens = 10e14
        if (elsymb, z1) == ('Fe', 17): critical_dens = 10e15
        if (elsymb, z1) == ('Fe', 19): critical_dens = 10e14

        if critical_dens != {}:
            print("Critical density found, will vary at dens =", critical_dens)
            dens = critical_dens
        elif critical_dens == {}:
            print("No critical density found, exiting.")
            exit()

    if vary == 'both': vary_list = ['exc', 'A']
    else: vary_list = [vary]

    if isinstance(errors, float): errors = [errors]

    uppers, lowers = [int(x[0]) for x in trans_list], [int(x[1]) for x in trans_list]

    inputs, values = set_up(Z, z1, Te, dens)
    table = values.get('table')

    # set up output tables
    wavelengths, orig = [],[]
    for x in trans_list:
        i = numpy.where([(table['Upper'] == x[0]) & (table['Lower'] == x[1])])[1][0]
        wavelengths.append(table[i]['Lambda'])
        orig.append(table[i]['Epsilon_orig'])
    t = Table([uppers, lowers, wavelengths, orig], names=('Upper', 'Lower', 'Lambda', 'Epsilon'))
    t2 = Table([uppers, lowers, wavelengths, orig], names=('Upper', 'Lower', 'Lambda', 'Epsilon'))
    t3 = Table(names=('Upper', 'Lower', 'Lambda', 'Epsilon'))
    t4 = Table(names=('Upper', 'Lower', 'Lambda', 'Epsilon'))
    t5 = Table(names=('Upper', 'Lower', 'Lambda', 'Epsilon'))
    t6 = Table(names=('Upper', 'Lower', 'Lambda', 'Epsilon'))

    for delta_r in errors:
        for transition in trans_list:
            print("\nVarying the transition", str(transition[0])+'->'+str(transition[1]), 'by', str(delta_r*100)+'%')
            for rate_type in vary_list:
                if rate_type == 'exc': #vary exc rate
                    extras={'process':'exc', 'delta_r': delta_r, 'transition':transition[::-1], 'transition_2':[],
                         'wavelen':wavelen, 'Te_range':(Te/10, Te*10), 'dens_range':(1,10e16),
                            'corrthresh':0.0, 'e_signif':0.0, 'pop_fraction':{}}
                    dR = 'dE/dR exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                    dE = 'dE/E exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                    inputs, values, exc_transition = set_up(Z, z1, Te, dens, extras=extras)
                    new_inputs, new_values = vary_exc(inputs, values, exc_transition)
                    table, new_table, inputs, results = get_tables(new_inputs, new_values)
                elif rate_type == 'A': #vary A value
                    extras = {'process': 'A', 'delta_r': delta_r, 'transition': transition, 'transition_2': [],
                               'wavelen': wavelen, 'Te_range': (Te / 10, Te * 10), 'dens_range': (1, 10e16),
                              'corrthresh': 0.0, 'e_signif': 0.0, 'pop_fraction':{}}
                    dR = 'dE/dR A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                    dE = 'dE/E A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                    inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                    new_inputs, new_values = vary_a(inputs, values, transition)
                    table, new_table, inputs, results = get_tables(new_inputs, new_values)

                partial_deriv, frac_E, linear_partial_deriv, linear_frac_E = [], [], [], []
                # get derivatives for the strong lines varied
                for x in t:
                    i = numpy.where([(new_table['Upper'] == x['Upper']) & (new_table['Lower'] == x['Lower'])])[1][0]
                    deriv, frac = new_table[i]['dE/dR'], new_table[i]['|dE/E|']
                    frac_E.append(frac)
                    partial_deriv.append(deriv)

                t[dR], t2[dE] = partial_deriv, frac_E

                # check for wavelength and emissivity cutoffs:
                if e_cutoff != {}:
                    ab = pyatomdb.atomdb.get_abundance()[Z]
                    table = table[e_cutoff <= table['Epsilon_orig'] * ab]
                    new_table = new_table[e_cutoff <= new_table['Epsilon_orig'] * ab]
                if wavelen != {}:
                    table = table[(wavelen[0] <= table['Lambda']) & (table['Lambda' <= wavelen[1]])]
                    new_table = new_table[(wavelen[0] <= new_table['Lambda']) & (new_table['Lambda' <= wavelen[1]])]

                #get partial derivs for lines in "interesting lines" and linear lines tables:
                partial_deriv, frac_E, linear_partial_deriv, linear_frac_E = [], [], [], []
                if len(t4) != 0:
                    for x in t4:
                        i = numpy.where([(new_table['Upper'] == x['Upper']) & (new_table['Lower'] == x['Lower'])])[1][0]
                        frac_E.append(new_table[i]['|dE/E|'])
                        partial_deriv.append(new_table[i]['dE/dR'])

                if len(t6) != 0:
                    for x in t6:
                        i = numpy.where([(new_table['Upper'] == x['Upper']) & (new_table['Lower'] == x['Lower'])])[1][0]
                        linear_frac_E.append(new_table[i]['|dE/E|'])
                        linear_partial_deriv.append(new_table[i]['dE/dR'])

                #now search through tables for other lines affected and add to tables
                new_table2 = new_table[(corrthresh <= new_table['|dE/E|'])]
                print("Affected lines:\n", new_table2)
                for x in new_table2:
                    if 0.98 <= x['|dE/E|'] and x['|dE/E|'] <= 1.02:
                        if (x['Upper'], x['Lower']) not in (t6 and trans_list):
                            row = [x['Upper'], x['Lower'], x['Lambda'], x['Epsilon_orig']] + ['0'] * (len(t6.colnames) - 4)
                            t5.add_row(row)
                            t6.add_row(row)
                        linear_partial_deriv.append(x['dE/dR'])
                        linear_frac_E.append(x['|dE/E|'])
                    else:
                        if (x['Upper'], x['Lower']) not in (t4 and trans_list):
                            row = [x['Upper'], x['Lower'], x['Lambda'], x['Epsilon_orig']] + ['0'] * (len(t4.colnames) - 4)
                            t3.add_row(row)
                            t4.add_row(row)
                        partial_deriv.append(x['dE/dR'])
                        frac_E.append(x['|dE/E|'])
                print("Partial deriv:", partial_deriv, '\ndE', frac_E)
                # interesting = [x for x in new_table2 if ((0.98 <= x['|dE/E|']) and (x['|dE/E|'] <= 1.02)) == False]
                # new = [x for x in interesting if (x['Upper'], x['Lower']) not in t4]
                # new_table2 = [x for x in new_table2 if (x['Upper'], x['Lower']) not in trans_list]
                # df = new_table2.to_pandas()
                #
                # interesting = df[numpy.logical_and(0.98 <= df['|dE/E|'], df['|dE/E|'] <= 1.02) == False]
                # new = interesting[numpy.logical_and(interesting['Upper'] not in t4['Upper'], interesting['Lower'] not in t4['Lower'])]
                # print(new)
                # rows = [row['Upper', 'Lower', 'Lambda', 'Epsilon_orig'] + ['0'] * (len(t4.colnames) - 4) for i,row in new.iterrows()]
                # if len(rows) > 0:
                #     t3.add_row(rows)
                #     t4.add_row(rows)
                # partial_deriv += list(new['dE/dR'])
                # frac_E += list(new['|dE/E|'])
                #
                # linear = [x for x in new_table2 if (0.98 <= x['|dE/E|']) and (x['|dE/E|'] <= 1.02)]
                # new = [x for x in linear if (x['Upper'], x['Lower']) not in t6]
                # print(new)
                # rows = [x[:,:4] + ['0'] * (len(t6.colnames) - 4) for x in new]
                # if len(rows) > 0:
                #     t5.add_row(rows)
                #     t6.add_row(rows)
                # linear_partial_deriv += list(new['dE/dR'])
                # linear_frac_E += list(new['|dE/E|'])

                t3[dR], t4[dE] = partial_deriv, frac_E
                t5[dR], t6[dE] = linear_partial_deriv, linear_frac_E

    #sort by epsilon descending
    try:
        for x in [t3,t4,t5,t6]:
            x.sort(['Epsilon'], reverse=True)
    except TypeError:
        for x in [t3,t4,t5,t6]:
            x.sort('Epsilon')
            x = x[::-1]

    #replace row separating lines changed with "other affected lines" with row of --
    t2.add_row([0] * len(t2.colnames), mask=[True]*len(t2.colnames))
    
    #if t6 is not empty, add separating line of -- between "interesting lines" t4 and "linear lines" t6
    if len(t6) != 0:
        t4.add_row([0] * len(t2.colnames), mask=[True]*len(t2.colnames))

    #join t2 and t4 and t6
    joined = vstack([t2, t4, t6])

    #write csv files
    fname = element+' '+str(z1)+' partial deriv '
    number = get_file_num(fname, '.csv')
    joined.write(fname + str(number) + '.csv', format='csv')
    print("Wrote to", fname + str(number) + '.csv')

def wrapper_get_partial_deriv(list, errors, type, dens={}):
    for (Z, z1) in list:
        for error in errors:
            if dens == {}: dens=[1]
            for ne in dens:
                if type(Z) == str:
                    get_partial_deriv(pyatomdb.atomic.elsymb_to_Z(Z), z1, [error], type, ne)
                elif type(Z) == int:
                    get_partial_deriv(Z, z1, [error], type, ne)

def find_peak_Te(Z, z1, unit='K'):
    """ Returns temperature where ion peaks.
    Unit is str of 'K' or 'keV'.
    Default unit is K."""

    # find temperature where ion fraction peaks
    Telist = numpy.logspace(4, 9, 1251)
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for temp_z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=temp_z1, extrap=True,settings=False, datacache=d)
        ionlist[:, temp_z1 - 1] = iontmp
        reclist[:, temp_z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)
    idx = numpy.argmax(eqpopn[:, z1-1])
    peak_Te = Telist[idx]

    if unit == 'keV':
        peak_Te = peak_Te/11604525.0061657

    return peak_Te

 #in progress I think
def get_ionfrac(Z, Te, delta_r, varyir=False, Teunit='K', z1=[]):
    """
      Calculates ionization fraction of a all ions or
      specific z1 at a given Te from CSD after varying all
      rate coefficients (varyir=False) or just one
      rate coefficient (varyir='i' or 'r').
      Assumes ionization equilibrium.

      Parameters
      ----------
      Z : int
        atomic number of element (e.g. 6 for carbon)
      Te : float or list
        electron temperature
      varyir : bool
        varies all rate coefficients if False
        varies z1 rate if set to 'i' or 'r'
      Teunit : str
        'K' or 'keV'
      z1 : int
        if provided, z+1 of ion (e.g. 5 for O V). If omitted, returns ionization
        fraction for all ions of element. If varyir != False, must provide z1

      Returns
      -------

      dictionary of min, median, and max ionization fractions
       of all ions or specific z1 at Te from Monte Carlo CSD

      """
    Telist = numpy.logspace(4,9,1251)

    #check what type of CSD varying
    if varyir == False:
        avg, low, high = monte_carlo_csd(Z, max_ionize, max_recomb)
    elif (varyir != False) and (z1 == []):
        print("Must specify z1 to change rate for, exiting.")
        exit()
    else:
        avg, low, high = get_new_popns(Z, Telist, delta_r, vary=varyir)

    #check what temperature
    if isinstance(Te, int): Te = [Te]
    if unit.lower() == 'keV':
        Te = [x/11604525.0061657 for x in Te]

    median, min, max = numpy.zeros([len(Te), Z+1]), numpy.zeros([len(Te), Z+1]), numpy.zeros([len(Te), Z+1])

    for i in range(len(Te)):
        for tmp_z1 in range(Z+1):
            median[i, tmp_z1-1] = numpy.interp(Te, Telist, avg[:, tmp_z1-1])
            min[i, tmp_z1-1] = numpy.interp(Te, Telist, low[:, tmp_z1 - 1])
            max[i, tmp_z1-1] = numpy.interp(Te, Telist, high[:, tmp_z1 - 1])
    if z1 != []:
        ret = {'min': min[:, z1-1], 'median': median[:, z1-1], 'max': max[:, z1-1]}
    else:
        ret = {'min': min, 'median': median, 'max': max}
    return ret

def affected_lines(Z, z1, up, lo, Te, dens, vary, delta_r, corrthresh=1e-4, e_signif=1e-20, wavelen={}, makefiles=True, observed_only=False):
    """ Varies Z, z1 'exc' or 'A' by delta_r at specified Te and dens
    and returns table of affected lines with original epsilon, dE/dR and dE/E,
    sorted by greatest change in emissivity dE/E to smallest.
    Writes table to csv file if makefiles=True."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    extras = {'process':vary, 'delta_r':delta_r,'transition':(up, lo), 'transition_2': [],
            'wavelen': {}, 'Te_range':{},'dens_range': {},'corrthresh':0, 'e_signif':0.0, 'pop_fraction':{}}
    if vary == 'exc':
        extras.update({'transition': (lo, up)})
        print("Transition is", extras['transition'])
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_exc(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)
    elif vary == 'A':
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_a(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)

    try:
      new_table.sort(['|dE/E|', 'Epsilon_orig'], reverse=True)
    except TypeError:
      new_table.sort(['|dE/E|', 'Epsilon_orig'])
      new_table = new_table[::-1]

    new_table = new_table[new_table['|dE/E|'] >= corrthresh]
    if wavelen != {}:
        new_table = new_table[wavelen[0] < new_table['Lambda']]
        new_table = new_table[new_table['Lambda'] < wavelen[1]]
    if e_signif != {}:
        new_table = new_table[new_table['Epsilon_orig'] >= e_signif]

    if observed_only == True:
        observed = Table(names=[x for x in new_table.colnames])
        ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
        in_range = ladat[1].data
        for row in new_table:
            for x in in_range:
                if (x['UPPER_LEV'], x['LOWER_LEV']) == (row['Upper'], row['Lower']):
                    lambda_obs, lambda_theory = x['WAVE_OBS'], x['WAVELEN']
                    if numpy.isnan(lambda_obs) == False: observed.add_row([row[x] for x in new_table.colnames])

    if makefiles == True:
        fname = element + ' ' + str(z1) + ' (' + str(up) +',' +str(lo) + ') affected lines '
        number = get_file_num(fname, '.csv')
        if observed_only == True: observed.write(fname+str(number)+'.csv', format='csv')
        else: new_table.write(fname+str(number)+'.csv', format='csv')
        print("Wrote to", fname + str(number) + '.csv')

    print(new_table)
    return new_table

def new_epsilons(Z, z1, up, lo, Te, dens, vary, delta_r, e_signif=1e-20, wavelen={}, makefiles=False):
    """ Varies Z, z1 'exc' or 'A' by delta_r at specified Te and dens
    and returns table of new emissivities for lines with original epsilon,
    minimum and maximum epsilons from the varied rate. Table is
    sorted by greatest epsilon to smallest.
    Writes table to csv file if makefiles=True."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    extras = {'process':vary, 'delta_r':delta_r,'transition':(up, lo), 'transition_2': [],
            'wavelen': {}, 'Te_range':{},'dens_range': {},'corrthresh':0, 'e_signif':e_signif}
    if vary == 'exc':
        extras.update({'transition': (lo, up)})
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_exc(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)
    elif vary == 'A':
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_a(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)

    try:
      table.sort('Epsilon_orig', reverse=True)
    except TypeError:
      table.sort('Epsilon_orig')
      table = table[::-1]

    if wavelen != {}:
        table = table[wavelen[0] < table['Lambda']]
        table = table[table['Lambda'] < wavelen[1]]
    if e_signif != {}:
        table = table[table['Epsilon_orig'] >= e_signif]

    if makefiles == True:
        fname = 'new epsilons for ' + element + str(z1) + ' '
        number = get_file_num(fname, '.csv')
        table.write(fname + str(number) + '.csv', format='csv')
        print("Wrote to", fname + str(number) + '.csv')

    print(table)
    return table

def get_peak_abund(Z, delta_ionize, delta_recomb, z1=[]):
    """ Returns dictionary of min, median, max peak abundance values from Monte Carlo CSD"""
    avg, low, high = monte_carlo_csd(Z, delta_ionize, delta_recomb)
    median, min, max = numpy.zeros([Z+1]), numpy.zeros([Z+1]), numpy.zeros([Z+1])
    for tmp_z1 in range(1, Z+1):
        peak_idx = numpy.argmax(avg[:, tmp_z1-1])
        median[tmp_z1-1], min[tmp_z1-1], max[tmp_z1-1] = avg[peak_idx, tmp_z1-1], low[peak_idx, tmp_z1-1], high[peak_idx, tmp_z1-1]
    if z1 != []:
        ret = {'median': median[z1-1], 'min': min[z1-1], 'max': max[z1-1]}
    else:
        ret = {'median': median, 'min': min, 'max': max}
    return ret

def find_max_error_csd(Z, z1, delta_ionize, delta_recomb):     ### in progress
    """ Finds temperature where CSD range from +/- delta_r is the largest"""
    median, min, max = monte_carlo_csd(Z, delta_ionize, delta_recomb)
    Telist = numpy.logspace(4,9,1251)
    percent_error, i = 0,0

    for x,y,z in zip(median[:, z1-1], min[:, z1-1], max[:, z1-1]):
        while (z-y)/x >= percent_error:
            percent_error: (z-y)/x
            i+=1

    ret = {'avg abund': median[i, z1-1], 'temp': '%.3E' % Decimal(Telist[i]), 'frac error': percent_error}
    return ret

def find_delta_temp(Z, frac, delta_r, vary='ir', z1={}, unit='K'):
    """Unit can be K or log. Returns dict of delta_Te values.
    and delta_Te values. Frac is fraction of peak abundance. """
    if isinstance(frac, int) or isinstance(frac, float): frac = [frac]
    Te_change = numpy.zeros([len(frac), Z + 1])
    new_temps = find_new_temp(Z, frac, delta_r, vary=vary, z1=z1, unit=unit)
    for i in range(len(frac)):
        for z1 in range(1, Z+1):
            for a,b,c in zip(ret['orig'][i,z1-1], ret['min'][i,z1-1], ret['max'][i,z1-1]):
                if unit.lower() == 'kev': Te_change[i, z1-1] = ((c-b)/a)/11604525.0061657
                elif unit.lower() == 'k': Te_change[i, z1-1] = (c-b)/a
    if z1 == {}: ret = {'unit': unit, 'frac':frac, 'delta_Te':Te_change[:, z1-1]}
    else: ret = {'unit': unit, 'frac':frac, 'delta_Te':Te_change}
    return ret

# def find_new_temp(Z, frac, delta_ionize={}, delta_recomb={}, vary='ir', z1={}, unit='K'):
#     """ Returns dict of orig/median, min, and max temp values at specified
#     fractional ionic concentration (as a fraction of peak). To vary 'i' or 'r'
#     rate, need to specify z1. To vary all rate coefficients, set vary
#     to 'ir.' Default is vary 'ir' with Monte Carlo CSD and
#     return dict of new value arrays indexed by z1-1. Frac can be int
#     (single frac) or list. Dict returned is of arrays with frac
#     number of rows, z1 columns (i.e. if z1 specified, returns array
#     of values for only z1, otherwise no z1 specified, gives new temps
#     for all ions).
#     """
#     if isinstance(frac, float) or isinstance(frac, int): frac = [frac]
#     orig_temps = numpy.zeros([len(frac), Z + 1])
#     min_temps, max_temps = numpy.zeros([len(frac), Z + 1]), numpy.zeros([len(frac), Z + 1])
#     Telist = numpy.logspace(4,9,51)
#     if (vary == 'i'):
#         eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, vary, delta_ionize)
#     elif (vary == 'r'):
#         eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, vary, delta_recomb)
#     else: eqpopn, pospopn, negpopn = monte_carlo_csd(Z, delta_ionize, delta_recomb, show=False)
#
#     for i in range(len(frac)):
#         for z1 in range(1, Z+1):
#             peak, index = numpy.max(eqpopn[:, z1 - 1]), numpy.argmax(eqpopn[:, z1 - 1])
#             orig_temps[i, z1-1] = numpy.interp(frac[i] * peak, eqpopn[:, z1-1][::-1], Telist[::-1])
#             min_temps[i, z1-1] = numpy.interp(frac[i] * peak, negpopn[:, z1 - 1][::-1], Telist[::-1])
#             max_temps[i, z1-1] = numpy.interp(frac[i] * peak, pospopn[:, z1 - 1][::-1], Telist[::-1])
#
#     if z1 == {}: ret = {'orig': orig_temps, 'min': min, 'max': max}
#     else: ret = {'orig': orig_temps[:, z1-1], 'min': min_temps[:, z1-1], 'max': max_temps[:, z1-1]}
#     return ret

def find_new_abund(Z, Te, delta_r, vary, unit='K', z1={}):
    """ Find orig/median and min and max values from new CSD.
    Te can be int (single temp) or list [] of multiple.
    Te unit is 'K' or 'keV'. Vary is either 'i' or 'r' (but
    must specify z1 ion to change rate for) or 'ir' to vary
    all rate coefficients for Z via Monte Carlo CSD.
    Regardless of vary, can specify z1 to return dict of only z1
    abundance arrays, otherwise returns dict of arrays
    for all ions. Arrays are length Te. Delta_r is error.
    Returns: dict = {'orig', 'min', max'} for vary = 'i' or 'r'
            dict = {'median', 'min', 'max'} for vary = 'ir'
    """

    if isinstance(Te, int): Te = [Te]
    if unit.lower() == 'keV': Te = [x*11604525.0061657 for x in Te]
    Telist = numpy.logspace(4,9,51)

    if (vary == 'i') or (vary == 'r'):  #CSD from only changing z1 rate
        orig, min, max = numpy.zeros(len(Telist), Z+1), numpy.zeros(len(Telist), Z+1), numpy.zeros(len(Telist), Z+1)
        if z1 == {}:
            print("Error: need to specify which z1 to change rate for. Exiting")
            exit()
        eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, vary, delta_r)
        for i in range(len(Te)):
            pop_fraction = pyatomdb.apec._solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
            for z1 in range(1, Z+1):
                orig[i, z1-1] = pop_fraction[z1 - 1]
                min[i, z1-1] = numpy.interp(Te[i], Telist, negpopn[:, z1 - 1])
                max[i, z1-1] = numpy.interp(Te[i], Telist, pospopn[:, z1 - 1])
        if z1 == {}:    #return full arrays
            ret = {'orig': orig, 'min': min, 'max': max}
        else:   #return only z1
            ret = {'orig': orig[z1-1], 'min': min[z1-1], 'max': max[z1-1]}

    elif vary == 'ir':  #Monte Carlo CSD
        median, min, max = numpy.zeros(len(Telist), Z+1), numpy.zeros(len(Telist), Z+1), numpy.zeros(len(Telist), Z+1)
        avg, low, high = monte_carlo_csd(Z, delta_r, delta_r)
        for i in range(len(Te)):
            for z1 in range(1, Z+1):
                median[i, z1 - 1] = numpy.interp(Te[i], Telist, avg[:, z1-1])
                min[i, z1 - 1] = numpy.interp(Te[i], Telist, low[:, z1-1])
                max[i, z1 - 1] = numpy.interp(Te[i], Telist, high[:, z1-1])
        if z1 == {}:    #return full arrays
            ret = {'orig': median, 'min': min, 'max': max}
        else:   #return only z1
            ret = {'orig': median[z1-1], 'min': min[z1-1], 'max': max[z1-1]}

    return ret

def find_abund_change(Z, Te, delta_r, vary, z1={}, unit='K'):
    """ Finds fractional abundance change from varying either 'i' or 'r'
    rate for specified z1 or vary all 'ir' rate coefficients for Z.
    Te can be int (single temp) or list. Delta_r is error.
    Te unit can be K or keV."""

    if isinstance(Te, int): Te = [Te]
    if unit.lower() == 'kev': [x * 11604525.0061657 for x in Te] #get temps in K
    abund_change = numpy.zeros([len(Te),Z+1])

    for i in range(len(Te)): #for each specified temp, find change
        abund = find_new_abund(Z, Te[i], delta_r, vary)
        for z1 in range(1, Z+1):
            abund_change[i, z1-1] = (abund['max']-abund['min'])/abund['orig']

    if z1 == {}: return abund_change
    else: return abund_change[:, z1-1]

def temps_of_interest(Z, z1, Telist=numpy.logspace(4,9,1251)):
    """Returns dict of temperatures in K and labels of where ion is at
    1%, 10%, peak fractional abundance, and back down to 10% and 1%."""

    eqpopn, fracs = get_orig_popn(Z, Telist), [1e-2, 0.1]
    peak = numpy.argmax(eqpopn[:, z1 - 1])
    increasing_temps = numpy.interp(fracs, eqpopn[:peak, z1 - 1], Telist[:peak])
    decreasing_temps = numpy.interp(fracs, eqpopn[peak:, z1 - 1][::-1], Telist[peak:][::-1])
    temp_bins = numpy.array([increasing_temps[0], increasing_temps[1], Telist[peak], decreasing_temps[1], decreasing_temps[0]])
    return temp_bins

def line_sensitivity(Z, z1, up, lo, vary, errors, trans_list, temps={}, pop_fraction={}, dens=1, trans_labels={}, show=True, savefig=True):
    """ Calculates fractional change in emissivity dE/E for
    specified Z (element), z1 (ion charge+1), up, lo line.
    Vary either 'exc' or 'A' for each transition in trans_list
    by each error in errors. Default density is 1 and default
    temps are temps where Z, z1 ionic fraction is 1%, 10%, peak,
    back down to 10% and 1%. Temps can be 'peak' for just the peak
    temperature or a single value.

    Z : int
        element
    z1 : int
        ion charge +1
    up : int
        upper level of transition to look at
    lo : int
        lower level of transition to look at
    vary : str
        type of rate to vary: 'exc' or 'A'
    errors : list or tuple
        list of errors on rate
    trans_list : list or tuple
        list of (up, lo) for transitions to vary
    temps : list or tuple
        temperatures to calculate new emissivity for
        default is 5 temps of interest
    pop_fraction : array
        perturbed pop fraction i.e. array outputted
        from monte_carlo_csd, otherwise use unperturbed
    dens : float
        density, default is 1
    trans_labels : list or tuple
        list of legend labels for each transition
        (should be same length as trans_list)
    show : bool
        True to show to screen
    savefig : bool
        save figure (False for wrapper routines using this)
    """

    if plot_settings['use_latex_font'] == True:
        try:
            setup_text_plots(usetex=True)
        except:
            setup_text_plots(usetex=False)
            print("Using regular matplotlib font for plotting, could not find an installation of LaTeX")

    print("Pop fraction is:", pop_fraction)
    if temps == 'peak': temps = [find_peak_Te]
    if isinstance(temps, int) or isinstance(temps, float): temps = [temps]
    elif temps == {}: temps = temps_of_interest(Z, z1, numpy.logspace(4,9,1251))

    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)

    if trans_list == {}:
        if (Z, z1) == (26,17):
            trans_list = [(2,1), (3,1), (5,1), (17, 1), (23,1), (27,1)]
        if (Z, z1) == (26, 19):
            trans_list = [(53,1), (68,1), (71,1), (74, 1), (76, 4)]
        if Z-z1 == 1:
            trans_list = [(2,1), (5,1), (6,1), (7,1), (13,1)]
            print("Varying He-like lines.")
        if Z==z1:
            print("Varying H-like lines.")
            trans_list=[(3,1), (4,1), (6,1), (7,1)]

    l = str(up) + '->' + str(lo)
    plt.figure(figsize=(10,8))
    plt.xlabel('Fractional error')
    plt.ylabel('Fractional error in emissivity')
    plt.title(element + ' ' + ion + ' ' + l + ' sensitivity')

    matrix, legend_labels = numpy.zeros([len(trans_list), len(temps), len(errors)]), []
    clist = get_cmap(len(temps) + 1, name=plot_settings['cmap'])    #color is for # of temps
    markers = ['o', 'v', 's', 'P', '^', '2']    #marker is for each transition
    temp_str = ['%.1E' % Decimal(x) for x in temps]

    #vary each transition in trans_list separately
    for a, transition in enumerate(trans_list):
        up, lo = transition[0], transition[1]
        
        #add legend label
        if trans_labels == {}:
            legend_labels.append(Line2D([0], [0], marker=markers[a], label=str(transition[0]) + '->' + str(transition[1])))
        else:
            legend_labels.append(Line2D([0], [0], marker=markers[a], label=trans_labels[a]))
            
        #loop through each temperature
        for b, Te in enumerate(temps):
            legend_labels.append(Line2D([0], [0], color=clist(b), label=temp_str[b]))
            
            #for each temperature, vary rate by the specified errors
            for c, delta_r in enumerate(errors):
                print("Varying rate by", delta_r*100, "%")
                extras = {'process': vary, 'delta_r': delta_r, 'transition': [], 'transition_2': [],
                       'wavelen': (10, 20), 'Te_range': {}, 'dens_range': {}, 'corrthresh': 0.0, 'e_signif': 0.0, 'pop_fraction':pop_fraction}
                if vary == 'exc':
                    extras['transition'] = (lo, up)
                    inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                    new_inputs, new_values = vary_exc(inputs, values, transition)
                elif vary == 'A':
                    extras['transition'] = (up, lo)
                    inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                    new_inputs, new_values = vary_a(inputs, values, transition)
                table, new_table, inputs, results = get_tables(new_inputs, new_values)

                i = numpy.where([(new_table['Upper'] == up) & (new_table['Lower'] == lo)])[1]
                if new_table[i]['|dE/E|'] < 0.00001: matrix[a,b,c] = 0
                else: matrix[a,b,c] = new_table[i]['|dE/E|']
            print(len(temps)-b-1, "temps left")
            plt.plot(errors, matrix[a,b,:], color=clist(b), linestyle='--', marker=markers[a])

    plt.legend(handles=legend_labels, fontsize='xx-small', loc='upper left')
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)
    plt.savefig(element+' '+ion+' '+str(up)+'-'+str(lo)+' sensitivity.pdf')
    if show == True:
        plt.show()
    return matrix

def get_orig_popn(Z, Telist):
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for temp_z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=temp_z1, extrap=True, settings=False, datacache=d)
        ionlist[:, temp_z1 - 1] = iontmp
        reclist[:, temp_z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)
    return eqpopn

def find_csd_change(Z, z1, delta_ionize, delta_recomb={}, Te={}, frac={}, varyir={}, printout=True):
    """ Find the change in CSD at specified temp or ion fraction.
    Default is applying random error from Gaussian distribution
    to all rates with a maximum error/sigma = delta_r.
    If varyir= 'i' or 'r', will gives changes for varying i or r
    rate of only the z1 ion and will use delta_ionize as error on rate.
    Te in K. Frac is abundance as a fraction of peak, i.e. half peak is 0.5.
    z1 can either be an int (single ion) or a list [] of multiple.
    If you don't want the routine to print out the orig and new values,
    set printout=False. If error is not specific to type of rate coeff
    (i.e. same error for ionization and recombination), leave delta_recomb empty.

    Returns: list of percent_change_Te and/or percent_change_abund
    where percent change is the difference from +/- error as a fraction of original value."""

    percent_change_Te, percent_change_abund = [], []
    abund_errors, Te_errors = [], []

    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)

    if type(z1) == int:
        ions = [z1]
    else:
        ions = z1

    if varyir != {}:
        #vary z1 rate only
        for z1 in ions:
            if varyir == 'i':
                print("Varying ionization rate for", element, str(z1-1), "+ ion only...\n")
            elif varyir == 'r':
                print("Varying recombination rate for", element, str(z1 - 1), "+ ion only...\n")

            eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, varyir, delta_ionize)

            #if frac abundance specified, find temp change in new CSD at that fraction
            if frac != {}:
                peak = numpy.max(eqpopn[:, z1 - 1])
                index = numpy.argmax(eqpopn[:, z1 - 1])

                # interpolate
                min_temp = numpy.interp(frac * peak, negpopn[index:, z1-1][::-1], Telist[index:][::-1])
                max_temp = numpy.interp(frac * peak, pospopn[index:, z1-1][::-1], Telist[index:][::-1])

                orig_dt = max_temp - min_temp
                orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)

                if printout == True:
                    print(element, str(z1-1), "+ original temperature at", frac, "*peak is", Telist[index], "(Log(T) = ", numpy.log10(Telist[index]))
                    print("\n New temperature range in K is", min_temp, "to", max_temp, "i.e. dT =", orig_dt)
                    print("\n Log(T) range is", numpy.log10(min_temp), "to", numpy.log10(max_temp), "i.e. dLog(T) =", orig_dex_change)
                    print("\n")
                percent_change_Te.append(abs(orig_dt/Telist[index])/2)

            #if temperature given, find abundance change in new CSD at that temperature
            if Te != {}:
                #original abundance
                pop_fraction = pyatomdb.apec._solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
                orig = pop_fraction[z1-1]

                min = numpy.interp(Te, Telist, negpopn[:, z1 - 1])
                max = numpy.interp(Te, Telist, pospopn[:, z1 - 1])

                #find change
                dA = max - min
                log_dA = -numpy.log10(max) - (-numpy.log10(min))

                if printout == True:
                    print(element, str(z1 - 1), "+ original abundance at", "%.2E" % Decimal(Te), "K =", orig, "(-Log10(X) =", -numpy.log10(orig), ')\n')
                    print("New abundance range is", min, "to", max, "i.e. d_abund =", dA)
                    print("\n(-Log10(X) new range is", -numpy.log10(min), "to", -numpy.log10(max), '), i.e. d_log_abund =', log_dA)
                    print("\n")
                percent_change_abund.append(abs(dA/orig)/2)
                #errors_dicts.append({'+': max-orig, '-': orig-min})

    #if varyir not specified, monte carlo rates
    else:
        print("\n Varying all rate coefficients for", element, '\n')
        eqpopn, pospopn, negpopn = monte_carlo_csd(Z, delta_ionize, delta_recomb)

        for z1 in ions:
            # if frac abundance specified, find temp change in new CSD at that fraction
            if frac != {}:
                peak = numpy.max(eqpopn[:, z1 - 1])
                index = numpy.argmax(eqpopn[:, z1 - 1])

                # interpolate
                min_temp = numpy.interp(frac * peak, negpopn[index:, z1 - 1][::-1], Telist[index:][::-1])
                max_temp = numpy.interp(frac * peak, pospopn[index:, z1 - 1][::-1], Telist[index:][::-1])

                orig_dt = max_temp - min_temp
                orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)

                if printout == True:
                    print(element, str(z1 - 1), "+ original temperature at", frac, "*peak is", Telist[index], "(Log(T) = ",
                    numpy.log10(Telist[index]))
                    print("\n New temperature range in K is", min_temp, "to", max_temp, "i.e. dT =", orig_dt)
                    print("\n Log(T) range is", numpy.log10(min_temp), "to", numpy.log10(max_temp), "i.e. dLog(T) =",
                          orig_dex_change)
                    print("\n")
                percent_change_Te.append(abs(orig_dt / Telist[index])/2)
                #errors_dicts.append({'+': max - orig, '-': orig - min})

            # if temperature given, find abundance change in new CSD at that temperature
            if Te != {}:
                print("z1 =", z1, "and z1+1 =", z1+1)
                orig = numpy.interp(Te, Telist, eqpopn[:, z1-1])
                min = numpy.interp(Te, Telist, negpopn[:, z1 - 1])
                max = numpy.interp(Te, Telist, pospopn[:, z1 - 1])
                print("orig:", orig, "min:", min, "max:", max)

                # find change
                dA = max - min
                log_dA = -numpy.log10(max) - (-numpy.log10(min))
                percent_error = (abs(dA)/orig)/2

                if printout == True:
                    print('\n', element, str(z1 - 1), "+ original abundance at", "%.2E" % Decimal(Te), "K =", orig,
                          "(-Log10(X) =", -numpy.log10(orig), ')\n')
                    print("New abundance range is", min, "to", max, "i.e. d_abund =", dA, "for orig value =", orig)
                    print("\n(-Log10(X) new range is", -numpy.log10(min), "to", -numpy.log10(max), '), i.e. d_log_abund =',
                          log_dA)
                    print("So abundance =", orig, "+", max-orig, "and -", orig-min)
                    print("\n")
                #percent_change_abund.append(abs(dA / orig))
                percent_change_abund.append(percent_error)
                abund_errors.append(percent_error)
                #abund_error_dicts.append({'+': max - orig, '-': orig - min})

    if (Te != {}) & (frac == {}):
        print("Fractional change in ionic concentration:", percent_change_abund)
        return abund_errors
    elif (Te == {}) & (frac != {}):
        print("Fractional change in temperatures:", percent_change_Te)
        return percent_change_Te
    else:
        print("Fractional change in ionic concentrations:", percent_change_abund, "\nFractional change in temperatures:", percent_change_Te)
        return percent_change_Te, percent_change_abund

def error_analysis(Z, z1, up, lo, Te, dens, errors, linewidth=2.5, filename={}, use_max_csd=True, cascade_error=0.15, direct_exc_error=0.15):    #in progress
    """ Generates error analysis PDF of CSD and emissivity errors
    for Z, z1 and transition (up, lo) at Te in K. Linewidth is in eV.
    Errors is either a single fractional error, i.e. 0.1 for 10%,
    or a list [] of 2 fractional error, [ionization, recombination].
    If one error specified, uses for both ionization and recombination.
    If filename for output pdf has spaces, they will be replaced
    with underscores (_), i.e. 'O 7 analysis' makes O_7_analysis.pdf
    Default file generated is ErrorAnalysis.pdf"""

    print("Linewidth is", linewidth, "eV")
    if isinstance(errors, (tuple, list)): delta_ionize, delta_recomb = errors[0], errors[1]
    else:
        delta_ionize = errors
        delta_recomb = errors
        print("Using", delta_ionize, "as ionize error and", delta_recomb, "as recomb error")
    print("Using", direct_exc_error, "as error on excitation rate from 1->"+str(up), "and", cascade_error, "as error on cascades")
    if filename == {}: filename = 'ErrorAnalysis'
    if '' in filename:
        filename = filename.split()
        filename = '_'.join(filename)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    name = element + ' ' + ion

    # get LV data
    lvdata = pyatomdb.atomdb.get_data(Z, z1, 'LV', datacache=d)
    lvdata = lvdata[1].data
    up_config, lo_config = lvdata['ELEC_CONFIG'][up - 1], lvdata['ELEC_CONFIG'][lo - 1]

    # get emissivity and multiply by elemental abundance
    ab = pyatomdb.atomdb.get_abundance()[Z]  # abundance of element Z
    inputs, values = set_up(Z, z1, Te, dens)
    table = values.get('table')
    i = numpy.where([(table['Upper'] == up) & (table['Lower'] == lo)])
    emiss = table[i]['Epsilon_orig']*ab

    # get A values
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
    in_range = ladat[1].data
    i = numpy.where([(in_range['UPPER_LEV'] == up) & (in_range['LOWER_LEV'] == lo)])
    lambda_obs, lambda_theory = in_range[i]['WAVE_OBS'], in_range[i]['WAVELEN']
    print("Lambda obs", lambda_obs, "lambda theory", lambda_theory)
    ref_obs, ref_theory = in_range[i]['WV_OBS_REF'], in_range[i]['WAVE_REF']
    A_val, A_ref = in_range[i]['EINSTEIN_A'], in_range[i]['EIN_A_REF']

    # check if references are bibcodes or just text
    citations, refs = {'obs': '', 'theory': '', 'A': ''}, [ref_obs, ref_theory, A_ref]
    for x, y in zip(refs, list(citations.keys())):
        if " " not in x:
            citations[y] = '\\href{https://ui.adsabs.harvard.edu/abs/' + x + '/abstract}{' + x + '})'
        else:
            citations[y] = x

    # find temp of peak line emissivity at low dens
    kT = Te / 11604525.0061657
    Telist = numpy.linspace(Te/10, Te*10, 51)
    s = pyatomdb.spectrum.CIESession()
    ldata = s.return_line_emissivity(Telist, Z, z1, up, lo, teunit='K')
    peak_emiss, peak_Te = 0, 0
    for i in range(len(Telist)):
        while ldata['epsilon'][i] > peak_emiss:
            peak_emiss = ldata['epsilon'][i]
            peak_Te = Telist[i]

    #now get emiss at peak_Te at specified density
    inputs, values = set_up(Z, z1, peak_Te, dens)
    table = values.get('table')
    i = numpy.where([(table['Upper'] == up) & (table['Lower'] == lo)])
    peak_emiss = table[i]['Epsilon_orig'] * ab

    #compare line emiss to peak emiss
    if emiss < peak_emiss: compare_emiss = 'below'
    elif emiss > peak_emiss: compare_emiss = 'above'
    elif emiss == peak_emiss: compare_emiss = 'at'

    ######## get emissivity contributions #######
    exc, cascade, cascade_dict, other_cascades, recombination, ionization, all_cascades = get_contributions(Z, z1, up, lo, Te, dens)

    # CSD uncertainty
    pop_fraction, ions = pyatomdb.apec._solve_ionbal_eigen(Z, Te, teunit='K'), []
    for i in range(1, Z+2):
        if pop_fraction[i-1] >= 0.001:     #get ions with abund > 0.1%
            ions.append(i)

    print("Finding CSD errors for ions:", ions)
    str_ions = [element + '$^{+' + str(x - 1) + '}$' for x in ions]
    CSD = [round(pop_fraction[x-1]*100,2) for x in ions]
    CSD_percent_errors= find_csd_change(Z, ions, delta_ionize, delta_recomb, Te=Te, printout=False)  # MC
    CSD_errors = [int(round(x*y,2)) for x,y in zip(CSD_percent_errors, CSD)]

    # get energies and flux uncertainty
    ionization_energy = pyatomdb.atomdb.get_ionpot(Z, z1, datacache=d)
    ionization_energy /= 1000  # in keV

    if numpy.isnan(lambda_obs):
        excitation_energy = (12.398 / lambda_theory) * 1000
        print("\nCalculating excitation energy with theoretical wavelength.")
    else:
        excitation_energy = (12.398 / lambda_obs) * 1000  # in eV
        print("\nCalculating excitation energy with observational wavelength.")

    if kT < excitation_energy: compare_exc = 'below'
    elif kT > excitation_energy: compare_exc = 'above'
    elif kT == excitation_energy: compare_exc = 'at'

    if kT < ionization_energy: compare_ionize = 'below'
    elif kT > ionization_energy: compare_ionize = 'above'
    elif kT == ionization_energy: compare_ionize = 'at'

    # check for DR sats within 2.5 eV with epsilon > 1e-20
    DR_list, sat_blend_flux = {}, 0
    a, b = pyatomdb.apec.calc_satellite(Z, z1 - 1, Te)
    if numpy.isnan(lambda_obs) == False:
        en = 12.398425 / lambda_obs
    else:
        en = 12.398425 / lambda_theory
    emin, emax = en - (linewidth/1000), en + (linewidth/1000)
    wvmin, wvmax = 12.398425 / emax, 12.398425 / emin
    a = a[(a['lambda'] >= wvmin) & (a['lambda'] <= wvmax)]  # filter for energy range
    a['epsilon'] *= ab * pop_fraction[z1 - 1]
    a = a[a['epsilon'] > 1e-20]
    sat_blend_flux = sum(a['epsilon'])
    new = sum(a['epsilon'])
    print("Comparing sum", sat_blend_flux, "and sum*pop fraction", new)
    DR_list.update({str(len(a)): element + ' ' + pyatomdb.atomic.int_to_roman(z1)})

    other_lines, total_blend_flux = {}, sat_blend_flux
    lines = s.return_linelist(kT,[wvmin,wvmax])

    #filter other lines for epsilon > 1% of transition's flux
    lines['Epsilon'] *= ab
    lines = lines[lines['Epsilon'] >= 0.01*emiss]
    print("ADDITIONAL LINES IN THIS LINEWIDTH:", lines)
    for i in range(len(lines[:7])):
        if (lines[i]['UpperLev'], lines[i]['LowerLev']) not in zip(a['upperlev'], a['lowerlev']):
            el, z_1 = pyatomdb.atomic.Ztoelsymb(int(lines[i]['Element'])), pyatomdb.atomic.int_to_roman(int(lines[i]['Ion']))
            total_blend_flux += lines[i]['Epsilon']
            if el+' '+z_1 in other_lines:
                num = other_lines.get(el+' '+z_1)
                other_lines.update({el+' '+z_1: num+1})
            else:
                other_lines.update({el+' '+z_1: 1})
    print(other_lines)

    #final error calculation            #still in progress
    in_range, matrix = values['in_range'], values['matrix']
    exc_init, exc_final, exc_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=False, do_ec=True, do_ir=False,
                                                                do_pc=False, do_ai=False, datacache=d)
    i = numpy.where([(exc_init == 0) & (exc_final == up - 1)])[1][0]
    new_emiss, edit_matrix = [], matrix.copy()
    for which in ['min', 'max']:
        if which == 'min': use_max_direct, use_max_a = False, False #minimize emissivity of (up, lo)
        else: use_max_direct, use_max_a = True, True    #maximize emissivity of (up, lo)
        new_matrix, original_vals, cascade_list = vary_a_cascades(Z, z1, up, lo, Te, dens, cascade_error, in_range, edit_matrix)
        values['matrix'] = new_matrix
        #vary exc rate from ground to upper level of transition
        old_rate = exc_rates[i]
        new_rate = [(1 - direct_exc_error) * old_rate, (1 + direct_exc_error) * old_rate]
        if use_max_direct == True: x = new_rate[-1]
        else: x = new_rate[0]
        new_matrix[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        new_matrix[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term
        new_matrix[0][:] = 1.0
        #recalculate line emissivity with perturbed CSD from ionization and recombination errors
        min, max = mc_popn(Z, delta_ionize, delta_recomb, Te, type='f')
        if use_max_csd == True:
            extras = {'process': {}, 'delta_r': {}, 'transition': (up, lo), 'transition_2': [], 'npnts': 2,
                     'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': 10e-5,
                      'e_signif': 0.0, 'pop_fraction': max}
        else:
            extras = extras = {'process': {}, 'delta_r': {}, 'transition': (up, lo), 'transition_2': [], 'npnts': 2,
                     'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': 10e-5,
                      'e_signif': 0.0, 'pop_fraction': min}
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        table = inputs['table']
        i = numpy.where([(table['Upper'] == up) & (table['Lower'] == lo)])
        new_emiss.append(table[i]['Epsilon_orig'] * ab)
        #restore a values from cascades
        restore_a_values(original_vals, cascade_list)

    max_new_emiss, min_new_emiss = new_emiss[-1], new_emiss[0]
    print("Min new", min_new_emiss, "max new", max_new_emiss, "compared to orig=", emiss)
    final_error = ((abs(max_new_emiss - min_new_emiss)/emiss)/2)

    # write to tex file
    with open(filename+'.tex', 'w') as file:
        file.write('\\documentclass[11pt]{article}\n\n')
        file.write('\\usepackage{hyperref}\n\n')
        file.write('\\begin{document}\n\n')
        file.write('\\begin{center}\n{\Large Uncertainty Analysis for ' + name + ' ' + str(up) + '$\\rightarrow$' + str(lo) + '\\\ \n')
        file.write('at ' + '%.2E' % Decimal(Te) + ' K and ' + '%.E' % Decimal(dens) + ' cm$^{-3}$ \\\ \n')
        file.write('\n\\large AtomDB v3.0.9; pyAtomDB error analysis v0.1}\n')
        file.write('\\end{center}\n\n')

        file.write('\\noindent \\begin{tabular}{lll}\n')
        file.write('\\hline Transition & ' + str(up_config) + ' $\\rightarrow$ ' + str(lo_config) + ' & \\\ \n')

        if numpy.isnan(lambda_obs) == False:
            file.write('$\\lambda_{\\rm obs}$ & ' + str(lambda_obs) + ' \\AA & ' + citations['obs'] + ' \\\ \n')

        file.write('$\\lambda_{\\rm th}$ & ' + str(lambda_theory) + ' \\AA & ' + citations['theory'] + ' \\\ \n')
        file.write('Einstein A & ' + '%.2E' % Decimal(float(A_val)) + ' s$^{-1}$ & ' + citations['A'] + ' \\\ \n')
        file.write('$\Lambda$(' + '%.2E' % Decimal(Te) + 'K) & ' + '%.3E' % Decimal(
            emiss) + ' ph cm$^3$ s$^{-1}$ & \\\ \\hline \n')
        file.write('\\end{tabular}\n\n')

        file.write('\\vskip 0.2in \n')
        file.write('\\noindent{\\bf Flux uncertainty} \n')
        file.write('\\vskip 0.1in \n')
        file.write(
            '\\noindent The temperature is ' + compare_emiss + ' the peak line emissivity, and the plasma temperature (' + str(
                round(kT, 3)) + ' keV) ')
        file.write('is ' + compare_exc + ' the excitation energy (' + str(round(excitation_energy / 1000, 2)) + ' keV) ')
        file.write('and ' + compare_ionize + ' the ' + element + '$^{+' + str(z1 - 1) + '}$ ionization energy (' + str(
            round(ionization_energy, 3)) + ' keV). ')
        file.write('Excitations will be dominated by hard-to-calculate resonances (see XXX). \n')

        file.write('\\vskip 0.1in \n')
        file.write('\\noindent Charge state distribution: \\\ \n')
        file.write('\\noindent \\begin{tabular}{ll}\n')
        for ion, frac, error in zip(str_ions, CSD, CSD_errors):
            file.write(ion + ' & ' + str(frac) + '\\% $\\pm$ ' + str(error) + '\\% \\\ \n')
        file.write('\\end{tabular}\n')
        file.write('\\vskip 0.1in \n')

        file.write('\\noindent Detailed breakdown: \\\ \n')
        file.write('\\noindent \\begin{tabular}{lll} \n')
        file.write('\\hline Direct Excitation & & ' + str(round(exc * 100, 1)) + '\\% \\\ \n')

        for key, val in cascade_dict.items():
            file.write('Cascade & ' + key + ' & ' + val + '\\\ \n')

        file.write('Other cascades & & ' + str(round(other_cascades*100, 1)) + '\\% \\\ \n')
        file.write('Ionization & & ' + str(round(ionization * 100, 1)) + '\\% \\\ \n')
        file.write('Recombination & & ' + str(round(recombination * 100, 1)) + '\\% \\\ \\hline \n')
        file.write('\\end{tabular} \n')
        file.write('\\vskip 0.1in \n')

        if delta_recomb != {}:
            file.write('Using ionization $\\Delta$R ' + str(delta_ionize*100) + ' and recombination $\\Delta$R '
            + str(delta_recomb*100) + '\\%, estimated total error is ' + str(round(final_error * 100)) + '\\% \\\ \n')
        else:
            file.write('Using max rate uncertainty ' + str(delta_ionize*100) + '\\%, estimated total error is ' + str(
                round(final_error * 100, 1)) + '\\% \\\ \n')
        file.write('\\indent ** {\\it Error is approximated by summing individual errors in quadrature.} \\\ \n')
        file.write('\\vskip 0.1in\n')

        if len(DR_list) != 0:
            file.write('\\noindent{\\bf Line blends within $\\pm$ ' + str(linewidth) + 'eV} \n')
            file.write('\\vskip 0.1in \n')
            for key, val in DR_list.items():
                file.write('\\indent ' + key + ' DR satellites: ' + val + '\\\ \n')
            file.write(
                '\\indent Satellite blend flux: ' + str(round((sat_blend_flux / emiss) * 100, 1)) + '\\% (' + '%.2E' % Decimal(
                    sat_blend_flux) + ' ph cm$^3$ s$^{-1}$) \\\ \n')
        else:
            file.write('\\noindent{\\bf No DR satellites within $\\pm$ ' + linewidth + 'eV} \n')

        if len(other_lines) > 0:
            file.write('\\indent Other lines: \\\ \n')
            for key, val in other_lines.items():
                file.write('\\indent ' + str(val) + ' lines: ' + key + '\\\ \n')
            file.write('\\indent Total blended flux: ' + str(round((total_blend_flux/emiss)*100,1)) + '\\% (' + '%.2E' % Decimal(
                        total_blend_flux) + ' ph cm$^3$ s$^{-1}$) \\\ \n')

        if len(DR_list) != 0:
            file.write('\\indent Note: DR lines have highly uncertain fluxes in general; 300\\% in some cases. ' +
                       '$\lambda$\ also uncertain.  Need to compare to Badnell for totals. \\\ \n')

        file.write('\\end{document}')

    os.system("pdflatex " + filename+'.tex')
    #remove extra files
    for ext in ['.log', '.out', '.aux', '.tex']:
        os.remove(filename+ext)

def ion_sensitivity(Z, errors, z1_plot={}, z1_interesting={}, distribution='f', show_legend=True, show=True):
    """ Finds fractional change in ion peak abundances from
    varying CSD with random errors from Gaussian or flat distribution
    using each delta_r from errors as sigma.
    If z1 empty, calculates for all ions."""

    Telist, element,  = numpy.logspace(4, 9, 1251), pyatomdb.atomic.Ztoelsymb(Z)
    if z1_plot == {}:
        z1_list = [z1 for z1 in range(1,Z+2)]
    elif z1_plot != {} and isinstance(z1_plot, int): z1_list = [z1_plot]
    else: z1_list = list(z1_plot)

    ions = [element + ' +'] + [element + ' ' + str(z1 - 1) + '+' for z1 in range(2, Z + 2)]
    values = Table()
    values['Ion'] = ions

    #plot peak concentration change for all z1 specified
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Fractional change on rate coefficients')
    ax1.set_ylabel('Fractional change in peak ionic concentration')
    clist = get_cmap(Z + 2)

    #if z1_interesting, make second plot
    if z1_interesting != {}:
        fig2, ax2 = plt.subplots()
        ax2.set_ylabel('Fractional change in ionic concentration')
        ax2.set_xlabel('Temperature (K)')
        clist2 = get_cmap(len(errors) + 2)

    for i, max_error in enumerate(errors):
        median, min, max = monte_carlo_csd(Z, max_error, max_error, runs=1000, makefiles=False, show=False, Te_range=[1e4,1e9])
        percent, numbers_only = [], []
        for z1 in range(1, Z+2):
            #find ion abundances at peak abundance
            peak = numpy.max(median[:, z1 - 1])
            peak_idx = numpy.argmax(median[:, z1 - 1])
            min_p = min[peak_idx, z1 - 1]
            max_p = max[peak_idx, z1 - 1]
            val = abs(((max_p - min_p) / peak)/2)
            if val < 0.01: val = 0
            numbers_only.append(round(val,2))
            val = numpy.format_float_positional(val*100, precision=2)
            percent.append(val+' %')

            if z1_interesting != {} and z1 == z1_interesting:
                idx = numpy.where(median[:, z1 - 1] >= 10e-3)[0]
                temps, median2, min2, max2 = Telist[idx], median[idx, :], min[idx, :], max[idx, :]
                plot_errors = ((max2[:, z1 - 1] - min2[:, z1 - 1]) / median2[:, z1 - 1]) / 2
                ax2.semilogx(temps, plot_errors, color=clist2(i), label="{}%".format(int(max_error * 100)))

        values['+/- '+str(round(max_error*100,2))+'%'] = numbers_only

    if z1_interesting != {}:
        ax2ticks = ax2.get_yticks()
        ax2.set_yticks(ax2ticks, numpy.round(ax2ticks, 2)[0])
        peak_Te = find_peak_Te(Z, z1_interesting)
        ax2.axvline(peak_Te, linestyle='--', color='grey', alpha=0.4)

    ticks = numpy.arange(0, numpy.max(numbers_only) + 0.05, 0.05)
    ax1.set_yticks(ticks)
    i = numpy.where(numpy.array(errors) < numpy.max(numbers_only) + 0.05)[0]
    ax1.plot(numpy.array(errors)[i], numpy.array(errors)[i], color='k', linestyle='--')

    for z1 in z1_list:
        v = [values[col][z1-1] for col in values.colnames[1:]]
        ax1.plot(errors, v, label=ions[z1-1], color=clist(z1-1))

    if show_legend == True:
        ax1.legend(fontsize='xx-small', loc=plot_settings['legendloc'])
        if z1_interesting != {}: ax2.legend(fontsize='xx-small', loc=plot_settings['legendloc'])

    fig.savefig(element+' peak ion concentration change.pdf')
    if z1_interesting != {}:
        fig2.savefig(element+' '+str(z1_interesting-1)+'+ sensitivity.pdf')

    if show == True:
        plt.show()
        plt.close('all')

def get_levels(Z, z1, wavelengths, relative_diff=10e-5):
    """ Returns list of levels (up, lo) used by ATOMDB
    for specified wavelengths. Will compare to observed
    wavelength if available, theoretical otherwise.
    Relative difference specifies the delta_lambda amount
    acceptable, default is 10e-5.

    Z : int
        element
    z1 : int
        ion charge +1
    wavelengths : list or tuple
        list of wavelengths in A
    relative_diff : float
        maximum difference in lambda for matching"""

    inputs, values = set_up(Z, z1, 3e6, 1)
    table = values.get('table')
    levels = []
    for wave in wavelengths:
        i = numpy.where(math.isclose(table['Lambda'], wave, rel_tol=relative_diff) == True)
        print("Lambda =", wave, "is:", table[i]['Upper'], table[i]['Lower'])
        levels.append((table[i]['Upper'], table[i]['Lower']))
    return levels

def get_wavelengths(Z, z1, transitions):
    """ Returns list of wavelengths in A
    for specified transitions.

    Z : int
        element
    z1 : int
        ion charge +1
    transitions : list or tuple
        list of transitions (up, lo)"""

    inputs, values = set_up(Z, z1, 3e6, 1)
    table = values.get('table')
    waves = []
    for line in transitions:
        i = numpy.where((table['Upper'] == line[0]) & (table['Lower'] == line[1]))
        print("Transition", line[0], "->", line[1], "has lambda =", table[i]['Lambda'])
        waves.append(table[i]['Lambda'])
    return waves

def get_strong_lines(Z, z1, Te, dens, cutoff=1e-13, wavelen=(1,30)):
    """ Returns list of (up, lo) strong lines for the ion
    at specific Te and dens given an emissivity cut off
    (with abundance correction factored in)and wavelength range.
    Default cut off is epsilon = 1e-13
    and default wavelength range is 1-30 A."""
    strong_lines = []
    inputs, values = set_up(Z, z1, Te, dens)
    table = values.get('table')
    ab = pyatomdb.atomdb.get_abundance()[Z]
    i = numpy.where((table['Epsilon_orig']*ab >= cutoff) & (table['Lambda'] in range(wavelen[0], wavlen[1])))
    strong_lines = [(int(row['Upper']), int(row['Lower'])) for row in table[i]]
    return strong_lines

def get_new_popns(Telist, Z, z1_test, varyir, delta_r):
    factor = delta_r
    ionlist, reclist = numpy.zeros([len(Telist), Z]), numpy.zeros([len(Telist), Z])

    # get the rates
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True, datacache=d)
        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

    # copy this rates
    iontmp, rectmp = ionlist * 1.0, reclist * 1.0

    # multiply rates by + factor
    if varyir.lower() == 'r':
        rectmp[:, z1_test - 1] *= (1 + factor)
    elif varyir.lower() == 'i':
        iontmp[:, z1_test - 1] *= (1 + factor)
    pospopn = solve_ionrec(Telist, iontmp, rectmp, Z)

    # here I am introducing a division to get you back to the original value before we start again
    if varyir.lower() == 'r':
        rectmp[:, z1_test - 1] /= (1 + factor)
        rectmp[:, z1_test - 1] *= (1 - factor)
    elif varyir.lower() == 'i':
        iontmp[:, z1_test - 1] /= (1 + factor)
        iontmp[:, z1_test - 1] *= (1 - factor)
    negpopn = solve_ionrec(Telist, iontmp, rectmp, Z)

    return eqpopn, pospopn, negpopn

def ion_ratio(Z, z1, z2, max_ionize, max_recomb={}, Te_range=[1e4,1e9]):

    if plot_settings['latex_font'] == True:
        try:
            setup_text_plots(usetex=True)
        except:
            setup_text_plots(usetex=False)
            print("Using regular matplotlib font for plotting, could not find an installation of LaTeX")

    if max_recomb == {}: max_recomb = max_ionize
    median, min, max = monte_carlo_csd(Z, max_ionize, max_recomb, Te_range=Te_range, runs=500)
    element, ion1, ion2 = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1), pyatomdb.atomic.int_to_roman(z2)
    name = element + ' ' + ion1 + '/' + ion2

    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    ax2.set_xlabel('Temperature (K)', fontsize=plot_settings['fontsize'])
    ax1.set_ylabel(name)
    ax2.set_ylabel('Fractional Error', fontsize=plot_settings['fontsize'])
    for axis in [ax1, ax2]:
        for which in ['x', 'y']:
            axis.tick_params(axis=which, labelsize=plot_settings['ticksize'])
    Telist = numpy.logspace(numpy.log10(Te_range[0]), numpy.log10(Te_range[1]), 1251)

    real, low, high, temp_bins = [], [], [], []
    for x, y, z, Te in zip(max, min, median, Telist):
        if z[z1-1] and z[z2-1] >= 10e-2:
            real.append(z[z1-1]/z[z2-1])
            low.append(y[z1-1]/x[z2-1])
            high.append(x[z1-1]/y[z2-1])
            temp_bins.append(Te)
    ax1.semilogx(temp_bins, real, color='b')
    ax1.fill_between(temp_bins, low, high, color='b', alpha=0.5)

    error = [abs(a-b)/z for a,b,z in zip(high, low, real)]
    ax2.semilogx(temp_bins, error, color='b')
    ax2.set_yticks(numpy.arange(0, numpy.max(error), 4))
    ax2.set_yticklabels([str(x) for x in numpy.arange(0, numpy.max(error), 4)])
    plt.subplots_adjust(hspace=0)
    fig.savefig(element + ' ' + ion1 + '-' + ion2 +' ratio.pdf')
    plt.show()
    plt.close()

def get_level_pop(Z, z1, Te, dens, all=False):
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, datacache=d)
    levpop = pyatomdb.apec.solve_level_pop(init, final, rates, False)

    if all:
        return init, final, rates, levpop
    else:
        return levpop

def vary_level_pop(init, final, rates, up, lo, new_rate):
    #vary rate of up, lo transition
    for i in range(len(rates)):
        if (init[i],final[i]) == (up, lo): rates[i] = new_rate

    #recalculate level populations
    levpop = pyatomdb.apec.solve_level_pop(init, final, rates, False)
    return levpop

def get_contributions(Z, z1, up, lo, Te, dens):

    # Te is in Kelvin, kT is in keV
    kT = Te*pyatomdb.const.KBOLTZ

    ## Step 0: get the ion balance and element abundance

    abund = pyatomdb.atomdb.get_abundance()[Z]
    ion_frac = pyatomdb.apec.return_ionbal(Z,Te, datacache=d)

    ## Step 1 gather all the rates, solve for level population ##

    init, final, rates, levpop = get_level_pop(Z, z1, Te, dens, all=True)

    levpop*=ion_frac[z1-1]*abund

    # level populations of adjacent ions
    levpop_m1 = get_level_pop(Z, z1-1, Te, dens)*ion_frac[z1-2]*abund
    levpop_p1 = get_level_pop(Z, z1+1, Te, dens)*ion_frac[z1]*abund

    # populations of this ion due to adjacent ions
    levpop_m1_casc = pyatomdb.apec.calc_ioniz_popn(levpop_m1, Z, z1, z1-1, Te, dens, datacache=d, do_xi=True)
    levpop_p1_casc = pyatomdb.apec.calc_recomb_popn(levpop_p1, Z, z1, z1+1, Te, dens, 0.0,0.0, datacache=d)

    # so... here we go.

    # within the same ion:
    rates_adj = rates * levpop[init]
    i = numpy.where(rates > 0)[0] # filter out the negatives, which are the diagonal terms
    rates_adj_filt=rates_adj[i]
    init_filt = init[i]
    final_filt = final[i]

    # total "out" rate
    rate_out_total = sum(rates_adj_filt[init_filt==up-1])
    # total "in" rate
    rate_in_total = sum(rates_adj_filt[final_filt==up-1])

    print("Excitation driven: total rate in %e; total rate out %e"%(rate_in_total, rate_out_total))

    # now look for major contributors
    init_la, final_la, rate_la=pyatomdb.apec.gather_rates(Z, z1, Te, dens, datacache=d,\
                               do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                               do_ir=False)
    i_check_la_out = numpy.where( (init_la==up-1) & (final_la==up-1))[0]
    init_la = init_la[rate_la > 0]
    final_la = final_la[rate_la > 0]
    rate_la = rate_la[rate_la > 0]

    init_ec, final_ec, rate_ec=pyatomdb.apec.gather_rates(Z, z1, Te, dens, datacache=d,\
                               do_la=False, do_ai=False, do_ec=True, do_pc=False,\
                               do_ir=False)

    # again, filter out the -ve values
    init_ec = init_ec[rate_ec > 0]
    final_ec = final_ec[rate_ec > 0]
    rate_ec = rate_ec[rate_ec > 0]

    print("population of level of interest:")
    print("  from excitation   : %e"%(levpop[up-1]))
    print("  from ionization   : %e"%(levpop_m1_casc[up-1]))
    print("  from recombination: %e"%(levpop_p1_casc[up-1]))
    print("  total             : %e"%(levpop_p1_casc[up-1]+levpop_m1_casc[up-1]+levpop[up-1]))

    # the total rate out for the recombination and ionization are only from la rates, and rate out == rate in, so...

    total_la_out = sum(rate_la[init_la==up-1])  # the total rate out of the level, unadjusted for population

    rate_ion = levpop_m1_casc[up-1]*total_la_out
    rate_recomb=levpop_p1_casc[up-1]*total_la_out
    print("Population rate from ionization %e"%(rate_ion))
    print("Population rate from recombination %e"%(rate_recomb))

    # for excitation this was pre-calculated
    print("Population rate from excitation %e"%(rate_in_total)) # this already has the population in it
    total_pop_rate = rate_ion+ rate_recomb +rate_in_total

    print("Total population rate for level %i: %e"%(up, total_pop_rate))

    # now find the rates we care about
    rate_la *= levpop[init_la]
    rate_ec *= levpop[init_ec]

    i_in_ec = numpy.where(final_ec==up-1)[0]
    i_in_la = numpy.where(final_la==up-1)[0]

    sum_in_ec = sum(rate_ec[i_in_ec])
    sum_in_la = sum(rate_la[i_in_la])

    i_out_ec = numpy.where(init_ec==up-1)[0]
    i_out_la = numpy.where(init_la==up-1)[0]

    sum_out_ec = sum(rate_ec[i_out_ec])
    sum_out_la = sum(rate_la[i_out_la])

    k = numpy.where(rate_la[i_in_la] > 0.0001*(total_pop_rate))[0]
    cascade_dict, cascade_frac = {}, []
    all_cascades = []
    for kk in k:
        all_cascades.append((init_la[i_in_la[kk]], final_la[i_in_la[kk]]))

        if rate_la[i_in_la[kk]] / total_pop_rate > 0.01:
            cascade_dict.update(
                {str(init_la[i_in_la[kk]]) + '$\\rightarrow$' + str(final_la[i_in_la[kk]]): str(
                    round((rate_la[i_in_la[kk]] / total_pop_rate) * 100,1)) + '\\%'})
            cascade_frac.append(rate_la[i_in_la[kk]] / total_pop_rate)

    print("rates into level by process:")
    print(" Excitation : %e  %f%%"%( sum_in_ec, 100*sum_in_ec/total_pop_rate))
    print(" Cascade : %e  %f%%"%( sum_in_la, 100*sum_in_la/total_pop_rate))
    print(" Recombination : %e  %f%%"%( rate_recomb, 100*rate_recomb/total_pop_rate))
    print(" Ionization : %e  %f%%"%( rate_ion, 100*rate_ion/total_pop_rate))

    #return contribution fractions to transition rate
    exc = sum_in_ec/total_pop_rate
    cascade = sum_in_la/total_pop_rate
    recombination = rate_recomb/total_pop_rate
    ionization = rate_ion/total_pop_rate
    other_cascades = cascade - sum(cascade_frac)
    print("Other cascades = cascade - sum(cascade_frac):")
    print(other_cascades, cascade, sum(cascade_frac))
    return exc, cascade, cascade_dict, other_cascades, recombination, ionization, all_cascades

def mc_popn(Z, max_ionize, max_recomb, Te, runs=100, distribution='f'):
    """
    Z : int
        element
    max_ionize : float or dict
        max ionization errors
        float if same error for all ions
        dict of {'z1': delta_r}
    max_recomb : float or dict
        max recombination errors
        float if same error for all ions
        dict of {'z1': delta_r}
    distribution : str
        'f' for flat distribution
        'g' for gaussian distribution

    Returns the minimum and max ion population
    for element Z at specified Te (in K) after
    Monte Carlo calculations to vary ionization
    and recombination rates."""

    mc_popn = numpy.zeros([runs, Z + 1])
    random_ion, random_rec = numpy.zeros([Z]), numpy.zeros([Z])

    if isinstance(max_ionize, dict):
        print("Filling in any empty z1 errors as avg of adjacent ions with supplied ionization errors")
        ionize_errors = get_errors(Z, max_ionize)
    elif isinstance(max_ionize, float):
        print("Using fractional error =", max_ionize, "for all ionization rates.")
        ionize_errors = [max_ionize] * (Z)
    if isinstance(max_recomb, dict):
        print("Filling in any empty z1 errors as avg of adjacent ions with supplied recombination errors")
        recomb_errors = get_errors(Z, max_recomb)
    elif isinstance(max_recomb, float):
        print("Using fractional error =", max_recomb, "for all recombination rates.")
        recomb_errors = [max_recomb] * (Z)

    ionlist, reclist = numpy.zeros([1, Z]), numpy.zeros([1, Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Te, False, Z=Z, z1=z1, extrap=True, datacache=d)
        ionlist[0, z1 - 1] = iontmp
        reclist[0, z1 - 1] = rectmp

    for run in range(runs):
        # get random ionization errors for each ion
        for i, max_error in enumerate(ionize_errors):
            if distribution.lower() == 'g':
                lower, upper = -2 * max_error, 2 * max_error
                mu, sigma = 0, max_error
                random_ion[i] = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=1)
            elif distribution.lower() == 'f':
                random_ion[i] = numpy.random.uniform(low=-max_error, high=max_error, size=1)
        # get random recombination errors
        for i, max_error in enumerate(recomb_errors):
            # use average if error not supplied for that ion
            if distribution.lower() == 'g':
                lower, upper = -2 * max_error, 2 * max_error
                mu, sigma = 0, max_error
                random_rec[i] = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=1)
            elif distribution.lower() == 'f':
                random_rec[i] = numpy.random.uniform(low=-max_error, high=max_error, size=1)

        # vary rates by a different random number for each ion
        iontmp = ionlist * 1.0
        rectmp = reclist * 1.0
        for z1, rand_ion, rand_rec in zip(range(1, Z + 1), random_ion, random_rec):
            rectmp[0, z1 - 1] *= (1 + random_rec[z1 - 1])
            iontmp[0, z1 - 1] *= (1 + random_ion[z1 - 1])
        newpopn = solve_ionrec([Te], iontmp, rectmp, Z)

        for z1 in range(1, Z + 2):
            mc_popn[run, z1 - 1] = newpopn[:, z1 - 1]

    # find mean and standard dev from all runs for each ion at each temp
    min, max = numpy.zeros([Z + 1]), numpy.zeros([Z + 1])
    for z1 in range(1, Z + 2):
        pop = mc_popn[:, z1 - 1]
        pop_list = numpy.sort(pop)
        min_index = int(0.16 * runs - 1)
        max_index = int(0.84 * runs - 1)
        min[z1 - 1] = pop_list[min_index]
        max[z1 - 1] = pop_list[max_index]

    return min, max

def calc_recomb_popn(which_transition, x, old_rate, levpop, Z, z1, z1_drv,T, dens, drlevrates, rrlevrates,\
                     newsettings2={}, settings=False, datacache=False, dronly=False,\
                     rronly=False):
  """
  Calculate the level population of a recombined ion

  Parameters
  ----------
  levpop: array(float)
    Level populations, already taking into account elemental abundance
    and ion fraction of z1_drv
  Z: int
  z1: int
  z1_drv: int
  T: electron temperature (K)
  dens: electron density (cm^-3)
  drlevrates: array(float)
    Rates into each level from DR calculations
  rrlevrates: array(float)
    Rates into each level from RR calculations

  Returns
  -------
  array(float)
    Level population
  """
  import scipy.sparse as sparse
  from scipy.sparse.linalg import spsolve

  # levpop at this point should alread have the corrected abundance in
  # there

  lvdat = pyatomdb.atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)
  if not lvdat:
    nlev = 1
    levpop = numpy.zeros(1, dtype=float)
    return levpop
  nlev = len(lvdat[1].data)
  Tarr, dummy = pyatomdb.util.make_vec(T)

  if nlev > pyatomdb.const.NLEV_NOSPARSE:

  # sort the levels
    aidat = pyatomdb.atomdb.get_data(Z,z1,'AI', settings=settings, datacache=datacache)
    if aidat:
      ailev = numpy.array(pyatomdb.util.unique(aidat[1].data['level_init']))-1
      nailev = len(ailev)
      isbound = numpy.ones(nlev, dtype=bool)
      isbound[ailev]=False
    else:
      nailev = 0
      isbound = numpy.ones(nlev, dtype=bool)

    recombrate = numpy.zeros(nlev, dtype=float)

  # get the recomb data

    irdat = pyatomdb.atomdb.get_data(Z, z1, 'IR',  settings=settings, datacache=datacache)

    for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
      if ir['TR_TYPE'] in ['RR','DR','XR']:
        recrate = pyatomdb.atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)*dens
        if not (numpy.isfinite(recrate)):
          pass
        else:
          recombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]

    maxlev = numpy.where(recombrate > 0)[0]
    if len(maxlev) == 0:  # so recombrate has all the influxes
      return numpy.zeros(nlev)
    maxlev=maxlev[-1]
    matrixB = recombrate

    matrixA = {}
    matrixA['init'], matrixA['final'], matrixA['rate']=\
    pyatomdb.apec.gather_rates(Z, z1, T, dens, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

    matrixB *= -1
    nlev = len(matrixB)

    # remove all matrix A from and to level 0
    i = (matrixA['final']>0) & (matrixA['init']>0)
    matrixA['final'] = matrixA['final'][i]
    matrixA['init'] = matrixA['init'][i]
    matrixA['rate'] = matrixA['rate'][i]

    # remove all matrix A from and to high levels
    i = (matrixA['final']<maxlev+1) & (matrixA['init']<maxlev+1)
    matrixA['final'] = matrixA['final'][i]
    matrixA['init'] = matrixA['init'][i]
    matrixA['rate'] = matrixA['rate'][i]

    # subtract 1 from the levels
    matrixA['init']-= 1
    matrixA['final']-= 1

    A  = numpy.zeros([maxlev,maxlev])
    for i in range(len(matrixA['final'])):
      A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]

    #update matrix A with new rate
    initial_lev, final_lev = which_transition[0], which_transition[1]
    A[final_lev - 1, initial_lev - 1] += (x - old_a)  # off diagonal term
    A[initial_lev - 1, initial_lev - 1] -= (x - old_a)  # diagonal term

    # if doing line ratio, get newsettings2
    if newsettings2 != {}:
        which_transition, x, old_rate = [newsettings2.get(k) for k in newsettings2]
        # update matrixA with new rate
        initial_lev, final_lev = which_transition[0], which_transition[1]
        A[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        A[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    levpop_this = numpy.zeros(len(matrixB))

    if sum(matrixB[1:] < 0):
      levpop_this[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])

  else:

    rrrecombrate = numpy.zeros(nlev, dtype=float)
    drrecombrate = numpy.zeros(nlev, dtype=float)
    irdat = pyatomdb.atomdb.get_data(Z, z1, 'IR', settings=settings, datacache=datacache)

    havedrrate=False
    haverrrate=False
    for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
      if ir['TR_TYPE'] in ['RR','XR']:

        recrate = pyatomdb.atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)

        rrrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*dens

        if ((ir['TR_TYPE'] in ['RR','XR']) & (ir['level_final']>1)):
          haverrrate=True
      if ir['TR_TYPE'] in ['DR','XD']:
        recrate = pyatomdb.atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)
        drrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*dens
        if ((ir['TR_TYPE'] in ['DR','XD']) & (ir['level_final']>1)):
          havedrrate=True

    if havedrrate:
      tmpdrlevrates=0.0
    else:
      tmpdrlevrates=drlevrates
    if haverrrate:
      tmprrlevrates=0.0
    else:
      tmprrlevrates=rrlevrates
    if isinstance(rrlevrates, numpy.ndarray):
      sumrrlevrates = sum(rrlevrates)
    else:
      sumrrlevrates = 0.0
    if isinstance(drlevrates, numpy.ndarray):
      sumdrlevrates = sum(drlevrates)
    else:
      sumdrlevrates = 0.0

    matrixB = rrrecombrate+drrecombrate+tmpdrlevrates+tmprrlevrates
    if dronly:
      matrixB = drrecombrate+tmpdrlevrates
    if rronly:
      matrixB = rrrecombrate+tmprrlevrates


    matrixA = numpy.zeros([nlev,nlev],dtype=float)


    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', settings=settings, datacache=datacache)

    matrixA_in = {}
    matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
    pyatomdb.apec.gather_rates(Z, z1, T, dens, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

    datacache={}
    for i in range(len(matrixA_in['init'])):
      matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]

    # update matrix A with new rate
    initial_lev, final_lev = which_transition[0], which_transition[1]
    matrixA[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
    matrixA[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    # if doing line ratio, get newsettings2
    if newsettings2 != {}:
        which_transition, x, old_rate = [newsettings2.get(k) for k in newsettings2]
        # update matrixA with new rate
        initial_lev, final_lev = which_transition[0], which_transition[1]
        matrixA[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        matrixA[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    # solve unless matrixB ==0
    if sum(matrixB[1:])>0:
      matrixB = -1*matrixB
      levpop_this = pyatomdb.apec.calc_cascade_population(matrixA, matrixB)
    else:
      levpop_this = numpy.zeros(nlev)

  return levpop_this

def calc_ioniz_popn(which_transition, x, old_rate, levpop, Z, z1, z1_drv,T, Ne, newsettings2={}, settings=False, \
                    datacache=False, do_xi=False):
  """
  Calculate the level population due to ionization into the ion

  Parameters
  ----------
  levpop: array(float)
    The level population of the parent ion. Should already have abundance
    and ion fraction built in.
  Z: int
  z1: int
  z1_drv: int
  T: float
  Ne: float
  settings: dict
  datacache: dict
  do_xi: bool
    Include collisional ionization

  Returns
  -------
  levpop_out: array(float)
    The level populations of the Z,z1 ion
  """

  # levpop at this point should alread have the corrected abundance in
  # there

  # This calculates the population of levels of z1, given a previous
  # ion's (z1-1) population of levpop.
  import scipy.sparse as sparse
  from scipy.sparse.linalg import spsolve

  #print("Starting calc_ioniz_popn at %s"%(time.asctime()))
  lvdat = pyatomdb.atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)

  # if we have no lv data, ignore.
  if not pyatomdb.util.keyword_check(lvdat):
    nlev = 1
    return numpy.array([0.0])
  nlev = len(lvdat[1].data)

  # get populating rate from previous ion
  aidat = pyatomdb.atomdb.get_data(Z, z1-1, 'AI', settings=settings, datacache=datacache)
  ionizrateai=numpy.zeros(nlev, dtype=float)
  ionizrateir=numpy.zeros(nlev, dtype=float)

  if aidat:
    tmp_pop = levpop[aidat[1].data['level_init']-1]
    for iai in range(len(aidat[1].data)):
      ionizrateai[aidat[1].data['level_final'][iai]-1] += \
               tmp_pop[iai]*aidat[1].data['auto_rate'][iai]

  if do_xi:

    irdat = pyatomdb.atomdb.get_data(Z, z1-1, 'IR', settings=settings, datacache=datacache)
    ionpot = float(irdat[1].header['ionpot'])
    if z1 >1:
      lvdatm1 = pyatomdb.atomdb.get_data(Z, z1-1, 'LV', settings=settings, datacache=datacache)
  # go through each excitation, have fun

    for iir, ir in enumerate(irdat[1].data):
      if ir['TR_TYPE'] in ['XI']:
        Te =  numpy.array([T])
        ionrate=pyatomdb.atomdb.get_maxwell_rate(Te, irdat, iir, lvdatm1, \
                                     lvdatap1=lvdat, ionpot=ionpot)
        ionizrateir[ir['level_final']-1] += levpop[ir['level_init']-1]*\
                                       ionrate


  ionizrate=ionizrateir+ionizrateai
  matrixB = ionizrate

  # save some time if there is nothing to ionize.

  if sum(matrixB[1:]) ==0:
    levpop_this = numpy.zeros(len(matrixB))
    return levpop_this

  maxlev = numpy.where(matrixB > 1e-40)[0]
  if len(maxlev)==0:
    print("No significant ionization found")
    popn = numpy.zeros(len(matrixB))
    return popn
  maxlev=maxlev[-1]

  matrixA_in={}
  matrixA_in['init'], matrixA_in['final'], matrixA_in['rate'] = \
   pyatomdb.apec.gather_rates(Z, z1, T, Ne, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=True, do_ec=False, do_pc=False,\
                 do_ir=False)

  i = (matrixA_in['init']<=maxlev) & (matrixA_in['final']<=maxlev)
  matrixA_in['init']=matrixA_in['init'][i]
  matrixA_in['final']=matrixA_in['final'][i]
  matrixA_in['rate']=matrixA_in['rate'][i]

  # fix the rates
  for i in range(len(matrixA_in['init'])):
    if matrixA_in['init'][i]==matrixA_in['final'][i]:
      if matrixA_in['rate'][i] >=0.0:
        matrixA_in['rate'][i] -=1e10

  if (maxlev <= pyatomdb.const.NLEV_NOSPARSE):
    matrixA = numpy.zeros([maxlev+1,maxlev+1], dtype=float)

    for i in range(len(matrixA_in['init'])):
      matrixA[matrixA_in['final'][i], matrixA_in['init'][i]] += matrixA_in['rate'][i]

    #update matrixA with new rate
    initial_lev, final_lev = which_transition[0], which_transition[1]
    matrixA[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
    matrixA[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    # if doing line ratio, get newsettings2
    if newsettings2 != {}:
        which_transition, x, old_rate = [newsettings2.get(k) for k in newsettings2]
        # update matrixA with new rate
        initial_lev, final_lev = which_transition[0], which_transition[1]
        A[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        A[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    # bug-u-fix
    for i in range(1, maxlev):
      if matrixA[i,i] >= 0:
        matrixA[i,i]=-1e10

    popn = numpy.zeros(nlev)

    matrixB*=-1

    try:
      popn[1:maxlev] = numpy.linalg.solve(matrixA[1:maxlev,1:maxlev], matrixB[1:maxlev])
    except numpy.linalg.linalg.LinAlgError:
      "EEK ERROR!"
      raise

  else:
    matrixA={}
    matrixB *= -1
    nlev = len(matrixB)

    if sum(matrixB)>=0:
      return numpy.zeros(len(matrixB))

    # remove ground level
    i = (matrixA_in['init']>0) & (matrixA_in['final']>0)

    matrixA['init'] = matrixA_in['init'][i]
    matrixA['final'] = matrixA_in['final'][i]
    matrixA['rate'] = matrixA_in['rate'][i]

    i = (matrixA['init']<=maxlev+1) & (matrixA['final']<=maxlev+1)

    matrixA['init'] = matrixA['init'][i]
    matrixA['final'] = matrixA['final'][i]
    matrixA['rate'] = matrixA['rate'][i]


    # subtract 1 from the levels
    matrixA['init']-=1
    matrixA['final']-=1

    A = numpy.zeros([maxlev, maxlev])
    for i in range(len(matrixA['final'])):
      A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]

    # update matrixA with new rate
    initial_lev, final_lev = which_transition[0], which_transition[1]
    A[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
    A[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    #if doing line ratio, get newsettings2
    if newsettings2 != {}:
        which_transition, x, old_rate = [newsettings2.get(k) for k in newsettings2]
        # update matrixA with new rate
        initial_lev, final_lev = which_transition[0], which_transition[1]
        A[final_lev - 1, initial_lev - 1] += (x - old_rate)  # off diagonal term
        A[initial_lev - 1, initial_lev - 1] -= (x - old_rate)  # diagonal term

    popn = numpy.zeros(len(matrixB))
    popn[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])

    popn_bak = popn*1.0

    popn[popn<0] = 0.0

  return popn

def line_sensitivity_csd(Z, z1, up, lo, vary, errors, trans_list, max_ionize, max_recomb, temps={}, trans_labels={}, show=True, show_legend=True, use_max_csd=True):
    """ Plots sensitivity of a line to perturbation
    in another transition and CSD errors.

    Z : int
        element
    z1 : int
        ion charge +1
    up : int
        upper level of sensitive transition
    lo : int
        lower level of sensitive transition
    vary : str
        'exc' or 'A'
    errors : list or tuple
        list of fractional errors to change rate by
    trans_list : list or tuple
        list of transitions (up, lo) to vary
        (can be just one)
    max_ionize : float or dict
        either a single fractional error or dict
        of {z1: error} pairs to change ionization rates by
    max_recomb : float or dict
        either a single fractional error or dict
        of {z1: error} pairs to change recomb rates by
    temps : list or tuple
        list of temperatures to calculate emissivity change at
        default is "5 temperatures of interest"
    show : bool
        True to show to screen, fig is saved regardless
    show_legend : bool
        True to label fig with legend (default)
    use_max_csd : bool
        use max values of perturbed CSD (default)
        set to False to use min values for lower limit
    """
    if temps == {}: temps = temps_of_interest(Z, z1)
    cmap = get_cmap(len(temps)+1, name=plot_settings['cmap'])
    fig, ax2 = plt.subplots()
    ax2.set_xlabel('Fractional error')
    ax2.set_ylabel('Fractional change in emissivity')
    temp_str, legend_labels = ['%.1E' % Decimal(x) for x in temps], []
    markers = ['o', 'v', 's', 'P', '^', '2']  # marker is for each transition

    # vary each transition in trans_list separately
    for a, transition in enumerate(trans_list):
        up, lo = transition[0], transition[1]

        # add legend label
        if trans_labels == {}:
            legend_labels.append(Line2D([0], [0], marker=markers[a], label=str(transition[0]) + '->' + str(transition[1])))
        else:
            legend_labels.append(Line2D([0], [0], marker=markers[a], label=trans_labels[a]))

        for b, Te in enumerate(temps):
            legend_labels.append(Line2D([0], [0], color=cmap(b), label=temp_str[b]))
            min, max = mc_popn(Z, max_ionize, max_recomb,Te, distribution='f')
            if use_max_csd == True: #use max possible csd values
                matrix2 = line_sensitivity(Z,z1,up,lo,vary,errors, trans_list, temps=Te, pop_fraction=max, show=False, savefig=False)
            else: # use min csd values
                matrix2 = line_sensitivity(Z, z1, up, lo, vary, errors, trans_list, temps=Te, pop_fraction=min, show=False, savefig=False)
            matrix3 = line_sensitivity(Z,z1,up,lo,vary,errors, trans_list, temps=Te, show=False, savefig=False)
            new_e, none = matrix2[0,0,:], matrix3[0,0,:]
            ax2.plot(errors, new_e, color=cmap(b))
            ax2.plot(errors, none, '--', color=cmap(b), marker='o')

    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    element, ion = pyatomdb.atomic.Ztoelsymb(Z), pyatomdb.atomic.int_to_roman(z1)
    fname = element + ' ' + ion + ' ' + str(up) + '-' + str(lo) + ' sensitivity with csd'
    num = get_file_num(fname, '.pdf')
    fig.suptitle("%d->%d sensitivity" % (up, lo))
    if show_legend == True:
        ax2.legend(fontsize=plot_settings['legendsize'], loc=plot_settings['legendloc'], handles=legend_labels)
    fig.savefig(fname+str(num)+'.pdf')
    if show == True:
        plt.show()

