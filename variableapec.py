"""
This module contains methods for looking at emission and line ratios
as a function of varying atomic data from the AtomDB files. Requires
PyAtomdB and Python 3."""

# Keri Heuer
# Version 2.0, March 25, 2020

import matplotlib.pyplot as plt, matplotlib.ticker as mtick, scip.stats as stats, \
pyatomdb, numpy, pickle, pathlib, csv, os, errno, haslib, requests, urllib.request, \
urllib.parse, urllib.error, time, subprocess, shutil, wget, glob, datetime, ftplib, \
pathlib, collections, operator, requests, matplotlib.pylab as pylab
from io import StringIO
from astropy.io import fits
from astropy.table import Table, Column
from matplotlib.offsetbox import AnchoredText
from decimal import Decimal
from matplotlib import cm
try:
  import astropy.io.fits as pyfits
except:
  import pyfits

def ionize(Z, z1, Te, dens, la_matrix, in_range, pop_fraction):
    
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1-1, Te, dens, do_la= True, \
                            do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)
    
    lvdat = pyatomdb.atomdb.get_data(Z, z1-1, 'LV', datacache=d)
    lvdat = lvdat[1].data
    nlev = len(lvdat)
    drv_matrix = numpy.zeros((nlev,nlev))
    drv_B = numpy.zeros(nlev)
    
    #populate full CR matrix by summing rates for all processes
    for i in range(len(init)):
        x = final[i]
        y = init[i]
        drv_matrix[x][y] += rates[i]
        
    #set up and solve CR matrix for level populations
    drv_matrix[0][:] = 1.0
    drv_B[0]=1.0
    drv_lev_pop = numpy.linalg.solve(drv_matrix,drv_B)
    
    z1_drv = z1-1
    #print(len(drv_lev_pop))
    ion_B = pyatomdb.apec.calc_ioniz_popn(drv_lev_pop, Z, z1, z1_drv, Te, dens, datacache=d)
    neg_ion_B = [-x for x in ion_B]
    neg_ion_B[0]=1.0
    ion_levpop = numpy.linalg.solve(la_matrix,neg_ion_B)

    ion_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    
    for i in range(len(in_range)):
        ion_linelist['lambda'][i] = in_range['WAVELEN'][i]
        ion_pop_level = in_range['UPPER_LEV'][i]
        ion_linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*ion_levpop[ion_pop_level-1]
        
    ion_linelist['epsilon']=[e*pop_fraction[z1_drv-1] for e in ion_linelist['epsilon']]
    
    return ion_linelist['epsilon']
            
def recombine(Z, z1, Te, dens, la_matrix, in_range, pop_fraction):
    # So, this is all completely correct...except:
    #  if z1+1 = Z+1, then you have no electrons, so the level population
    # for the recombining ion is meaningless. In this case, just set it
    # to be all in the ground state (whatever that would mean) and move on

    print(" I am now calling gather_rates for Z=%i, z1=%i" % (Z, z1 + 1))

    if z1 < Z:
        init, final, rates = pyatomdb.apec.gather_rates(Z, z1 + 1, Te, dens, do_la=True, \
                                                        do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)

        lvdat = pyatomdb.atomdb.get_data(Z, z1 + 1, 'LV', datacache=d)
        lvdat = lvdat[1].data
        nlev = len(lvdat)
        drv_matrix = numpy.zeros((nlev, nlev))
        drv_B = numpy.zeros(nlev)

        # populate full CR matrix by summing rates for all processes
        for i in range(len(init)):
            x = final[i]
            y = init[i]
            drv_matrix[x][y] += rates[i]

        # set up and solve CR matrix for level populations
        drv_matrix[0][:] = 1.0
        drv_B[0] = 1.0
        drv_lev_pop = numpy.linalg.solve(drv_matrix, drv_B)

    else:
        print("in here")
        # declare 1 level, fully populated
        nlev = 1
        drv_lev_pop = numpy.ones(1, dtype=float)

    z1_drv = z1 + 1
    recomb_B = pyatomdb.apec.calc_recomb_popn(drv_lev_pop, Z, z1, z1_drv, Te, dens, drlevrates=0, rrlevrates=0,
                                              datacache=d)
    neg_recomb_B = [-x for x in recomb_B]
    neg_recomb_B[0] = 1.0
    recomb_levpop = numpy.linalg.solve(la_matrix, neg_recomb_B)

    recomb_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))

    # calculate emissivity from recombination
    for i in range(len(in_range)):
        recomb_linelist['lambda'][i] = in_range['WAVELEN'][i]
        recomb_pop_level = in_range['UPPER_LEV'][i]
        recomb_linelist['epsilon'][i] = in_range['EINSTEIN_A'][i] * recomb_levpop[recomb_pop_level - 1]

    # multiply by fraction of plasma from z1_drv
    recomb_linelist['epsilon'] = [e * pop_fraction[z1_drv - 1] for e in recomb_linelist['epsilon']]

    return recomb_linelist['epsilon']

#def set_up(Z, z1, Te, dens, process, delta_r, transition, transition_2, \
           #npnts, wavelen, Te_range, dens_range, corrthresh, e_signif):
def set_up(Z, z1, Te, dens, extras={}):
    
    """ Uses inputs from check_sensitivity(), see routine for variable definitions as they are same throughout.
    
    Set_up() gets original rates, radiative collision matrix, level populations, line lists, 
    and intensities and returns these in a dictionary called values. Creates sensitivity tables.
    
    Returns inputs in dictionary called inputs, returns the variable transition."""

    # create data cache
    global d
    d = {}

    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= True, \
                            do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)
        
    lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV')
    lvdat = lvdat[1].data
    nlev = len(lvdat)
                   
    matrix = numpy.zeros((nlev,nlev))
    B = numpy.zeros(nlev)
    #populate full CR matrix by summing rates for all processes
    for i in range(len(init)):
        x = final[i]
        y = init[i]
        matrix[x][y] += rates[i]
        
    #set up and solve CR matrix for level populations
    matrix[0][:] = 1.0
    B[0]=1.0
    lev_pop = numpy.linalg.solve(matrix,B)

    #convert level populations into line lists & intensities for excitation only
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
    in_range = ladat[1].data
        
    linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype')) 
    
    upper_level= numpy.zeros(len(in_range))
    lower_level = numpy.zeros(len(in_range))
    
    for i in range(len(in_range)):
        linelist['lambda'][i] = in_range['WAVELEN'][i]
        pop_level = in_range['UPPER_LEV'][i]
        linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*lev_pop[pop_level-1]
        upper_level[i] = in_range['UPPER_LEV'][i]
        lower_level[i] = in_range['LOWER_LEV'][i]
        
    #resolve CR matrix "A" using LA only (for cascading emission from ionization/recombination)
    la_init, la_final, la_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= True, \
                            do_ec=False, do_ir=False, do_pc=False, do_ai=False, datacache=d)
    la_matrix = numpy.zeros((nlev,nlev))
    for i in range(len(la_init)):
        x = la_final[i]
        y = la_init[i]
        la_matrix[x][y] += la_rates[i]
    la_matrix[0][:] = 1.0
    
    #find fraction of each ion in plasma
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
    
    #set up complete line list (emiss only due to excitation at this point)
    full_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    full_linelist['lambda'] = linelist['lambda']
    full_linelist['epsilon']= linelist['epsilon']
    
    #now add emissivity from ionization and recombination to excitation linelist (depending on z1)
    if z1 == 1: #skip ionization
        recomb_emiss = recombine(Z, z1, Te, dens, la_matrix, in_range, pop_fraction)
        full_linelist['epsilon'] += recomb_emiss
    elif z1 == Z+1: #skip recombination
        ion_emiss = ionize(Z, z1, Te, dens, la_matrix, in_range, pop_fraction)
        full_linelist['epsilon'] += ion_emiss
    else: #do both
        recomb_emiss = recombine(Z, z1, Te, dens, la_matrix, in_range, pop_fraction)
        ion_emiss = ionize(Z, z1, Te, dens, la_matrix, in_range, pop_fraction)
        full_linelist['epsilon'] += recomb_emiss
        full_linelist['epsilon'] += ion_emiss
    
    #set up sensitivity & partial derivatives tables
    table = Table([full_linelist['lambda'], upper_level, lower_level, full_linelist['epsilon']], \
        names=('Lambda', 'Upper', 'Lower', 'Epsilon_orig'))
    new_table = Table([full_linelist['lambda'], upper_level, lower_level], names=('Lambda', 'Upper', 'Lower'))

    #save variables 
    inputs = {}
    inputs.update({'Z':Z, 'z1':z1, 'Te':Te, 'dens':dens})
    if extras != {}:
        process, delta_r, transition, transition_2, npnts, wavelen, \
        Te_range, dens_range, corrthresh, e_signif = [extras.get(k) for k in extras]
        inputs.update({'process': process, 'delta_r': delta_r,
        'transition': transition, 'transition_2': transition_2,
        'npnts': npnts, 'wavelen': wavelen, 'Te_range': Te_range, 'dens_range': dens_range,
         'corrthresh': corrthresh, 'e_signif': e_signif})
    values = {}
    values.update( {'matrix': matrix, 'B': B, 'in_range': in_range, 'linelist': linelist, 'table': table, 'new_table': new_table} )

    if extras != {}:
        return inputs, values, transition
    else: return inputs, values

def vary_a(inputs, values, which_transition):
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the A value for by the specified delta_r
    in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("****************************************************")
    print("in vary_a for " + str(which_transition) + ", calculating rate for Te=%e and dens=%e" %(Te, dens))
    print("****************************************************")
    
    initial_lev = which_transition[0]
    final_lev = which_transition[1]
    
    old_a = in_range['EINSTEIN_A'][(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)][0]
    a_index = numpy.where([(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)])[1][0]

    table['Epsilon_orig'].unit='A'

    if old_a == 0:
        old_a = 1e-40
    
    #vary A values
    min_a = 1-delta_r
    max_a = 1+delta_r
    new_a = [min_a*old_a, max_a*old_a]
    print("old A value =", old_a)
    print("new A values =", new_a)
    q_max = new_a[-1]
    q_min = new_a[0]
        
    index=1
    for x in new_a:
        #update LA data for specified transition
        in_range['EINSTEIN_A'][a_index] = x
        
        #get new CR matrix and resolve level pops with new A
        frac = str(round(x/old_a,2))
        new_matrix = matrix.copy()
        new_linelist=linelist.copy()
        new_matrix[final_lev-1, initial_lev-1] += (x-old_a)   #off diagonal term
        new_matrix[initial_lev-1, initial_lev-1] -= (x-old_a)   #diagonal term
        new_matrix[0][:] = 1.0
        
        #find new level populations and get new epsilon values from excitation
        new_lev_pop = numpy.linalg.solve(new_matrix,B)

        for i in range(len(in_range)):
            new_linelist['lambda'][i] = in_range['WAVELEN'][i] #wavelength will be same
            pop_level = in_range['UPPER_LEV'][i]
            new_linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*new_lev_pop[pop_level-1]

        # resolve CR matrix "A" using LA only (for cascading emission from ionization/recombination)
        la_init, la_final, la_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=True, \
                                        do_ec=False, do_ir=False, do_pc=False, do_ai=False, datacache=d)
        lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV', datacache=d)
        lvdat = lvdat[1].data
        nlev = len(lvdat)
        la_matrix = numpy.zeros((nlev, nlev))
        for i in range(len(la_init)):
            x = la_final[i]
            y = la_init[i]
            la_matrix[x][y] += la_rates[i]

        #update la_matrix with new A value
        la_matrix[final_lev-1, initial_lev-1] += (x-old_a)
        la_matrix[initial_lev-1, initial_lev-1] -= (x-old_a)

        la_matrix[0][:] = 1.0

        # find fraction of each ion in plasma
        pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

        # now add emissivity from ionization and recombination to full linelist (depending on z1)
        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, la_matrix,
                                     in_range, pop_fraction)
            new_linelist['epsilon'] += recomb_emiss
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, la_matrix,
                               in_range, pop_fraction)
            new_linelist['epsilon'] += ion_emiss
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, la_matrix,
                                     in_range, pop_fraction)
            ion_emiss = ionize(Z, z1, Te, dens, la_matrix,
                               in_range, pop_fraction)
            new_linelist['epsilon'] += recomb_emiss
            new_linelist['epsilon'] += ion_emiss
        
        #update sensitivity table 
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' A')
        table.add_columns([new_col])
            
        index+=1
        
    values = {}
    values.update( {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, \
                    'q_max': q_max, 'q_min': q_min} )
    
    return inputs, values
    
def vary_exc(inputs, values, which_transition):
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the excitation rate for by the specified
    delta_r in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("***********************************************************************************")
    print("in vary_exc for " + str(which_transition) + ", calculating rate for Te=%e and dens=%e" %(Te, dens))
    print("***********************************************************************************")
    
    exc_init, exc_final, exc_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= False, \
                        do_ec=True, do_ir=False, do_pc=False, do_ai=False, datacache=d)
    
    initial_lev=which_transition[0]
    final_lev=which_transition[1]

    for a,b,c in zip(exc_init, exc_final, exc_rates):
        if (a,b) == (which_transition[0]-1, which_transition[1]-1):
            old_rate = c

    try:
        if old_rate == 0:
            old_rate = 1e-40
    except UnboundLocalError:
        print("Could not find transition", which_transition, " - please check input transition levels")

    table['Epsilon_orig'].unit='orig rate'
        
    #vary rate
    min_rate = 1-delta_r
    max_rate = 1+delta_r
    new_rate = [min_rate*old_rate, max_rate*old_rate]
    print("old rate =", old_rate)
    print("new rates =", new_rate)
    q_max = new_rate[-1]
    q_min = new_rate[0]

    index=1
    for x in new_rate:
        #loop through varied rates, get new matrix and resolve level pops
        frac = str(round(x/old_rate,2))
        new_matrix = matrix.copy()
        new_linelist=linelist.copy()
        new_matrix[final_lev-1, initial_lev-1] += (x-old_rate)   #off diagonal term
        new_matrix[initial_lev-1, initial_lev-1] -= (x-old_rate)   #diagonal term
        new_matrix[0][:] = 1.0
        
        #find new level populations and get new epsilon values
        new_lev_pop = numpy.linalg.solve(new_matrix,B)
        
        for i in range(len(in_range)):
            new_linelist['lambda'][i] = in_range['WAVELEN'][i] #wavelength will be same
            pop_level = in_range['UPPER_LEV'][i]
            new_linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*new_lev_pop[pop_level-1]

        # resolve CR matrix "A" using LA only (for cascading emission from ionization/recombination)
        la_init, la_final, la_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=True, \
                            do_ec=False, do_ir=False, do_pc=False, do_ai=False, datacache=d)
        lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV', datacache=d)
        lvdat = lvdat[1].data
        nlev = len(lvdat)
        la_matrix = numpy.zeros((nlev, nlev))
        for i in range(len(la_init)):
            x = la_final[i]
            y = la_init[i]
            la_matrix[x][y] += la_rates[i]
        la_matrix[0][:] = 1.0

        # find fraction of each ion in plasma
        pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

        # now add emissivity from ionization and recombination to full linelist (depending on z1)
        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, la_matrix,
                                     in_range, pop_fraction)
            new_linelist['epsilon'] += recomb_emiss
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, la_matrix,
                               in_range, pop_fraction)
            new_linelist['epsilon'] += ion_emiss
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, la_matrix,
                                     in_range, pop_fraction)
            ion_emiss = ionize(Z, z1, Te, dens, la_matrix,
                               in_range, pop_fraction)
            new_linelist['epsilon'] += recomb_emiss
            new_linelist['epsilon'] += ion_emiss
        
        #update sensitivity table 
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' rate')
        table.add_columns([new_col])

        index +=1
    
    values = {}
    values.update( {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, \
                    'q_max': q_max, 'q_min': q_min} )
    
    return inputs, values

def get_tables(inputs, values):
    
    """ Inputs and values are dictionaries outputted by either vary_exc() or vary_a().
    
    Calculates dE/dR as the difference in min and max emissivities divided by the difference in changed rates,
    dE/dE_orig as the difference in min and max emissivities divided by the original emissivity, and sorts
    the sensitivity table by dE/dE_orig. 
    
    If the variables corrthresh and e_signif were specified by the user, it will filter the
    sensitivity table for values with dE/dE_orig greater than corrthresh and/or for lines
    with a dE/dR greater than e_signif.
    
    Prints the sensitivity table.
    
    Returns table and new_table and the dictionaries inputs and results."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    table, new_table, new_linelist, q_max, q_min = [values.get(k) for k in values]

    epsilon_max = table['Epsilon_'+str(npnts)]
    
    ratio_col = Column(name = 'Epsilon_max/Epsilon_min', data = (epsilon_max/table['Epsilon_1']))
    table.add_columns([ratio_col])
    table['Epsilon_max/Epsilon_min'].unit = None
        
    #add partial derivative, dE/dR
    delta_epsilon = Column(name = 'dE/dR', data = (epsilon_max-table['Epsilon_1'])/(q_max-q_min))
    new_table.add_columns([delta_epsilon])
    new_table['dE/dR'].unit = None
    
    orig_eps = Column(name = 'Epsilon_orig', data = table['Epsilon_orig'])
    new_table.add_columns([orig_eps])
    
    #add "correlation factor" dE/dE_orig
    epsilon_corr = Column(name = 'dE/dE_orig', data = (epsilon_max-table['Epsilon_1'])/table['Epsilon_orig'])
    abs_epsilon_corr = Column(name = '|dE/dE_orig|', data = [abs(val) for val in epsilon_corr])
    new_table.add_columns([abs_epsilon_corr])
    try:
      new_table.sort('|dE/dE_orig|', reverse=True)
    except TypeError:
      new_table.sort('|dE/dE_orig|')
      new_table = new_table[::-1]
    
    #apply filters
    if corrthresh != 0.0:     #only show lines whose "epsilon correlation" >= than specified value
        new_table = new_table[new_table['|dE/dE_orig|'] >= corrthresh]
    elif e_signif != 0.0:    #only show lines with partial epsilon/partial rate derivative is >= specified value
        new_table = new_table[new_table['dE/dR'] >= eps_der]
    
    #update output dictionary
    results={}
    results.update({'inputs': inputs,
                    'wavelength': numpy.array(new_table['Lambda']),\
                    'upper': numpy.array(new_table['Upper']), \
                    'lower': numpy.array(new_table['Lower']), \
                    'dE/dR': numpy.array(new_table['dE/dR']), \
                    'epsilon_orig': numpy.array(new_table['Epsilon_orig']),\
                    '|dE/dE_orig|': numpy.array(new_table['|dE/dE_orig|']), \
                    'min_eps': numpy.array(table['Epsilon_1']), \
                    'max_eps': epsilon_max})

    return table, new_table, inputs, results

def plot_nei(inputs):
    
    """ Gets the non-equilibrium line emissivities for each driving ion and plots. Inputs
    is a dictionary outputted by set_up()."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    #get nei line emissivities and plot
    if process == 'A':
        emiss = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition[0], transition[1])
    elif process == 'exc':
        emiss = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition[1], transition[0])
        
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(transition[0])+'-'+str(transition[1])
    name=element+'_'+ion+'_'+process+'_'+line+'_'

    axs[0,0].set_xlabel('Temperature in keV', fontsize=12)
    axs[0,0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    axs[0,0].tick_params(labelsize=12)
    
    for key in emiss.keys():
        if isinstance(key, numpy.int32) == True:
            if (int(key) < z1):
                number=str(int(key)-1)
                label='Ionization of '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label)
            if (int(key) == z1):
                number=str(int(key)-1)
                label='Excitation from '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label, linewidth=3)
            if (isinstance(key, numpy.int32) == True) & (int(key) > z1):
                number=str(int(key)-1)
                label='Recombination of '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label)

    axs[0,0].legend(fontsize = 'x-small')
    
def wrapper_plot_nei(inputs):   #get nei line emissivities and plot
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    if process == 'A':
        emiss = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition[0], transition[1])
    elif process == 'exc':
        emiss = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition[1], transition[0])
    orig_text = '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1])
        
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(transition[0])+'-'+str(transition[1])
    name=element+'_'+ion+'_'+process+'_'+line+'_'
    
    axs[0,0].set_xlabel('Temperature in keV', fontsize=12)
    axs[0,0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    axs[0,0].tick_params(labelsize=12)
    anchored_text = AnchoredText(orig_text, loc='upper right', frameon=False)
    axs[0,0].add_artist(anchored_text)
    
    for key in emiss.keys():
        if isinstance(key, numpy.int32) == True:
            if (int(key) < z1):
                number=str(int(key)-1)
                label='Ionization of '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label)
            if (int(key) == z1):
                number=str(int(key)-1)
                label='Excitation from '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label, linewidth=3)
            if (isinstance(key, numpy.int32) == True) & (int(key) > z1):
                number=str(int(key)-1)
                label='Recombination of '+element+' '+number+'+'
                axs[0,0].loglog(emiss['Te'], emiss[key], label=label)
            
    axs[0,0].legend(fontsize = 'x-small')
    
    #now repeat for transition_2
    if process == 'A':
        emiss_2 = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition_2[0], transition_2[1])
    elif process == 'exc':
        emiss_2 = pyatomdb.spectrum.get_nei_line_emissivity(Z, z1, transition_2[1], transition_2[0])

    text = '{0:g}'.format(transition_2[0]) + '->' + '{0:g}'.format(transition_2[1])
    
    axs[0,1].set_xlabel('Temperature in keV', fontsize=12)
    axs[0,1].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    axs[0,1].tick_params(labelsize=12)
    anchored_text = AnchoredText(text, loc='upper right', frameon=False)
    axs[0,1].add_artist(anchored_text)
    
    for key in emiss_2.keys():
        if isinstance(key, numpy.int32) == True:
            if (int(key) < z1):
                number=str(int(key)-1)
                label='Ionization of '+element+' '+number+'+'
                axs[0,1].loglog(emiss_2['Te'], emiss_2[key], label=label)
            if (int(key) == z1):
                number=str(int(key)-1)
                label='Excitation from '+element+' '+number+'+'
                axs[0,1].loglog(emiss_2['Te'], emiss_2[key], label=label, linewidth=3)
            if (isinstance(key, numpy.int32) == True) & (int(key) > z1):
                number=str(int(key)-1)
                label='Recombination of '+element+' '+number+'+'
                axs[0,1].loglog(emiss_2['Te'], emiss_2[key], label=label)
            
    axs[0,1].legend(fontsize = 'x-small')
    
def plot_sensitivity(inputs, new_table):
    
    """ Inputs is a dictionary outputted by set_up() and new_table is the sensitivity
    table outputted by get_tables().
    
    Plots the lines affected by the parameter(s) changed, including the original transition.
    Currently plotting for wavelengths between 10-20 Angstroms (as set by mask1 and mask2 below).""" 
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(transition[0])+'-'+str(transition[1])
    name=element+'_'+ion+'_'+process+'_'+line+'_'

    #set up plot of "sensitive epsilons"
    axs[0,1].set_xlabel('Wavelength ($\AA$)', fontsize=12)
    axs[0,1].set_ylabel('% Emissivity Change', fontsize=12)
    axs[0,1].tick_params(labelsize=12)
    
    #filter data for significance
    cutoff_data = new_table[new_table['|dE/dE_orig|'] >= corrthresh]
    if wavelen != {}:
        cutoff_data = cutoff_data[(wavelen[0] < cutoff_data['Lambda']) & (cutoff_data['Lambda'] < wavelen[1])]
    
    #plot wavelength vs. % emissivity change
    axs[0,1].semilogy(cutoff_data['Lambda'], cutoff_data['|dE/dE_orig|'], linestyle='', marker='o', label=transition)
    
    #label each point w/ transition
    transition_labels=[]
    for x in cutoff_data:
        if process == 'A':
            transition_name = '{0:g}'.format(x['Upper'])+'->'+'{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name = '{0:g}'.format(x['Lower'])+'->'+'{0:g}'.format(x['Upper'])
        transition_labels.append(transition_name)
    for (x,y,label) in zip(numpy.array(cutoff_data['Lambda']),numpy.array(cutoff_data['|dE/dE_orig|']),transition_labels):
        axs[0,1].annotate(label, xy=(x,y))

    axs[0,1].legend(fontsize='x-small')
    
def wrapper_plot_sensitivity(inputs, new_table1, new_table2):
    
    """ Results is a dictionary.
    
    Plots lines affected by the changed parameter(s) including the original transition,
    for the multiple transitions specified, in a different color for each transition.
    
    By default, will plot for wavelengths between 10-20 Angstroms, unless a wavelength
    range is specified by wavelen=()."""

    #plot sensitive epsilons for first transition
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]

    axs[1, 0].set_xlabel('Wavelength ($\AA$)', fontsize=12)
    axs[1, 0].set_ylabel('% Emissivity Change', fontsize=12)
    axs[1, 0].tick_params(labelsize=12)

    # filter data for significance
    cutoff_data = new_table1[new_table1['|dE/dE_orig|'] >= corrthresh]
    cutoff_data = cutoff_data[(wavelen[0] < cutoff_data['Lambda'])]
    cutoff_data = cutoff_data[cutoff_data['Lambda'] < wavelen[1]]

    # plot wavelength vs. % emissivity change
    axs[1, 0].semilogy(cutoff_data['Lambda'], cutoff_data['|dE/dE_orig|'], linestyle='', marker='o', label=transition)

    # label each point w/ transition
    transition_labels = []
    for x in cutoff_data:
        if process == 'A':
            transition_name = '{0:g}'.format(x['Upper']) + '->' + '{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name = '{0:g}'.format(x['Lower']) + '->' + '{0:g}'.format(x['Upper'])
        transition_labels.append(transition_name)
    for (x, y, label) in zip(numpy.array(cutoff_data['Lambda']), numpy.array(cutoff_data['|dE/dE_orig|']),
                             transition_labels):
        axs[1, 0].annotate(label, xy=(x, y))

    axs[1, 0].legend(fontsize='x-small')

    #now plot sensitive epsilons for second transition
    # Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    # npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs_2.get(n) for n in inputs_2]

    axs[1, 1].set_xlabel('Wavelength ($\AA$)', fontsize=12)
    axs[1, 1].set_ylabel('% Emissivity Change', fontsize=12)
    axs[1, 1].tick_params(labelsize=12)

    # filter data for significance
    cutoff_data2 = new_table2[new_table2['|dE/dE_orig|'] >= corrthresh]
    cutoff_data2 = cutoff_data2[(wavelen[0] < cutoff_data2['Lambda'])]
    cutoff_data2 = cutoff_data2[cutoff_data2['Lambda'] < wavelen[1]]

    # plot wavelength vs. % emissivity change
    axs[1, 1].semilogy(cutoff_data2['Lambda'], cutoff_data2['|dE/dE_orig|'], linestyle='', marker='o', label=transition_2)

    # label each point w/ transition
    transition_labels = []
    for x in cutoff_data2:
        if process == 'A':
            transition_name = '{0:g}'.format(x['Upper']) + '->' + '{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name = '{0:g}'.format(x['Lower']) + '->' + '{0:g}'.format(x['Upper'])
        transition_labels.append(transition_name)
    for (x, y, label) in zip(numpy.array(cutoff_data2['Lambda']), numpy.array(cutoff_data2['|dE/dE_orig|']),
                             transition_labels):
        axs[1, 1].annotate(label, xy=(x, y))

    axs[1, 1].legend(fontsize='x-small')

        
def run_line_diagnostics(table, inputs, values, which_transition, type, num={}):
    
    """ Table is table of transitions and emissivity values. Inputs and values are dictionaries and
    which_transition is a tuple specifying which transition to run line diagnostics for.
    
    Varies temperature and density separately and recalculates emissivitiy values. 
    
    Returns the dictionary line_diagnostics of arrays of the temperature and density bins used,
    the original, min, and max emissivities, the name of element and ion, label of what process varied
    and by how much, and which_transition."""

    #if number of points not specified by user, set defaults
    if num == {}:
        temp_num = 20
        dens_num = 8
    else:
        temp_num = num
        dens_num = num

    #read in inputs and values from set_up()
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]

    extras={}
    extras.update({'process': process, 'delta_r': delta_r,
                   'transition': transition, 'transition_2': transition_2,
                   'npnts': npnts, 'wavelen': wavelen, 'Te_range': Te_range, 'dens_range': dens_range,
                   'corrthresh': corrthresh, 'e_signif': e_signif})

    print("////////////////////////////////////////////////////////////////")
    print("Now running line diagnostics for transition:", which_transition)
    print("////////////////////////////////////////////////////////////////")
    
    table = values.get('table')
    new_table = values.get('new_table')
    full_linelist = values.get('full_linelist')
    q_max = values.get('q_max')
    q_min = values.get('q_min')
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(which_transition[0])+'-'+str(which_transition[1])
    name=element+'_'+str(z1)+'_'
    percentage=str(delta_r*100)+'%'
    
    line_diagnostics = {}

    #check type of diagnostics to run
    if type == 'both':
        temperature = True
        density = True
    if type == 'temp':
        temperature = True
        density = False
    if type == 'dens':
        temperature = False
        density = True

    if temperature == True:
        # vary temperature and recalculate emissivities
        print("\nChanging temperature now.\n")
        temp_bins = list(numpy.geomspace(Te_range[0], Te_range[1], num=temp_num))  # 20
        Te_eps_orig, Te_eps_min, Te_eps_max =[], [], []
        counter=0
        for temp_Te in temp_bins:
            print("Current temperature is:", temp_Te)
            Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, extras=extras)

            if process == 'A':
                Te_new_inputs, Te_new_values = vary_a(Te_inputs, Te_values, which_transition)
                Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
                for x in Te_table:
                    if (x['Upper'], x['Lower']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        Te_eps_orig.append(x['Epsilon_orig'])
                        Te_eps_min.append(x['Epsilon_1'])
                        Te_eps_max.append(x['Epsilon_'+str(npnts)])

            elif process == 'exc':
                Te_new_inputs, Te_new_values = vary_exc(Te_inputs, Te_values, which_transition)
                Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)

                for x in Te_table:
                    if (x['Upper'], x['Lower']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        Te_eps_orig.append(x['Epsilon_orig'])
                        Te_eps_min.append(x['Epsilon_1'])
                        Te_eps_max.append(x['Epsilon_'+str(npnts)])

            counter += 1
            print(str(temp_num-counter), 'temperatures left for', element, ion)
            #print('\nCurrent Te arrays:\n', Te_eps_orig, '\n', Te_eps_min, '\n', Te_eps_max, '\n')


    if density == True:
        #vary density and recalculate emissivities
        print("\nChanging density now.\n")

        dens_bins = list(numpy.geomspace(dens_range[0], dens_range[1], num=dens_num)) #8

        dens_eps_orig, dens_eps_min, dens_eps_max =[],[],[]
        counter=0
        for temp_dens in dens_bins:
            dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, extras=extras)
            if process == 'A':
                dens_new_inputs, dens_new_values = vary_a(dens_inputs, dens_values, which_transition)
                dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)

                for x in dens_table:
                    if (x['Lower'], x['Upper']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        dens_eps_orig.append(x['Epsilon_orig'])
                        dens_eps_min.append(x['Epsilon_1'])
                        dens_eps_max.append(x['Epsilon_'+str(npnts)])

            elif process == 'exc':
                dens_new_inputs, dens_new_values = vary_exc(dens_inputs, dens_values, which_transition)
                dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)

                for x in dens_table:
                    if (x['Upper'], x['Lower']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        dens_eps_orig.append(x['Epsilon_orig']/ temp_dens)
                        dens_eps_min.append(x['Epsilon_1']/temp_dens)
                        dens_eps_max.append(x['Epsilon_'+str(npnts)]/temp_dens)

            counter += 1
            print(str(dens_num-counter), 'densities left for', element, ion)
    
    if process == 'A':
        label = 'Range due to $\Delta $\pm$' + percentage + ' in A value'
    elif process == 'exc':
        label = 'Range due to $\Delta $\pm$' + percentage + ' in direct excitation rate'
        
    if type == 'temp':
        print("Ran temp diagnostics.")
        line_diagnostics.update({'temp_bins': numpy.asarray(temp_bins),
        'Te_eps_orig': numpy.asarray(Te_eps_orig),
        'Te_eps_min': numpy.asarray(Te_eps_min),
        'Te_eps_max': numpy.asarray(Te_eps_max),
        'name': name, 'label': label, 'transition': which_transition})
    elif type == 'dens':
        print("Ran dens diagnostics.")
        line_diagnostics.update({'dens_bins': numpy.asarray(dens_bins),
        'dens_eps_orig': numpy.asarray(dens_eps_orig),
        'dens_eps_min': numpy.asarray(dens_eps_min),
        'dens_eps_max': numpy.asarray(dens_eps_max),
        'name': name, 'label': label, 'transition': which_transition})
    elif type == 'both':
        print("Ran temp and dens diagnostics.")
        line_diagnostics.update({'temp_bins': numpy.asarray(temp_bins), \
        'dens_bins': numpy.asarray(dens_bins), \
        'Te_eps_orig': numpy.asarray(Te_eps_orig), \
        'Te_eps_min': numpy.asarray(Te_eps_min), \
        'Te_eps_max': numpy.asarray(Te_eps_max), \
        'dens_eps_orig': numpy.asarray(dens_eps_orig), \
        'dens_eps_min': numpy.asarray(dens_eps_min), \
        'dens_eps_max': numpy.asarray(dens_eps_max), \
        'name': name, 'label': label, \
        'transition': which_transition})

    return line_diagnostics

def wrapper_plot_line_diagnostics(inputs, line_diagnostics, line_diagnostics_2):
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    temp_bins, dens_bins, Te_eps_orig, Te_eps_min, Te_eps_max, dens_eps_orig, \
            dens_eps_min, dens_eps_max, name, label, transition = [line_diagnostics.get(k) for k in line_diagnostics]

    orig_text='{0:g}'.format(transition[0])+'->'+'{0:g}'.format(transition[1])

    #plot emissivity versus temperature
    axs2[0,0].tick_params(labelsize=12)
    axs2[0,0].set_xlabel('Temperature in K', fontsize=12)
    axs2[0,0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    anchored_text = AnchoredText(orig_text, loc='lower right', prop=dict(size=8), frameon=False)
    axs2[0,0].add_artist(anchored_text)
    axs2[0,0].loglog(temp_bins, Te_eps_orig, label='Original')
    axs2[0,0].legend(fontsize='x-small')
    axs2[0,0].fill_between(temp_bins, Te_eps_min, Te_eps_max, alpha=0.5, color='g', \
                     label="Range")
    axs2[0,0].legend(fontsize='x-small', loc='upper left')
    
    #plot emissivity versus density
    axs2[0,1].tick_params(labelsize=12)
    axs2[0,1].set_xlabel('Density in cm$^{-3}$', fontsize=12)
    axs2[0,1].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    anchored_text = AnchoredText(orig_text, loc='lower right', prop=dict(size=8), frameon=False)
    axs2[0,1].add_artist(anchored_text)
    axs2[0,1].loglog(dens_bins, dens_eps_orig, label='Original')
    axs2[0,1].legend(fontsize='x-small')
    axs2[0,1].fill_between(dens_bins, dens_eps_min, dens_eps_max, alpha=0.5, color='g', \
                        label='Range')
    axs2[0,1].legend(fontsize='x-small', loc='upper left')
    
    #now redefine variables from line_diagnostics for transition_2
    #and plot on separate subplots

    temp_bins_2, dens_bins_2, Te_eps_orig_2, Te_eps_min_2, Te_eps_max_2, dens_eps_orig_2, \
            dens_eps_min_2, dens_eps_max_2, name_2, label_2, transition2 = [line_diagnostics_2.get(n) for n in line_diagnostics_2]
    orig_text_2 = '{0:g}'.format(transition_2[0]) + '->' + '{0:g}'.format(transition_2[1])

    #plot emissivity versus temperature
    axs2[1,0].tick_params(labelsize=12)
    axs2[1,0].set_xlabel('Temperature in K', fontsize=12)
    axs2[1,0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    anchored_text = AnchoredText(orig_text_2, loc='lower right', prop=dict(size=8), frameon=False)
    axs2[1,0].add_artist(anchored_text)
    axs2[1,0].loglog(temp_bins_2, Te_eps_orig_2, label='Original')
    axs2[1,0].legend(fontsize='x-small')
    axs2[1,0].fill_between(temp_bins_2, Te_eps_min_2, Te_eps_max_2, alpha=0.5, color='g', \
                     label="Range")
    axs2[1,0].legend(fontsize='x-small', loc='upper left')
    
    #plot emissivity versus density
    #text = orig_text + ', ' + temp
    axs2[1,1].tick_params(labelsize=12)
    axs2[1,1].set_xlabel('Density in cm$^{-3}$', fontsize=12)
    axs2[1,1].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=12)
    anchored_text = AnchoredText(orig_text_2, loc='lower right', prop=dict(size=8), frameon=False)
    axs2[1,1].add_artist(anchored_text)
    axs2[1,1].loglog(dens_bins_2, dens_eps_orig_2, label='Original')
    axs2[1,1].legend(fontsize='x-small')
    axs2[1,1].fill_between(dens_bins_2, dens_eps_min_2, dens_eps_max_2, alpha=0.5, color='g', \
                        label='Range')
    axs2[1,1].legend(fontsize='x-small', loc='upper left')
    
def plot_line_diagnostics(inputs, line_diagnostics):
    
    """ Line_diagnostics is a dictionary outputted by run_line_diagnostics. It holds arrays of
    the temperature and density bins, the original, minimum, and maximum emissivities from varying
    temperature and then density, the name of the element and ion, label of which process varied
    and by how much, and the transition.
    
    Plots emissivity as a function of temperature and density for the specified single transition."""

    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]


    temp_bins, dens_bins, Te_eps_orig, Te_eps_min, Te_eps_max, dens_eps_orig, \
            dens_eps_min, dens_eps_max, name, label, transition = [line_diagnostics.get(k) for k in line_diagnostics]

    #plot emissivity versus temperature
    axs[1,0].tick_params(labelsize=12)
    axs[1,0].set_xlabel('Temperature in K', fontsize=14)
    axs[1,0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=14)
    axs[1,0].loglog(temp_bins, Te_eps_orig, label='Original')
    axs[1,0].legend(fontsize='x-small')
    axs[1,0].fill_between(temp_bins, Te_eps_min, Te_eps_max, alpha=0.5, color='g', \
                     label="Range")
    axs[1,0].legend(fontsize='x-small')
    
    #plot emissivity versus density
    axs[1,1].tick_params(labelsize=12)
    axs[1,1].set_xlabel('Density in cm$^{-3}$', fontsize=14)
    axs[1,1].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=14)
    axs[1,1].loglog(dens_bins, dens_eps_orig, label='Original')
    axs[1,1].legend(fontsize='x-small')
    axs[1,1].fill_between(dens_bins, dens_eps_min, dens_eps_max, alpha=0.5, color='g', \
                        label='Range')
    axs[1,1].legend(fontsize='x-small')
    
    
def run_line_ratio_diagnostics(line_diagnostics_1, line_diagnostics_2):
    
    """ Table1 and table2 are tables from individually run get_tables() on the two transitions
    specified by the user. Inputs1, values1, inputs2, values2, are  dictionaries holding the
    inputs and the sensitivity table/emissivity values for each transition respectively.
    
    Varies temperature and density separately for each transition and recalculates emissivity.
    
    Returns dictionary line_ratio_diagnostics containing arrays of the temperature and density bins
    for each transition, as well as the original, minimum, and maximum line ratios calculated
    from varying temperature and density independently."""

    line_ratio_diagnostics = {}
    
    temp_bins1, dens_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1, dens_eps_orig1, \
        dens_eps_min1, dens_eps_max1, name1, label1, transition1 = [line_diagnostics_1.get(k) for k in line_diagnostics_1]
    temp_bins2, dens_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2, dens_eps_orig2, \
        dens_eps_min2, dens_eps_max2, name2, label2, transition2 = [line_diagnostics_2.get(k) for k in line_diagnostics_2]

    line1=str(transition1[0])+'-'+str(transition1[1])
    line2=str(transition2[0])+'-'+str(transition2[1])
    line_ratio = line1+'_'+line2

    print("////////////////////////////////////////////////////////////////")
    print("Running diagnostics for line ratio", line_ratio)
    print("////////////////////////////////////////////////////////////////")

    Te_line_ratios = [x/y for x,y in zip(Te_eps_orig1, Te_eps_orig2)]
    Te_line_ratios_min = [x/y for x,y in zip(Te_eps_min1, Te_eps_max2)]
    Te_line_ratios_max = [x/y for x,y in zip(Te_eps_max1, Te_eps_min2)]
    dens_line_ratios = [x/y for x,y in zip(dens_eps_orig1, dens_eps_orig2)]
    dens_line_ratios_min = [x/y for x,y in zip(dens_eps_min1, dens_eps_max2)]
    dens_line_ratios_max = [x/y for x,y in zip(dens_eps_max1, dens_eps_min2)]
    
    line_ratio_diagnostics.update({'temp_bins1': numpy.asarray(temp_bins1), \
                                   'temp_bins2': numpy.asarray(temp_bins2), \
                                   'dens_bins1': numpy.asarray(dens_bins1), \
                                    'dens_bins2': numpy.asarray(dens_bins2), \
                                 'Te_line_ratios': numpy.asarray(Te_line_ratios),  \
                                 'Te_line_ratios_min': numpy.asarray(Te_line_ratios_min), \
                                 'Te_line_ratios_max': numpy.asarray(Te_line_ratios_max), \
                                 'dens_line_ratios': numpy.asarray(dens_line_ratios), \
                                 'dens_line_ratios_min': numpy.asarray(dens_line_ratios_min), \
                                 'dens_line_ratios_max': numpy.asarray(dens_line_ratios_max), \
                                   'ratio': line_ratio, 'name': name1, 'label': label1})

    t1 = Table([temp_bins1, Te_line_ratios_min, Te_line_ratios, Te_line_ratios_max],
               names=('temp', 'ratio min', 'ratio orig', 'ratio max'))
    t2 = Table([dens_bins1, dens_line_ratios_min, dens_line_ratios, dens_line_ratios_max],
               names=('dens', 'ratio min', 'ratio orig', 'ratio max'))
    for number in range(1, 20, 1):
        file = pathlib.Path(name1 + '_line_ratios_Te_' + str(number) + '.fits')
        if file.exists():
            continue
        else:
            t1.write(name1 + 'line_ratios_Te_' + str(number) + '.fits', format='fits')
            t2.write(name1 + 'line_ratios_dens_' + str(number) + '.fits', format='fits')
            break

    return line_ratio_diagnostics
    
def plot_line_ratio_diagnostics(inputs, line_ratio_diagnostics):
    
    """ Line_ratio_diagnostics is a dictionary from run_line_ratio_diagnostics() containing
    arrays of temperature and density bins, and line ratios (original, min, max) calculated
    from varying temperature and density.
    
    Plots the line ratios of emissivities for specified two transitions as a function of
    temperature and density."""

    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    #file = open('line_ratio_diagnostics', 'rb')
    #line_ratio_diagnostics = pickle.load(file)
    #file.close()
    
    temp_bins1, temp_bins2, dens_bins1, dens_bins2, Te_line_ratios, Te_line_ratios_min, \
               Te_line_ratios_max, dens_line_ratios, dens_line_ratios_min, \
               dens_line_ratios_max, ratio, name, label = \
               [line_ratio_diagnostics.get(k) for k in line_ratio_diagnostics]

    text = '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1]) + \
        ' / ' + '{0:g}'.format(transition_2[0]) + '->' + '{0:g}'.format(transition_2[1])

    #plot emissivity versus temperature
    axs2[2,0].tick_params(labelsize=12)
    axs2[2,0].set_xlabel('Temperature in K', fontsize=12)
    axs2[2,0].set_ylabel('Line Ratio\n'+text, fontsize=12)
    axs2[2,0].semilogx(temp_bins1, Te_line_ratios, label='Original')
    axs2[2,0].fill_between(temp_bins1, Te_line_ratios_min, Te_line_ratios_max, alpha=0.5, color='g', \
                     label="Range")
    axs2[2,0].legend(fontsize='x-small')
    
    #plot emissivity versus density
    axs2[2,1].tick_params(labelsize=12)
    axs2[2,1].set_xlabel('Density in cm$^{-3}$', fontsize=12)
    axs2[2,1].set_ylabel('Line Ratio\n'+text, fontsize=12)
    axs2[2,1].semilogx(dens_bins1, dens_line_ratios, label='Original')
    axs2[2,1].legend(fontsize='x-small')
    axs2[2,1].fill_between(dens_bins1, dens_line_ratios_min, dens_line_ratios_max, alpha=0.5, color='g', label='Range')
    axs2[2,1].legend(fontsize='x-small')

    plt.tight_layout()

def wrapper_check_sensitivity(transition_list):
    """
    Runs entire check_sensitivity() routine on a list [] of transitions

    :param transition_list:
    :return:
    """

    for transition in transition_list:
        print("Checking sensitivity for transition", transition)
        check_sensitivity(Z, z1, Te, dens, process, delta_r, \
            transition, transition_2, npnts, wavelen, corrthresh, e_signif)

def plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, type, A_lines={}, exc_lines={}, wavelen=(10,20), corrthresh=10e-5):

    """ Plots sensitive epsilons for multiple transitions all on one plot.
    A_lines and exc_ lines are dictionaries of ('name': (initial, final)) """

    if type == 'A':
        processes = ['A']
        set_lines = [A_lines]
    elif type == 'exc':
        processes = ['exc']
        set_lines = [exc_lines]
    elif type == 'both':
        processes = ['A', 'exc']
        set_lines = [A_lines, exc_lines]

    Te_range = (Te / 10, Te * 10)
    dens_range = (10e0, 10e16)



    for set, process in zip(set_lines, processes):
        plt.figure()
        plt.xlabel('Wavelength ($\AA$)', fontsize=14)
        plt.ylabel('% Emissivity Change', fontsize=14)
        plt.tick_params(labelsize=12)
        temp = '%.E' % Decimal(Te) + 'K'
        percentage = '{0:g}'.format(delta_r * 100) + '%'
        density = 'dens=' + str(dens)
        element = pyatomdb.atomic.Ztoelsymb(Z)
        ion = pyatomdb.atomic.int_to_roman(z1)
        if process == 'A':
            text = element + ' ' + ion + ' ' + ', ' + temp + ', ' + density + '\n' + \
                   ' A value $\Delta$ $\pm$' + percentage
        elif process == 'exc':
            text = element + ' ' + ion + ' ' + ', ' + temp + ', ' + density + '\n' + \
                   ' Direct excitation rate $\Delta$ $\pm$' + percentage

        #plt.title(text, fontsize=16)

        for k in set:
            transition = set.get(k)
            extras = {'process': process, 'delta_r': delta_r,
                      'transition': transition, 'transition_2': None, 'npnts': 2,
                      'wavelen': wavelen, 'Te_range': Te_range,
                      'dens_range': dens_range, 'corrthresh': corrthresh, 'e_signif': 0.0}
            inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
            if process == 'A':
                new_inputs, new_values = vary_a(inputs, values, transition)
            elif process == 'exc':
                new_inputs, new_values = vary_exc(inputs, values, transition)
            table, new_table, inputs, results = get_tables(new_inputs, new_values)

            # filter data for significance
            cutoff_data = new_table[new_table['|dE/dE_orig|'] >= corrthresh]
            cutoff_data = cutoff_data[(wavelen[0] < cutoff_data['Lambda'])]
            cutoff_data = cutoff_data[cutoff_data['Lambda'] < wavelen[1]]

            # plot wavelength vs. % emissivity change
            plt.semilogy(cutoff_data['Lambda'], cutoff_data['|dE/dE_orig|'], linestyle='', marker='o',
                         label=k)

            # label each point w/ transition
            transition_labels = []
            for x in cutoff_data:
                if process == 'A':
                    transition_name = '{0:g}'.format(x['Upper']) + '->' + '{0:g}'.format(x['Lower'])
                if process == 'exc':
                    transition_name = '{0:g}'.format(x['Lower']) + '->' + '{0:g}'.format(x['Upper'])
                transition_labels.append(transition_name)
            for (x, y, label) in zip(numpy.array(cutoff_data['Lambda']), numpy.array(cutoff_data['|dE/dE_orig|']),
                                     transition_labels):
                plt.annotate(label, xy=(x, y))

            plt.tight_layout()
            plt.legend(fontsize='x-small')

    file_name = 'sensitive lines '+process+' '+element + str(z1) + '_'
    for number in range(1, 20, 1):
        file = pathlib.Path(file_name + str(number) + '.pdf')
        if file.exists():
            continue
        else:
            plt.savefig(file_name+ str(number) + '.pdf')
            break

    plt.show()

def check_sensitivity(Z, z1, Te, dens, process, delta_r, transition, transition_2={}, \
            npnts=2, wavelen={}, Te_range={}, dens_range={}, corrthresh={}, e_signif={}):

    """
    Check emissivity sensitivity for specified element, ion, and transition.
    
    Parameters
    ----------
    Z: int
    The nuclear charge of the element
    
    z1 : int
    ion charge +1
    
    Te : float
    temperature (Kelvin)
    
    dens : float
    electron density (cm^-3)
    
    process : str
    specify rate to vary, i.e. 'A' or 'exc'
    
    delta_r : float
    delta of rate to vary, i.e. delta_r = 0.1 varies rate by +-0.1
    
    transition : tuple
    (upper, lower) transition to vary
    i.e. to vary 5->4 A value, transition = (5,4)
        
    transition_2 : tuple
    if a second transition is provided, will run line ratio diagnostics
    
    npnts : int
    number of points to calculate emissivity for. default is 2

    wavelen : tuple
    range of wavelengths to plot sensitive lines over in Angstroms
    default is (10-20)

    Te_range : tuple
    range of temperatures to run line diagnostics on
    default is (Te/10, Te*10)

    dens_range : tuple
    range of densities to run line diagnostics on
    default is (10e0, 10e16)
    
    corrthresh : float
    the minimum desired correlation threshold for epsilon, dE/dE_orig
    default is 10e-5
    
    e_signif : float
    the minimum value of the partial derivative of epsilon to rate
        i.e. significance of dE/dR
  
    Returns
    -------
    if one transition specified: inputs, new_table, line_diagnostics
    if line ratio specified: inputs, new_table1, new_table2, line_diagnostics_1, \
    line_diagnostics_2, line_ratio_diagnostics
    
    """
    #set defaults
    if npnts == {}: npnts = 2
    if transition_2 == {}: transition_2 = []
    if corrthresh == {}: corrthresh = 10e-5
    if e_signif == {}: e_signif = 0.0
    if wavelen == {}: wavelen = (10, 20)
    if Te_range == {}: Te_range = (Te/10, Te*10)
    if dens_range == {}: dens_range = (10e0, 10e16)

    extras = {'process':process, 'delta_r':delta_r,
             'transition':transition, 'transition_2': transition_2, 'npnts':2,
              'wavelen':wavelen, 'Te_range':Te_range,
                    'dens_range':dens_range,'corrthresh':corrthresh, 'e_signif':e_signif}

    print("Z=" + str(Z), "z1=" + str(z1), "Te=" + str(Te), "dens=" + str(dens), "process=" + str(process), \
          "delta_r=" + str(delta_r), "transition=" + str(transition), "transition2=" + str(transition_2), \
          "npnts=" + str(npnts), "wavelength range=" + str(wavelen), "temp range=" + str(Te_range), \
          "dens range=" + str(dens_range), "correlation threshold=" + str(corrthresh),
          "epsilon significance=" + str(e_signif))

    if transition_2 =={}:    #check sensitivity for a single transition
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        if process == 'A':
            new_inputs, new_values = vary_a(inputs, values, transition)
        elif process == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)
        
        line_diagnostics = run_line_diagnostics(table, inputs, values, transition, type='both', num={})

        #put all graphs on one plot
        global fig
        global axs
        fig, axs = plt.subplots(2,2,figsize=(10,6))
        temp = '%.E' % Decimal(Te) + 'K'
        percentage = '{0:g}'.format(delta_r * 100) + '%'
        density = 'dens=' + str(dens)
        element = pyatomdb.atomic.Ztoelsymb(Z)
        ion = pyatomdb.atomic.int_to_roman(z1)
        if process == 'A':
            text = element+' '+ion+' ' + ', '+temp + ', ' + density + '\n' + \
                   '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1]) + \
                        ' A value $\Delta$ $\pm$' + percentage
        elif process == 'exc':
            text = element+' '+ion+' ' + ', '+temp + ', ' + density + '\n' + \
                   '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1]) + \
                        ' Direct excitation rate $\Delta$ $\pm$' + percentage

        fig.suptitle(text, fontsize=16)
        plot_nei(inputs)
        plot_sensitivity(inputs, new_table)
        plot_line_diagnostics(inputs, line_diagnostics)
        fig.tight_layout()
        fig.subplots_adjust(top=0.86)
        
        #file = open('results_'+process, 'wb')
        #pickle.dump([results, line_diagnostics], file)
        #file.close()

        file_name = element + str(z1) + '_' + process + '_' + \
                    str(transition[0]) + '-' + str(transition[1])
        for number in range(1, 20, 1):
            file = pathlib.Path(file_name + '_data_' + str(number) + '.pkl')
            if file.exists():
                continue
            else:
                fig.savefig(file_name + '_diagnostics_' + str(number) + '.pdf')

                f = open(file_name + '_data_' + str(number) + '.pkl', 'wb')
                pickle.dump(all_data, f)
                f.close()
                break

        # print sensitivity table
        if process == 'exc':
            print("\nFor", element, ion + ", changed inputs:", "excitation rate from level", \
                  str(transition[0]) + "->" + str(transition[1]), percentage)
            print("Lines affected at Te=" + str(Te) + ", dens=" + str(dens))
            print(new_table)
        if process =='A':
            print("\nFor", element, ion + ", changed inputs:", "A value for level", \
                  str(transition[0]) + "->" + str(transition[1]), percentage)
            print("Lines affected at Te=" + str(Te) + ", dens=" + str(dens))
            print(new_table)

        plt.show()

        all_data = {'inputs':inputs, 'table':table, 'new_table':new_table, 'line_diagnostics':line_diagnostics}
        return all_data

    elif transition_2 != {}:  #calculate diagnostics for a line ratio
        #which_transition = transition
        #extras.update({'transition':which_transition})
        print("Doing calculations for", which_transition)
        inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
        if process == 'A':
            new_inputs, new_values = vary_a(inputs, values, transition)
        elif process == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, transition)

        which_transition = transition_2
        #transition_2 = None
        extras.update({'transition':transition_2, 'transition_2':None})
        print("Doing calculations for", which_transition, "and transition 2 is", transition_2)
        inputs_2, values_2, transition_2 = set_up(Z, z1, Te, dens, extras=extras)
        if process == 'A':
            new_inputs_2, new_values_2, = vary_a(inputs_2, values_2, transition_2)
        elif process == 'exc':
            new_inputs_2, new_values_2, = vary_exc(inputs_2, values_2, transition_2)

        table1, new_table1, inputs1, results1 = get_tables(new_inputs, new_values)
        table2, new_table2, inputs2, results2 = get_tables(new_inputs_2, new_values_2)

        #run line diagnostics for individual transitions
        line_diagnostics_1 = run_line_diagnostics(table1, inputs, values, transition, type='both', num={})
        line_diagnostics_2 = run_line_diagnostics(table2, inputs_2, values_2, transition_2, type='both', num={})
        
        #run line ratio diagnostics
        line_ratio_diagnostics = run_line_ratio_diagnostics(line_diagnostics_1, line_diagnostics_2)

        #set up first page of plots
        fig, axs = plt.subplots(2,2, figsize=(10,6))
        temp = '%.E' % Decimal(Te) + 'K'
        percentage = '{0:g}'.format(delta_r * 100) + '%'
        density = 'dens=' + str(dens)
        element = pyatomdb.atomic.Ztoelsymb(Z)
        ion = pyatomdb.atomic.int_to_roman(z1)
        if process == 'A':
            text = element + ' ' + ion + ' ' + ', ' + temp + ', ' + density + '\n' + \
                   ' A value $\Delta$ $\pm$' + percentage
        elif process == 'exc':
            text = element + ' ' + ion + ' ' + ', ' + temp + ', ' + density + '\n' + \
                   ' Direct excitation rate $\Delta$ $\pm$' + percentage


        #plot nei line emissitivities individually
        wrapper_plot_nei(inputs)
        #plot sensitivity for both transitions individually
        wrapper_plot_sensitivity(inputs, new_table1, new_table2)

        fig.suptitle(text, fontsize=16)
        fig.tight_layout()
        fig.subplots_adjust(top=0.86)

        #set up second page of plots
        global fig2
        global axs2
        fig2, axs2 = plt.subplots(3,2, figsize=(10,6))
        fig2.suptitle(text, fontsize=16)

        #plot line diagnostics for each transition individually
        wrapper_plot_line_diagnostics(inputs, line_diagnostics_1, line_diagnostics_2)
        #now plot line ratio diagnostics
        plot_line_ratio_diagnostics(inputs, line_ratio_diagnostics)

        fig2.tight_layout()
        fig2.subplots_adjust(top=0.84, wspace=0.54, hspace=0.59)

        #print sensitivity tables
        if process == 'exc':
            print("\nFor", element, ion + ", changed inputs:", "excitation rate from level", \
              str(transition[0]) + "->" + str(transition[1]), "+/-", percentage)
            print("Lines affected at Te=" + str(Te) + ", dens=" + str(dens))
            print(new_table1)
            print("\nFor", element, ion + ", changed inputs:", "excitation rate from level", \
                  str(transition_2[0]) + "->" + str(transition_2[1]), "+/-", percentage)
            print("Lines affected at Te=" + str(Te) + ", dens=" + str(dens))
            print(new_table2)
        if process == 'A':
            print("\nFor", element, ion + ", changed inputs:", "A value for level", \
                  str(transition[0]) + "->" + str(transition[1]), "+/-", percentage)
            print("Lines affected at Te=" + str(Te) + ", dens=" + str(dens))
            print(new_table1)
            print("\nFor", element, ion + ", changed inputs:", "A value for level", \
                  str(transition_2[0]) + "->" + str(transition_2[1]), "+/-", percentage)
            print("Lines affected at Te="+str(Te) + ", dens="+str(dens))
            print(new_table2)

        all_data = {'inputs':inputs, 'table1':table1, 'table2':table2, 'new_table1':new_table1, \
                    'new_table2':new_table2, 'line_diagnostics_1': line_diagnostics_1, \
                    'line_diagnostics_2':line_diagnostics_2, 'line_ratio_diagnostics': line_ratio_diagnostics}

        # save plots and all_data
        file_name = element + '_' + ion + '_' + process + '_' + \
                    str(transition[0]) + '-' + str(transition[1]) + '_' + \
                    str(transition_2[0]) + '-' + str(transition_2[1])
        for number in range(1, 20, 1):
            file = pathlib.Path('data_' + str(number) + '.pkl')
            if file.exists():
                continue
            else:
                fig.savefig(file_name + '_sensitivity_' + str(number) + '.pdf')
                fig2.savefig(file_name + '_diagnostics_' + str(number) + '.pdf')
                break

        # save plots and all_data
        file_name = element + '_' + ion + '_' + process + '_' + \
                    str(transition[0]) + '-' + str(transition[1]) + '_' + \
                    str(transition_2[0]) + '-' + str(transition_2[1])
        for number in range(1, 20, 1):
            file = pathlib.Path('data_' + str(number) + '.pkl')
            if file.exists():
                continue
            else:
                fig.savefig(file_name + '_sensitivity_' + str(number) + '.pdf')
                fig2.savefig(file_name + '_diagnostics_' + str(number) + '.pdf')

                f = open(file_name+'_data_'+str(number)+'.pkl', 'wb')
                pickle.dump(all_data, f)
                f.close()
                break

        plt.show()
        plt.close('all')

        return all_data

def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]

def blended_line_ratio(Z, z1, Te, dens, process, delta_r, transition_list, denom, type={}, num={}, Te_range={}, dens_range={}):
    #specify equation of blended line ratio
    #i.e. denom=1, then blended line ratio = [line 1 + line 2] / line 3
    #denm=2, then blended line ratio = line 1 / [line 2 + line 3]

    emiss1, emiss2, emiss3 = {}, {}, {}
    emiss_list = (1,2,3)
    transition_2, npnts, wavelen, corrthresh, e_signif = {}, 2, (10, 20), 10e-5, 0.0

    if type == 'temp':
        for transition, emiss in zip(transition_list, emiss_list):
            extras={'process':process, 'delta_r':delta_r,
            'transition':transition, 'transition_2':transition_2, 'npnts':npnts,
            'wavelen':wavelen, 'Te_range':Te_range, 'dens_range':dens_range,
            'corrthresh':corrthresh, 'e_signif':e_signif}
            print('calculating for', transition)

            inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
            if process == 'A':
                new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
            elif process == 'exc':
                new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)
            table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

            lines=[]
            for upper, lower, wavelen in zip(table['Upper'], table['Lower'], table['Lambda']):
                if (upper, lower) == transition:
                    lines.append(wavelen)

            diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type='temp', num=num)
            temp_bins, Te_eps_orig, Te_eps_min, Te_eps_max, name, label, transition = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}
            if emiss == 2:
                emiss2 = {'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}
            if emiss == 3:
                emiss3 = {'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}

        # add emissivities from line diagnostics
        # line ratio = [ line 1 + line 2 ] / line 3

        temp_bins1, Te_eps_min, Te_eps_max, Te_eps_orig = [emiss1.get(k) for k in emiss1]

        temp_bins2, Te_eps_min2, Te_eps_max2, Te_eps_orig2 = [emiss2.get(k) for k in emiss2]

        temp_bins3, Te_eps_min3, Te_eps_max3, Te_eps_orig3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            Te_total_min = Te_eps_min / (Te_eps_min2 + Te_eps_min3)
            Te_total_max = Te_eps_max / (Te_eps_max2 + Te_eps_max3)
            Te_total_orig = Te_eps_orig / (Te_eps_orig2 + Te_eps_orig3)
        elif denom==1:
            Te_total_min = (Te_eps_min + Te_eps_min2) / Te_eps_min3
            Te_total_max = (Te_eps_max + Te_eps_max2) / Te_eps_max3
            Te_total_orig = (Te_eps_orig + Te_eps_orig2) / Te_eps_orig3

        t1 = Table([temp_bins1, Te_total_min, Te_total_max, Te_total_orig], \
                   names=('temp', 'Te_min', 'Te_max', 'Te_orig'))
        for number in range(1, 20, 1):
            file = pathlib.Path(name + 'blended_ratio_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t1.write(name + 'blended_ratio_Te_' + str(number) + '.fits', format='fits')
                break

        # plot
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        ax2.set_xlabel('Log T(K)', fontsize=12)
        ax2.set_ylabel('New Ratio/Original', fontsize=12)
        blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
        ax1.ylabel(blended_ratio, fontsize=12)
        plt.semilogx(temp_bins, Te_total_orig, label='Original')
        plt.fill_between(temp_bins, Te_total_min, Te_total_max, label='Range', color='g', alpha=0.5)
        min = [x / y for x, y in zip(Te_total_min, Te_total_orig)]
        max = [x / y for x, y in zip(Te_total_max, Te_total_orig)]
        ax2.fill_between(temp_bins, min, max, color='g', alpha=0.5, label='Range')
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)
        ax2.axhline(y=1)
        ax1.legend()
        ax2.legend()

        for number in range(1, 10):
            outfile = pathlib.Path(element + ' ' + str(z1) + ' Te blended line ratio ' + str(number) + '.pdf')
            if outfile.exists():
                continue
            else:
                plt.savefig(element + ' ' + str(z1) + ' Te blended ratio' + str(number) + '.pdf')
                break

    elif type == 'dens':
        for transition, emiss in zip(transition_list, emiss_list):
            extras = {'process': process, 'delta_r': delta_r,
                      'transition': transition, 'transition_2': transition_2, 'npnts': npnts,
                      'wavelen': wavelen, 'Te_range': Te_range, 'dens_range': dens_range,
                      'corrthresh': corrthresh, 'e_signif': e_signif}
            print('calculating for', transition)

            inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)

            if process == 'A':
                new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
            elif process == 'exc':
                new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)

            table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

            lines = []
            for upper, lower, wavelen in zip(table['Upper'], table['Lower'], table['Lambda']):
                if process == 'exc':
                    if (lower, upper) == transition:
                        lines.append(wavelen)
                elif process == 'A':
                    if (upper, lower) == transition:
                        lines.append(wavelen)

            diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type='dens', num=num)

            dens_bins, dens_eps_orig, dens_eps_min, dens_eps_max, name, label, which_transition = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'dens_bins': dens_bins, 'dens_min': dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
            if emiss == 2:
                emiss2 = {'dens_bins': dens_bins, 'dens_min':dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
            if emiss == 3:
                emiss3 = {'dens_bins': dens_bins, 'dens_min':dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}

        # add emissivities from line diagnostics
        # line ratio = [ line 1 + line 2 ] / line 3

        dens_bins1, dens_eps_min, dens_eps_max, dens_eps_orig = [emiss1.get(k) for k in emiss1]
        dens_bins2, dens_eps_min2, dens_eps_max2, dens_eps_orig2 = [emiss2.get(k) for k in emiss2]
        dens_bins3, dens_eps_min3, dens_eps_max3, dens_eps_orig3 = [emiss3.get(k) for k in emiss3]

        dens_total_min = (dens_eps_min + dens_eps_min2) / dens_eps_min3
        dens_total_max = (dens_eps_max + dens_eps_max2) / dens_eps_max3
        dens_total_orig = (dens_eps_orig + dens_eps_orig2) / dens_eps_orig3

        t2 = Table([dens_bins1, dens_total_min, dens_total_max, dens_total_orig],
                   names=('dens', 'dens_min', 'dens_max', 'dens_orig'))
        t3 = Table([dens_bins1, dens_eps_min,dens_eps_max, dens_eps_orig, dens_bins2, dens_eps_min2, dens_eps_max2, dens_eps_orig2,
                    dens_bins3, dens_eps_min3, dens_eps_max3, dens_eps_orig3],
                   names=('dens bins1', 'dens min1', 'dens max1', 'dens orig1',
                          'dens bins2',  'dens min2', 'dens max2', 'dens orig2',
                          'dens bins3',  'dens min3', 'dens max3', 'dens orig3'))
        for number in range(1, 20, 1):
            file = pathlib.Path(name + 'blended_ratio_dens_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t2.write(name+ 'blended_ratio_dens_' + str(number) + '.fits', format='fits')
                t3.write(name +'blended_individual_emiss_dens_' + str(number) + '.fits', format='fits')
                break

        # plot
        blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'

        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        ax2.set_xlabel('Log Density (cm$^{-3}$)', fontsize=12)
        ax2.set_ylabel('New Ratio/Original', fontsize=12)
        ax1.ylabel(blended_ratio, fontsize=12)
        ax1.semilogx(dens_bins, dens_total_orig, label='Original')
        ax1.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='g', alpha=0.5)
        ax2.axhline(y=1)
        min = [x / y for x, y in zip(dens_total_min, dens_total_orig)]
        max = [x / y for x, y in zip(dens_total_max, dens_total_orig)]
        ax2.fill_between(dens_bins, min, max, color='g', alpha=0.5, label='Range')
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)
        ax1.legend()
        ax2.legend()

        for number in range(1, 10):
            outfile = pathlib.Path(element + ' ' + str(z1) + ' dens blended line ratio ' + str(number) + '.pdf')
            if outfile.exists():
                continue
            else:
                plt.savefig(element + ' ' + str(z1) + ' dens blended ratio' + str(number) + '.pdf')
                break

    elif type == 'both':
        temp_bins, dens_bins, Te_eps_orig, Te_eps_min, Te_eps_max, dens_eps_orig, \
        dens_eps_min, dens_eps_max, name, label, transition = [diagnostics.get(k) for k in diagnostics]

        if emiss == emiss1:
            emiss1 = {'temp_bins': temp_bins, 'dens_bins': dens_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max,
                      'Te_orig': Te_eps_orig, \
                      'dens_min': dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
        if emiss == emiss2:
            emiss2 = {'temp_bins': temp_bins, 'dens_bins': dens_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max,
                      'Te_orig': Te_eps_orig, \
                      'dens_min': dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
        if emiss == emiss3:
            emiss3 = {'temp_bins': temp_bins, 'dens_bins': dens_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max,
                      'Te_orig': Te_eps_orig, \
                      'dens_min': dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}

        # add emissivities from line diagnostics
        # line ratio = [ line 1 + line 2 ] / line 3

        temp_bins1, dens_bins1, Te_eps_min, Te_eps_max, Te_eps_orig, dens_eps_min, dens_eps_max, dens_eps_orig = \
            [emiss1.get(k) for k in emiss1]

        temp_bins2, dens_bins2, Te_eps_min2, Te_eps_max2, Te_eps_orig2, dens_eps_min2, dens_eps_max2, dens_eps_orig2 = \
            [emiss2.get(k) for k in emiss2]

        temp_bins3, dens_bins3, Te_eps_min3, Te_eps_max3, Te_eps_orig3, dens_eps_min3, dens_eps_max3, dens_eps_orig3 = \
            [emiss3.get(k) for k in emiss3]

        Te_total_min = (Te_eps_min + Te_eps_min2) / Te_eps_min3
        Te_total_max = (Te_eps_max + Te_eps_max2) / Te_eps_max3
        Te_total_orig = (Te_eps_orig + Te_eps_orig2) / Te_eps_orig3
        dens_total_min = (dens_eps_min + dens_eps_min2) / dens_eps_min3
        dens_total_max = (dens_eps_max + dens_eps_max2) / dens_eps_max3
        dens_total_orig = (dens_eps_orig + dens_eps_orig2) / dens_eps_orig3

        t1 = Table([temp_bins1, Te_total_min, Te_total_max, Te_total_orig], \
                   names=('temp', 'Te_min', 'Te_max', 'Te_orig'))
        t2 = Table([dens_bins1, dens_total_min, dens_total_max, dens_total_orig],
                   names=('dens', 'dens_min', 'dens_max', 'dens_orig'))
        for number in range(1, 20, 1):
            file = pathlib.Path(name + 'blended_ratio_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t1.write(name+ 'blended_ratio_Te_' + str(number) + '.fits', format='fits')
                t2.write(name+ 'blended_ratio_dens_' + str(number) + '.fits', format='fits')
                break

        # plot
        fig1 = plt.figure()
        fig1.xlabel('Log T(K)', fontsize=12)
        blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
        fig1.ylabel('Blended line ratio ' + blended_ratio, fontsize=12)
        fig1.semilogx(temp_bins, Te_total_orig, label='Original')
        fig1.fill_between(temp_bins, Te_total_min, Te_total_max, label='Range', color='g', alpha=0.5)

        fig2 = plt.figure()
        fig2.xlabel('Log Density (cm$^{-3}$)', fontsize=12)
        fig2.ylabel('Blended line ratio ' + blended_ratio, fontsize=12)
        fig2.semilogx(dens_bins, dens_total_orig, label='Original')
        fig2.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='g', alpha=0.5)

        for number in range(1,10):
            outfile = pathlib.Path(element+' '+str(z1)+' Te blended line ratio '+str(number)+'.pdf')
            if outfile.exists():
                continue
            else:
                fig1.savefig(element + ' ' + str(z1) + ' Te blended line ratio ' + str(number) + '.pdf')
                fig2.savefig(element + ' ' + str(z1) + ' dens blended line ratio ' + str(number) + '.pdf')
                break

    plt.show()
    plt.close('all')

def plot_blended_ratio(type, file):
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    ax1.set_ylabel('Blended line ratio', fontsize=12)
    ax2.set_ylabel('New Ratio/Original', fontsize=12)

    if type == 'dens':
        with fits.open(file) as hdul:
            data = hdul[1].data
            dens_bins, dens_total_min, dens_total_max, dens_total_orig = data['dens'], \
                    data['dens_min'], data['dens_max'], data['dens_orig']
        ax2.set_xlabel('Log Density (cm$^{-3}$)', fontsize=12)
        ax1.semilogx(dens_bins, dens_total_orig, label='Original')
        ax1.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='g', alpha=0.5)
        ax2.axhline(y=1)
        min = [x/y for x,y in zip(dens_total_min, dens_total_orig)]
        max = [x/y for x,y in zip(dens_total_max, dens_total_orig)]
        ax2.fill_between(dens_bins, min, max, color='g', alpha=0.5, label='Range')

    elif type == 'temp':
        with fits.open(file) as hdul:
            data = hdul[1].data
            temp_bins, Te_total_min, Te_total_max, Te_total_orig = data['temp'], \
                    data['Te_min'], data['Te_max'], data['Te_orig']
        ax2.set_xlabel('Log T(K)', fontsize=12)
        ax1.semilogx(temp_bins, Te_total_orig, label='Original')
        ax1.fill_between(temp_bins, Te_total_min, Te_total_max, label='Range', color='g', alpha=0.5)
        ax2.axhline(y=1)
        min = [x/y for x,y in zip(Te_total_min, Te_total_orig)]
        max = [x/y for x,y in zip(Te_total_max, Te_total_orig)]
        ax2.fill_between(temp_bins, min, max, color='g', alpha=0.5, label='Range')

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    ax1.legend(fontsize='xx-small')
    ax2.legend(fontsize='xx-small')
    a = file.split('.')
    a = a[0].split('_')
    name = a[0]+' '+a[1]+ ' blended ratio '
    for number in range(1,10):
        outfile = pathlib.Path(name + str(number) + '.fits')
        if outfile.exists():
            continue
        else:
            plt.savefig(name+str(number)+'.pdf')
            break
    plt.show()

def four_line_diagnostic(Z, z1, Te, dens, process, delta_r, Te_range={}, dens_range={}, type={}, num={}):
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    if Te_range == {}: Te_range = (Te / 10, Te * 10)
    if dens_range == {}: dens_range=(10e12, 10e25)
    if type == {}: type = 'both'
    if num == {}: num = 20

    extras={'process':process, 'delta_r':delta_r, 'transition':{}, 'transition_2':{},
            'npnts':2, 'wavelen':(10,20), 'Te_range':Te_range, 'dens_range':dens_range,
            'corrthresh':10e-5, 'e_signif':0.0}

    if (Z % 2) == 0:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 6), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (6, 1), 'i2': (5, 1)}
    else:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 4), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (4, 1), 'i2': (5, 1)}

    for line, transition in list.items():
        extras.update({'transition':transition})
        inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)

        if process == 'exc':
            new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)
        if process == 'A':
            new_inputs, new_values = variableapec.vary_a(inputs, values, transition)

        table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

        if line == 'r':
            r_diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type=type, num=num)
        elif line == 'f':
            f_diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type=type, num=num)
        elif line == 'i':
            i_diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type=type, num=num)
        elif line == 'i2':
            i2_diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type=type, num=num)

    # write diagnostics
    if type == 'temp':
        temp_bins1, Te_r_orig, Te_r_min, Te_r_max, name1, label1, transition1 = [r_diagnostics.get(k) for k in r_diagnostics]
        temp_bins2, Te_f_orig, Te_f_min, Te_f_max, name2, label2, transition2 = [f_diagnostics.get(k) for k in f_diagnostics]
        temp_bins3, Te_i_orig, Te_i_min, Te_i_max, name3, label3, transition3 = [i_diagnostics.get(k) for k in i_diagnostics]
        temp_bins4, Te_i_orig2, Te_i_min2, Te_i_max2, name4, label4, transition4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temp_bins1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + '_four_line_data_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table.write(element + '_' + str(z1) + 'four_line_data_Te_' + str(number) + '.fits', format='fits')
                break
    elif type == 'dens':
        dens_bins1, dens_r_orig, dens_r_min, dens_r_max, name1, label1, transition1 = [r_diagnostics.get(k) for k in r_diagnostics]
        dens_bins2, dens_f_orig, dens_f_min, dens_f_max, name2, label2, transition2 = [f_diagnostics.get(k) for k in f_diagnostics]
        dens_bins3, dens_i_orig, dens_i_min, dens_i_max, name3, label3, transition3 = [i_diagnostics.get(k) for k in i_diagnostics]
        dens_bins4, dens_i_orig2, dens_i_min2, dens_i_max2, name4, label4, transition4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([dens_bins1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2], names=
                       ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + '_four_line_data_dens_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table.write(element + '_' + str(z1) + 'four_line_data_dens_' + str(number) + '.fits', format='fits')
                break
    elif type == 'both':
        temp_bins1, dens_bins1, Te_r_orig, Te_r_min, Te_r_max, dens_r_orig, \
        dens_r_min, dens_r_max, name1, label1, transition1 = [r_diagnostics.get(k) for k in r_diagnostics]
        temp_bins2, dens_bins2, Te_f_orig, Te_f_min, Te_f_max, dens_f_orig, \
        dens_f_min, dens_f_max, name2, label2, transition2 = [f_diagnostics.get(k) for k in f_diagnostics]
        temp_bins3, dens_bins3, Te_i_orig, Te_i_min, Te_i_max, dens_i_orig, \
        dens_i_min, dens_i_max, name3, label3, transition3 = [i_diagnostics.get(k) for k in i_diagnostics]
        temp_bins4, dens_bins4, Te_i_orig2, Te_i_min2, Te_i_max2, dens_i_orig2, \
        dens_i_min2, dens_i_max2, name4, label4, transition4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temp_bins1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))
        table2 = Table([dens_bins1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2], names=
            ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
             'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + 'four_line_data_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table2.write(element + '_' + str(z1) + 'four_line_data_dens_' + str(number) + '.fits', format='fits')
                table.write(element + '_' + str(z1) + 'four_line_data_Te_' + str(number) + '.fits', format='fits')
                break

def g_ratio(Z, z1, process):
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    if (Z % 2) == 0:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 6), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (6, 1), 'i2': (5, 1)}
    else:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 4), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (4, 1), 'i2': (5, 1)}

    # retrieve diagnostics
    for number in range(20, 0, -1):
        file = pathlib.Path(element + '_' + str(z1) + '_four_line_data_Te_' + str(number) + '.fits')
        if file.exists():
            hdul = fits.open(element + '_' + str(z1) + '_four_line_data_Te_' + str(number) + '.fits')
            data = hdul[1].data
            temp_bins, r_orig, r_min, r_max, f_orig, f_min, f_max, i_orig, i_min, i_max, i2_orig, i2_min, i2_max = data['temp'],
            data['r orig'], data['r min'], data['r max'], data['f orig'], data['f min'], data['f max'], \
            data['i orig'], data['i min'], data['i max'], data['i2 orig'], data['i2 min'], data['i2 max']
            hdul.close()
            break
        else:
            continue

    # temp dependent g ratio
    fig, (ax_1, ax_2) = plt.subplots(nrows=2, sharex=True)
    name = element + ' ' + ion + ' G ratio'
    ax_1.set_xlabel('Log T(K)', fontsize=14)
    ax_1.set_ylabel(name, fontsize=14)

    # do math
    g_min = []
    g_orig = []
    g_max = []
    for rmin, r, rmax, fmin, f, fmax, imin, i, imax, i2min, i2, i2max \
            in zip(r_min, r_orig, r_max, f_min, f_orig, f_max, \
                   i_min, i_orig, i_max, i2_min, i2_orig, i2_max):
        orig = (f + i + i2) / r
        g_orig.append(orig)

        # error propagation for positive dg
        dr, df, di, di2 = (rmax - r), (fmax - f), (imax - i), (i2max - i2)
        r_term = (((-f - i - i2) / r ** 2) ** 2) * (dr ** 2)
        f_term = df ** 2 / r ** 2
        i_term = di ** 2 / r ** 2
        i2_term = di2 ** 2 / r ** 2
        dg = math.sqrt(r_term + f_term + i_term + i2_term)
        g_max.append(orig + dg)

        # error propagation for negative dg
        dr, df, di, di2 = (r - rmin), (f - fmin), (i - imin), (i2 - i2min)
        r_term = (((-f - i - i2) / r ** 2) ** 2) * (dr ** 2)
        f_term = df ** 2 / r ** 2
        i_term = di ** 2 / r ** 2
        i2_term = di2 ** 2 / r ** 2
        dg = math.sqrt(r_term + f_term + i_term + i2_term)
        g_min.append(orig - dg)

    ax_1.semilogx(temp_bins1, g_orig, label='Original')
    ax_1.fill_between(temp_bins, g_min, g_max, alpha=0.5, color='g', \
                      label="Range")
    min = [a / b for a, b in zip(g_min, g_orig)]
    max = [a / b for a, b in zip(g_max, g_orig)]
    ax_2.axhline(y=1, color='black')
    ax_2.fill_between(temp_bins, min, max, color='g', alpha=0.5, label='Range')
    ax_2.set_ylabel('New Ratio/Original', fontsize=14)
    ax_1.legend(fontsize='x-small')
    ax_2.legend(fontsize='x-small')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(element + '_' + ion + '_' + 'G ratio.pdf')

    plt.show()

def r_ratio(Z, z1, process):
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    if (Z % 2) == 0:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 6), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (6, 1), 'i2': (5, 1)}
    else:
        if process == 'exc':
            list = {'r': (1, 7), 'f': (1, 2), 'i': (1, 4), 'i2': (1, 5)}
        elif process == 'A':
            list = {'r': (7, 1), 'f': (2, 1), 'i': (4, 1), 'i2': (5, 1)}

    # retrieve diagnostics
    for number in range(20, 0, -1):
        file = pathlib.Path(element + '_' + str(z1) + '_four_line_data_dens_' + str(number) + '.fits')
        if file.exists():
            hdul = fits.open(element + '_' + str(z1) + '_four_line_data_dens_' + str(number) + '.fits')
            data = hdul[1].data
            dens_bins, r_orig, r_min, r_max, f_orig, f_min, f_max, i_orig, i_min, i_max, i2_orig, i2_min, i2_max = data['dens'],
            data['r orig'], data['r min'], data['r max'], data['f orig'], data['f min'], data['f max'], \
            data['i orig'], data['i min'], data['i max'], data['i2 orig'], data['i2 min'], data['i2 max']
            hdul.close()
            break
        else:
            continue
    # density dependent R ratio
    fig2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    fig2.subplots_adjust(hspace=0)
    name = element + ' ' + ion + ' R ratio'
    ax1.set_ylabel(name, fontsize=14)
    ax1.set_xlabel('Log Density (cm$^{-3}$)', fontsize=14)

    # do math
    R_orig = []
    R_min = []
    R_max = []
    for rmin, r, rmax, fmin, f, fmax, imin, i, imax, i2min, i2, i2max \
            in zip(r_min, r_orig, r_max, f_min, \
                   f_orig, f_max, i_min, i_orig, i_max, \
                   i2_min, i2_orig, i2_max):
        orig = f / (i + i2)
        R_orig.append(orig)

        # error propagation for positive dg
        dr, df, di, di2 = (rmax - r), (fmax - f), (imax - i), (i2max - i2)
        f_term = df ** 2 / (i + i2) ** 2
        i_term = di ** 2 * (-f / (i + i2) ** 2) ** 2
        i2_term = di2 ** 2 * (-f / (i + i2) ** 2) ** 2
        dR = math.sqrt(f_term + i_term + i2_term)
        R_max.append(orig + dR)

        # error propagation for negative dg
        dr, df, di, di2 = (r - rmin), (f - fmin), (i - imin), (i2 - i2min)
        f_term = df ** 2 / (i + i2) ** 2
        i_term = di ** 2 * (-f / (i + i2) ** 2) ** 2
        i2_term = di2 ** 2 * (-f / (i + i2) ** 2) ** 2
        dR = math.sqrt(f_term + i_term + i2_term)
        R_min.append(orig - dR)

    ax1.semilogx(dens_bins, R_orig, label='Original')
    ax1.fill_between(dens_bins1, R_min, R_max, alpha=0.5, color='g', label="Range")
    ax1.legend(fontsize='x-small')
    min = [a / b for a, b in zip(R_min, R_orig)]
    max = [a / b for a, b in zip(R_max, R_orig)]
    ax2.axhline(y=1, color='black')
    ax2.fill_between(dens_bins, min, max, color='g', alpha=0.5, label='Range')
    ax2.set_ylabel('New Ratio/Original', fontsize=14)
    ax2.set_xlabel('Log Density (cm$^{-3}$)')
    ax2.legend(fontsize='x-small')
    fig2.savefig(element + '_' + ion + '_' + 'R ratio.pdf')

    plt.show()

def solve_ionrec(Telist, ionlist, reclist, Z):
    popn = numpy.zeros([len(Telist), Z + 1])
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

        popn[ite, :] = c

    return popn

def find_temp_change(Z, z1, frac, delta_r):
  factor = delta_r
  Telist = numpy.logspace(4, 9, 1251)
  element = pyatomdb.atomic.Ztoelsymb(Z)
  z1_test = z1

  # Change ionization
  varyir = 'i'

  ionlist = numpy.zeros([len(Telist), Z])
  reclist = numpy.zeros([len(Telist), Z])

  # get original rates
  for z1 in range(1, Z + 1):
      iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

      ionlist[:, z1 - 1] = iontmp
      reclist[:, z1 - 1] = rectmp
  eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

  # copy rates to temp variable
  iontmp = ionlist * 1.0
  rectmp = reclist * 1.0

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

  # find peak temperatures
  peak = numpy.max(eqpopn[:, z1_test - 1])
  index = numpy.argmax(eqpopn[:, z1_test-1])
  next_peak = numpy.max(eqpopn[:, z1_test])
  next_index = numpy.argmax(eqpopn[:, z1_test])

  #interpolate
  min_temp = numpy.interp(frac * peak, negpopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_min_temp = numpy.interp(frac * next_peak, negpopn[:next_index, z1_test], Telist[:next_index])
  max_temp = numpy.interp(frac * peak, pospopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_max_temp = numpy.interp(frac * next_peak, pospopn[:next_index, z1_test], Telist[:next_index])

  print("min and max temp:", min_temp, max_temp)
  orig_dt = max_temp - min_temp
  print("diff in temp for", element, str(z1_test - 1) + '+ is', orig_dt)
  orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
  print("orig ion's dex change", orig_dex_change)
  print("next ion's min and max temp:", next_min_temp, next_max_temp)
  next_dt = next_max_temp - next_min_temp
  print("diff in temp for", element, str(z1_test) + '+ is', next_dt)
  next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)
  print("next ion's dex change", next_dex_change)

  dex_changes = {'Uncertainty': '+/-' + str(factor * 100) + '%', 'Ion frac': str(frac) + ' peak',
                 'Unit': 'dex', element + ' ' + str(z1_test - 1) + '+ (ionize)': str(orig_dex_change),
                 element + ' ' + str(z1_test) + '+ (ionize)': str(next_dex_change)}
  temp_changes = {'Uncertainty': ' ', 'Ion frac': ' ',
                  'Unit': 'Kelvin', element + ' ' + str(z1_test - 1) + '+ (ionize)': str(orig_dt),
                  element + ' ' + str(z1_test) + '+ (ionize)': str(next_dt)}

  # Change recombination
  varyir = 'r'

  ionlist = numpy.zeros([len(Telist), Z])
  reclist = numpy.zeros([len(Telist), Z])

  # get the rates
  for z1 in range(1, Z + 1):
      iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

      ionlist[:, z1 - 1] = iontmp
      reclist[:, z1 - 1] = rectmp
  eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

  # copy this rates
  iontmp = ionlist * 1.0
  rectmp = reclist * 1.0

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

  # find peak temperatures
  peak = numpy.max(eqpopn[:, z1_test - 1])
  index = numpy.argmax(eqpopn[:, z1_test-1])
  next_peak = numpy.max(eqpopn[:, z1_test])
  next_index = numpy.argmax(eqpopn[:, z1_test])

  #interpolate
  min_temp = numpy.interp(frac * peak, negpopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_min_temp = numpy.interp(frac * next_peak, negpopn[:next_index, z1_test], Telist[:next_index])
  max_temp = numpy.interp(frac * peak, pospopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_max_temp = numpy.interp(frac * next_peak, pospopn[:next_index, z1_test], Telist[:next_index])

  print("min and max temp:", min_temp, max_temp)
  orig_dt = max_temp - min_temp
  print("diff in temp for", element, str(z1_test - 1) + '+ is', orig_dt)
  orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
  print("orig ion's dex change", orig_dex_change)
  print("next ion's min and max temp:", next_min_temp, next_max_temp)
  next_dt = next_max_temp - next_min_temp
  print("diff in temp for", element, str(z1_test) + '+ is', next_dt)
  next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)
  print("next ion's dex change", next_dex_change)

  # write to csv file
  dex_changes.update({element + ' ' + str(z1_test - 1) + '+ (recomb)': str(orig_dex_change),
                      element + ' ' + str(z1_test) + '+ (recomb)': str(next_dex_change)})
  temp_changes.update({element + ' ' + str(z1_test - 1) + '+ (recomb)': str(orig_dt),
                      element + ' ' + str(z1_test) + '+ (recomb)': str(next_dt)})

  fieldnames = ['Uncertainty', 'Ion frac', 'Unit', element + ' ' + str(z1_test - 1) + '+ (ionize)',
                element + ' ' + str(z1_test) + '+ (ionize)',
                element + ' ' + str(z1_test - 1) + '+ (recomb)', element + ' ' + str(z1_test) + '+ (recomb)']

  import csv
  file = pathlib.Path('temp change ' + element + str(z1_test) + '.csv')
  if file.exists():
    with open('temp change ' + element + str(z1_test) + '.csv', mode='a+') as write_obj:
      dict_writer = csv.DictWriter(write_obj, fieldnames=fieldnames)
      dict_writer.writerow(dex_changes)
      dict_writer.writerow(temp_changes)

  else:
    with open('temp change ' + element + str(z1_test) + '.csv', mode='w') as csv_file:
      dict_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
      dict_writer.writeheader()
      dict_writer.writerow(dex_changes)
      dict_writer.writerow(temp_changes)

def wrapper_find_temp_change(list, frac_list, errors):
  for (Z, z1) in list:
    for frac in frac_list:
      for delta_r in errors:
        print("Finding temp shift at", frac, 'peak for', (Z, z1), 'and', str(delta_r*100), '% uncertainty')
        new_find_dex_change(Z, z1, frac, delta_r, Tebins)

def vary_csd(Z, z1, varyir, delta_r):
    #gets plot of new CSD for new rate. varyir = 'i' for changing ionization, 'r' for recomb
    #delta_r is fractional change, i.e. 10% = 0.10 delta_r

    z1_test = z1
    factor = delta_r
    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    # get the rates
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

    fig = plt.figure()
    # fig.show()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax)
    clist = []
    for z1 in range(1, Z + 2):
        # line,=ax.semilogx(Telist, eqpopn[:,z1-1], label=repr(z1))
        if z1 > 1:
            label = element + ' ' + str(z1 - 1) + '+'
        elif z1 == 1:
            label = element + ' +'
        line, = ax.semilogx(Telist, eqpopn[:, z1 - 1], label=label)
        clist.append(line.get_color())

    # copy this rates
    iontmp = ionlist * 1.0
    rectmp = reclist * 1.0

    # multiply rates by + factor
    if varyir.lower() == 'r':
        rectmp[:, z1_test - 1] *= (1 + factor)
    elif varyir.lower() == 'i':
        iontmp[:, z1_test - 1] *= (1 + factor)
    pospopn = solve_ionrec(Telist, iontmp, rectmp, Z)

    #reset temp variables
    iontmp = ionlist * 1.0
    rectmp = reclist * 1.0

    # multiply rates by - factor
    if varyir.lower() == 'r':
        rectmp[:, z1_test - 1] *= (1 - factor)
    elif varyir.lower() == 'i':
        iontmp[:, z1_test - 1] *= (1 - factor)
    negpopn = solve_ionrec(Telist, iontmp, rectmp, Z)

    for z1 in range(1, Z + 2):
        ax.fill_between(Telist, negpopn[:, z1 - 1], pospopn[:, z1 - 1], color=clist[z1 - 1], alpha=0.5)

    for z1 in range(Z + 1, 0, -1):
        # filter out low popn temperatures
        i = numpy.where(eqpopn[:, z1 - 1] > 1e-5)[0]
        if z1 > 1:
            label = element + ' ' + str(z1 - 1) + '+'
        elif z1 == 1:
            label = element + ' +'

        ax2.fill_between(Telist[i], (eqpopn[i, z1 - 1]) / (negpopn[i, z1 - 1]),
             (eqpopn[i, z1 - 1]) / (pospopn[i, z1 - 1]), color=clist[z1 - 1], label=label) # alpha=0.5)

    ax.legend(fontsize='xx-small')
    ax2.legend(fontsize='xx-small')

    fig.suptitle('CSD from new ' + varyir + ' rate on ' + str(z1_test - 1) + '+')
    ax.set_ylabel("Ion fraction")
    ax2.set_ylabel("Original/New Ion Fraction")
    ax2.set_xlabel("Log T(K)")
    fig.savefig(element+' '+str(z1_test-1)+'+ new CSD '+varyir+' '+str(factor*100)+'%.pdf')

    #do percent changes and slope plots
    fracs = [1e-2, 0.1]

    # find 5 temperatures of interest
    peak = numpy.argmax(eqpopn[:, z1_test - 1])
    increasing_temps = numpy.interp(fracs, eqpopn[:peak, z1_test-1], Telist[:peak])
    decreasing_temps = numpy.interp(fracs, eqpopn[peak:, z1_test-1][::-1], Telist[peak:][::-1])
    temp_bins = numpy.array([increasing_temps[0], increasing_temps[1], Telist[peak], decreasing_temps[1], decreasing_temps[0]])

    # find 5 temperatures of interest for neighbor ion
    peak = numpy.argmax(eqpopn[:, z1_test])
    increasing_temps = numpy.interp(fracs, eqpopn[:peak, z1_test], Telist[:peak])
    decreasing_temps = numpy.interp(fracs, eqpopn[peak:, z1_test][::-1], Telist[peak:][::-1])
    next_temp_bins = numpy.array([increasing_temps[0], increasing_temps[1], Telist[peak], decreasing_temps[1], decreasing_temps[0]])

    fig2, ax2 = plt.subplots(2, 2, figsize=(10, 6))  # percent change plot
    ax2[0, 0].set_xlabel('Fractional Error')
    ax2[0, 0].set_ylabel('New ' + element + str(z1_test - 1) + '+ /Original')
    ax2[1, 0].set_xlabel('T (K)')
    ax2[1, 0].set_ylabel('Slope')
    ax2[0, 1].set_ylabel('New ' + element + str(z1_test) + '+ /Original')
    ax2[1, 1].set_ylabel('Slope')
    ax2[1, 1].set_xlabel('T (K)')
    ax2[0, 1].set_xlabel('Fractional Error')

    percent_names = [' (1%)', ' (10%)', ' (peak)', ' (10%)', ' (1%)']
    slopes = []
    new_errors = [-0.20, -0.15, -0.10, +0.10, +0.15, +0.20]
    frac = numpy.zeros([len(temp_bins), len(new_errors)])
    next_frac = numpy.zeros([len(temp_bins), len(new_errors)])

    ionlist = numpy.zeros([len(temp_bins), Z])
    reclist = numpy.zeros([len(temp_bins), Z])
    next_ionlist = numpy.zeros([len(next_temp_bins), Z])
    next_reclist = numpy.zeros([len(next_temp_bins), Z])

    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(temp_bins, False, Z=Z, z1=z1, extrap=True)
        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

        next_iontmp, next_rectmp = pyatomdb.atomdb.get_ionrec_rate(next_temp_bins, False, Z=Z, z1=z1, extrap=True)
        next_ionlist[:, z1 - 1] = next_iontmp
        next_reclist[:, z1 - 1] = next_rectmp

    eqpopn = solve_ionrec(temp_bins, ionlist, reclist, Z)
    next_eqpopn = solve_ionrec(next_temp_bins, ionlist, reclist, Z)

    for i in range(len(new_errors)):
        delta_r = new_errors[i]

        # copy rates
        iontmp = ionlist * 1.0
        rectmp = reclist * 1.0

        # multiply rates by + factor
        if varyir.lower() == 'r':
            rectmp[:, z1_test - 1] *= (1 + delta_r)
        elif varyir.lower() == 'i':
            iontmp[:, z1_test - 1] *= (1 + delta_r)
        newpopn = solve_ionrec(temp_bins, iontmp, rectmp, Z)
        next_newpopn = solve_ionrec(next_temp_bins, iontmp, rectmp, Z)

        frac[:, i] = (newpopn[:, z1_test - 1]) / (eqpopn[:, z1_test - 1])
        next_frac[:, 1] = (next_newpopn[:, z1_test]) / (next_eqpopn[:, z1_test])

    slopes = []
    for i, percent in zip(range(len(temp_bins)), percent_names):
        Te = temp_bins[i]
        y = frac[i, :]
        label = numpy.format_float_scientific(Te, exp_digits=1, precision=1) + ' K' + percent
        ax2[0, 0].plot(new_errors, y, label=label, marker='.')
        slope, intercept = numpy.polyfit(new_errors, y, 1)
        slopes.append(slope)
        ax2[1, 0].semilogx(Te, slope, label=label, marker='.', zorder=2)
    ax2[1, 0].semilogx(temp_bins, slopes, color='black', linestyle='-', zorder=1)

    next_slopes = []
    for i, percent in zip(range(len(next_temp_bins)), percent_names):
        Te = next_temp_bins[i]
        y = next_frac[i, :]
        label = numpy.format_float_scientific(Te, exp_digits=1, precision=1) + ' K' + percent
        ax2[0, 1].plot(new_errors, y, label=label, marker='.')
        slope, intercept = numpy.polyfit(new_errors, y, 1)
        next_slopes.append(slope)
        ax2[1, 1].semilogx(Te, slope, label=label, marker='.', zorder=2)
    ax2[1, 1].semilogx(next_temp_bins, next_slopes, color='black', linestyle='-', zorder=1)

    for ax in (ax2[0, 0], ax2[0, 1], ax2[1, 0], ax2[1, 1]):
        ax.legend(fontsize='xx-small')

    plt.tight_layout()
    fig2.savefig(element+' '+str(z1_test-1)+'+ CSD %change '+varyir+ ' '+str(factor*100)+'%.pdf')

    plt.show()

def monte_carlo_csd(Z, max_error, runs, makefiles=False, plot=True):
    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    clist = get_cmap(Z+2)

    #set up plot
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlabel('Log T(K)', fontsize=12)
    plt.ylabel('Ion Fraction', fontsize=12)

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

    mc_popn = numpy.zeros([runs, Z+1, len(Telist)])
    random_ion = numpy.zeros([len(Telist), Z])
    random_rec = numpy.zeros([len(Telist), Z])
    for run in range(runs):
        lower, upper = -2*max_error, 2*max_error
        mu, sigma = 0, max_error
        random1 = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
        random2 = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
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

        for z1 in range(1,Z+2):
            mc_popn[run,z1-1,:] = newpopn[:,z1-1]

    if plot==True:
        #find mean and standard dev of each column
        for z1 in range(1, Z+2):
            median=[]
            middle, min, max = [], [], []
            for i in range(len(Telist)):
                pop = mc_popn[:,z1-1, i]
                pop_list = numpy.sort(pop)
                median.append(numpy.median(pop_list))
                min_index = int(0.16*runs-1)
                max_index = int(0.84*runs-1)
                min.append(pop_list[min_index])
                max.append(pop_list[max_index])

            if (z1 == 26) | (z1 == 17):
                label = element + ' ' + str(z1 - 1) + '+'
                plt.semilogx(Telist, median, label=label, color=clist(z1 - 1), linestyle='-')
            else:
                plt.semilogx(Telist, median, color=clist(z1-1), linestyle='-')

            plt.fill_between(Telist, min, max, color = clist(z1-1), alpha=0.4)

        plt.title('CSD with Monte Carlo error between +/-'+str(max_error*100)+'%')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol = 9, fontsize='xx-small')
        plt.savefig(element+' '+str(max_error*100)+'% Monte Carlo CSD.pdf')
        plt.close('all')

def wrapper_monte_carlo_csd(list, errors, runs, makefiles=False, plot=True):
    "List is [], errors is []"

    for Z in list:
        for max_error in errors:
            monte_carlo_csd(Z, max_error, runs, makefiles, plot)

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def extremize_csd(Z, max_error, runs, makefiles=False, plot=True):
    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    clist = get_cmap(Z + 2)

    # set up plot
    plt.figure()
    ax = plt.subplot(111)
    plt.xlabel('Log T(K)', fontsize=12)
    plt.ylabel('Ion Fraction', fontsize=12)

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

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
                label = element + ' +'
            else:
                label = element + ' ' + str(z1 - 1) + '+'
            plt.semilogx(Telist, median, label=label, color=clist(z1 - 1), linestyle='-')
            plt.fill_between(Telist, min, max, color=clist(z1 - 1), alpha=0.4)

        plt.title("Monte Carlo CSD, ionize + " + str(max_error * 100) + "%, recomb - " + str(max_error*100) + "%")
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol = 9, fontsize='xx-small')
        plt.savefig(element + ' i+ r- ' + str(max_error*100) +'% Monte Carlo CSD.pdf')
        plt.close('all')

    # case 2: ionize - error, recomb + error
    plt.figure()
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
                # pop = []
                # for run in range(runs):
                #    pop.append(mc_popn[run, z1-1, i])
                pop = mc_popn[:, z1 - 1, i]
                pop_list = numpy.sort(pop)
                median.append(numpy.median(pop_list))
                min_index = int(0.16 * runs - 1)
                max_index = int(0.84 * runs - 1)
                min.append(pop_list[min_index])
                max.append(pop_list[max_index])
            # print("length of Telist, means, devs", len(Telist), len(means), len(devs))
            if z1 == 1:
                label = element + ' +'
            else:
                label = element + ' ' + str(z1 - 1) + '+'
            plt.semilogx(Telist, median, label=label, color=clist(z1 - 1), linestyle='-')
            # min = [a - b for a, b in zip(means, devs)]
            # max = [a + b for a, b in zip(means, devs)]
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

    #ionlist, reclist = vary_rates(Z, Telist, ion_deltar=ion_deltar, rec_deltar=rec_deltar)

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

    fname = 'varied_eigen%s_%s.fits' % (atomic.Ztoelsymb(Z).lower(), filesuffix)
    hdulist.writeto(fname, checksum=True, clobber=True)

def update_eigen(Z, filename, session):
    """
    Updates eigenvector file data for element Z from file filename

    """
    d = pyatomdb.pyfits.open(filename)

    session.spectra.datacache['data']['misc']['EIGEN'][Z] = d
    session.spectra.datacache['datasums']['misc']['EIGEN'][Z] = d[1].header['DATASUM']

def new_line_emissivity(Z, z1, up, lo, Te, files):
  #gets plots of new emissivities from eigen files made from new ionization/recombination rates

  #files is a list [] of all eigen files to use for getting new emissivities
  #Te can either be a single temp in K (3e6), 'peak' for Te of peak ion fraction,
  #a list [1e6,3e6] or numpy array = numpy.logspace(4,9,50)

  tau = 1e13

  if (type(Te) != list) and (type(Te) != numpy.ndarray):
    print("Running for one temperature.")
    varied_emiss = numpy.zeros([1, len(files)+1])
    if Te == 'peak':
      Te = ratios.find_peak_Te(Z, z1, unit='keV')
    else:
      # convert input Te in K to keV
      Te = Te / 11604525.0061657
    s = pyatomdb.spectrum.NEISession(elements=[Z])
    orig = s.return_line_emissivity(Te, tau, Z, z1, up, lo, init_pop=Te)
    varied_emiss[0,0] = orig['epsilon']
    for f in range(1, len(files)+1):
      file = files[f-1]
      update_eigen(Z, file, s)
      ret = s.return_line_emissivity(Te, tau, Z, z1, up, lo, init_pop=Te)
      varied_emiss[0,f] = ret['epsilon']

    # now plot emissivity of upper lower over runs
    plt.figure()
    plt.title('Line emissivity '+'('+str(up)+'->'+str(lo)+') from MC random error')
    plt.xlabel('Runs')
    plt.ylabel('Emissivity')
    plt.plot(varied_emiss[0, :], marker='o')
    plt.show()

  elif (type(Te) == list) or (type(Te) == numpy.ndarray):
    varied_emiss = numpy.zeros([len(Te), len(files)+1])
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    ax2.set_xlabel('Log T (K)')
    ax1.set_ylabel('Emissivity')
    ax2.set_ylabel('% change in emissivity')
    print("Got a list of temps.")
    for i in range(len(Te)):
      temp_Te = Te[i]
      if temp_Te == 'peak':
        temp_Te = ratios.find_peak_Te(Z, z1, unit='keV')
      else:
        # convert input Te in K to keV
        temp_Te = temp_Te / 11604525.0061657
      print("for temp:", temp_Te, "we have:")
      s = pyatomdb.spectrum.NEISession(elements=[Z])
      orig = s.return_line_emissivity(temp_Te, tau, Z, z1, up, lo, init_pop=temp_Te)
      varied_emiss[i,0] = orig['epsilon']
      for f in range(1, len(files)+1):
        file = files[f-1]
        update_eigen(Z, file, s)
        ret = s.return_line_emissivity(temp_Te, tau, Z, z1, up, lo, init_pop=temp_Te)
        varied_emiss[i,f] = ret['epsilon']

    # for each run, plot emiss vs. temp
    for f in range(len(files)+1):
      percent = [((a-b)/b) for a,b in zip(varied_emiss[:, f], varied_emiss[:,0])]
      if f == 0:
        label='Original'
        ax1.semilogx(Te, varied_emiss[:, f], label=label, color='black', zorder=-1)
      else:
        label=str(f)
        ax1.semilogx(Te, varied_emiss[:, f], label=label)
        ax2.semilogx(Te, percent, label=label)

    fig.suptitle('Line emissivity '+'('+str(up)+'->'+str(lo)+') from MC random error')
    ax1.legend(fontsize='xx-small')
    ax2.legend(fontsize='xx-small')
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.show()
    plt.close('all')

def new_line_ratio(Z, z1, up1, lo1, up2, lo2, Te, files):
  #gets plot of new line ratio from eigen files made from new ionization and recombination rates
  #Te can be a single temp in K (3e6), 'peak' if want Te of peak ion fraction,
  # a list [1e6,3e6], or array = numpy.logspace(4,9,50)
  import ratios
  tau = 1e13

  if (type(Te) != list) and (type(Te) != numpy.ndarray):
    print("Running for one temperature.")
    line_ratio = numpy.zeros([1, len(files)+1])
    if Te == 'peak':
      Te = ratios.find_peak_Te(Z, z1, unit='keV')
    else:
      # convert input Te in K to keV
      Te = Te / 11604525.0061657
    s = pyatomdb.spectrum.NEISession(elements=[Z])
    orig1 = s.return_line_emissivity(Te, tau, Z, z1, up1, lo1, init_pop=Te)
    orig2 = s.return_line_emissivity(Te, tau, Z, z1, up2, lo2, init_pop=Te)
    line_ratio[0,0] = orig1['epsilon']/orig2['epsilon']
    for f in range(1, len(files)+1):
      file = files[f-1]
      update_eigen(Z, file, s)
      ret1 = s.return_line_emissivity(Te, tau, Z, z1, up1, lo1, init_pop=Te)
      ret2 = s.return_line_emissivity(Te, tau, Z, z1, up2, lo2, init_pop=Te)
      line_ratio[0,f] = ret1['epsilon']/ret2['epsilon']

    # now plot emissivity of upper lower over runs
    plt.figure()
    plt.title('Line ratio '+'('+str(up1)+'->'+str(lo1)+'/'+str(up2)+'->'+str(lo2)+') from MC shifted error')
    plt.xlabel('Runs')
    plt.ylabel('Line Ratio')
    plt.plot(line_ratio[0, :], marker='o')
    plt.show()

  elif (type(Te) == list) or (type(Te) == numpy.ndarray):
    line_ratio = numpy.zeros([len(Te), len(files)+1])
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    ax2.set_xlabel('Log T (K)')
    ax1.set_ylabel('Line Ratio')
    ax2.set_ylabel('% change in line ratio')
    print("Got a list of temps.")
    for i in range(len(Te)):
      temp_Te = Te[i]
      if temp_Te == 'peak':
        temp_Te = ratios.find_peak_Te(Z, z1, unit='keV')
      else:
        # convert input Te in K to keV
        temp_Te = temp_Te / 11604525.0061657
      print("for temp:", temp_Te, "we have:")
      s = pyatomdb.spectrum.NEISession(elements=[Z])
      orig1 = s.return_line_emissivity(temp_Te, tau, Z, z1, up1, lo1, init_pop=temp_Te)
      orig2 = s.return_line_emissivity(temp_Te, tau, Z, z1, up2, lo2, init_pop=temp_Te)
      line_ratio[i,0] = orig1['epsilon']/orig2['epsilon']
      print("Orig is", line_ratio[i,0])
      for f in range(1, len(files)+1):
        file = files[f-1]
        update_eigen(Z, file, s)
        ret1 = s.return_line_emissivity(temp_Te, tau, Z, z1, up1, lo1, init_pop=temp_Te)
        ret2 = s.return_line_emissivity(temp_Te, tau, Z, z1, up2, lo2, init_pop=temp_Te)
        line_ratio[i,f] = ret1['epsilon']/ret2['epsilon']
        print("New is", line_ratio[i,f])

    # for each run, plot emiss vs. temp
    for f in range(len(files)+1):
      percent = [((a-b)/b) for a,b in zip(line_ratio[:, f], line_ratio[:,0])]
      if f == 0:
        label='Original'
        ax1.semilogx(Te, line_ratio[:, f], label=label, color='black', zorder=-1)
      else:
        label=str(f)
        ax1.semilogx(Te, line_ratio[:, f], label=label)
        ax2.semilogx(Te, percent, label=label)

    fig.suptitle('Line ratio '+'('+str(up1)+'->'+str(lo1)+'/'+str(up2)+'->'+str(lo2)+') from MC shifted error')
    ax1.legend(fontsize='xx-small', ncol=2)
    ax2.legend(fontsize='xx-small', loc='upper right', ncol=2)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0, top=0.9)
    plt.show()
    plt.close('all')

def get_new_ionreclist(Z, Telist, ion_deltar={}, rec_deltar={}):
    #vary all ions of element Z, vary either ion rate, rec rate or both

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    #get original rates
    for z1 in range(1, Z + 1):
        iontmp, rectmp = atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True,
                                                settings=False)
        if ion_deltar != {}:
            #change ion rate
            iontmp = [x + (x * (ion_deltar / 100)) for x in iontmp]

        if rec_deltar != {}:
            #change rec rate
            rectmp = [x + (x * (rec_deltar / 100)) for x in rectmp]

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    return ionlist, reclist

def get_partial_deriv(Z, z1, errors, dens={}):
    element = pyatomdb.atomic.Ztoelsymb(Z)
    print("Calculating partial derivatives for", element, str(z1))

    #find temperature where ion fraction peaks
    peak_Te = find_peak_Te(Z, z1)

    if (Z, z1) == (26,17):
        trans_list = [(2,1), (3,1), (5,1), (17, 1), (23,1), (27,1)]
    if (Z, z1) == (26, 19):
        trans_list = [(53,1), (68,1), (71,1), (74, 1), (76, 4)]
    if Z-z1 == 1:
        trans_list = [(2,1), (5,1), (6,1), (7,1), (13,1)]
        print("Varying He-like lines.")
    if Z==z1:
        print("Varying H-like lines.")
        trans_list=[(3,1), (4,1)]

    #define critical densities for ions
    if dens=={}: dens=1
    elif dens=='critical':
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
            print("Critical density found:", critical_dens)
            dens = critical_dens
            print("Will vary rates at density", dens)
        elif critical_dens == {}:
            print("No critical density found, will do for low density")
            dens = 1

    uppers = [x[0] for x in trans_list]
    lowers = [x[1] for x in trans_list]

    #get orig emissivities so use any random transition for set up
    inputs, values = set_up(Z, z1, peak_Te, dens)
    table = values.get('table')

    # set up output table
    t = Table([uppers, lowers], names=('Upper', 'Lower'))
    wavelen, orig = [],[]
    for x in t:
        for y in table:
            if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                wavelen.append(y['Lambda'])
                orig.append(y['Epsilon_orig'])
    new1 = Column(name='Orig', data=orig)
    new2 = Column(name='Lambda', data=wavelen)
    t.add_columns([new2, new1])
    print("Table is set up", '\n', t)

    t2 = Table([uppers, lowers, wavelen, orig], names=('Upper', 'Lower', 'Lambda', 'Orig'))
    t3 = Table(names=('Upper', 'Lower', 'Lambda', 'Orig'))
    t4 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))
    t5 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))
    t6 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))

    counter1=1.1
    for delta_r in errors:
        counter2 = counter1
        for transition in trans_list:
            print("Varying the transition", transition, 'by', str(delta_r*100)+'%')

            #vary exc rate
            extras={'process':exc, 'delta_r': delta_r, 'transition':transition[::-1], 'transition_2':[],
                'npnts':2, 'wavelen':(10,20), 'Te_range':(peak_Te/10, peak_Te*10), 'dens_range':(1,10e16),
                    'corrthresh':0.0, 'e_signif':0.0}

            inputs, values, exc_transition = set_up(Z, z1, peak_Te, dens, extras=extras)
            new_inputs, new_values = vary_exc(inputs, values, exc_transition)
            q_max, q_min = new_values.get('q_max'), new_values.get('q_min')
            table, new_table, inputs, results = get_tables(new_inputs, new_values)

            partial_deriv, frac_E = [], []
            for x in t:
                for y in table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        epsilon_max = y['Epsilon_' + str(npnts)]
                        print("epsilon max =", epsilon_max, 'epsilon min=', y['Epsilon_1'], 'orig', y['Epsilon_orig'])
                        deriv = (epsilon_max - y['Epsilon_1'])/(q_max - q_min)
                        if abs(deriv) < 0.0001: deriv = 0.0
                        partial_deriv.append(deriv)
                        frac = (epsilon_max - y['Epsilon_1']) / y['Epsilon_orig']
                        if abs(frac) < 0.0001: frac = 0.0
                        frac_E.append(frac)

            name_1 = 'dE/dR exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
            name_2 = 'dE/E exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t.add_columns([deriv])
            t2.add_columns([frac])

            #remove rows of lines already checked separately because tables are ordered differently
            idx1=[]
            idx2 = []
            for x in trans_list:
                for i in range(len(new_table)):
                    if (table[i]['Upper'], table[i]['Lower']) == x:
                        idx1.append(i)
                for i in range(len(new_table)):
                    if (new_table[i]['Upper'], new_table[i]['Lower']) == x:
                        idx2.append(i)
            table.remove_rows(idx1)
            new_table.remove_rows(idx2)

            #now search through tables for other lines affected and add to tables
            for x in new_table:
                #interesting lines
                if 0.02 <= x['dE/dR'] <= 0.98:
                    print("Found interesting line with dE/dR = ", x['dE/dR'])
                    if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                        row = [x['Upper'], x['Lower'], x['Lambda'], x['Epsilon_orig']]+['0']*(len(t3.colnames)-4)
                        t3.add_row(row)
                        t4.add_row(row)
                        print("Adding this line, \n", t3, '\n', t4)
                #linear change lines
                elif 0.98 <= x['dE/dR'] <= 1.02:
                    print("Found linear change line with dE/dR = ", x['dE/dR'])
                    if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                        t5.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']]+['0']*(len(t5.colnames)-4))
                        t6.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))

            #get partial derivs for interesting lines
            partial_deriv, frac_E = [], []
            for x in t3:
                for y in new_table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        deriv, frac = y['dE/dR'], y['|dE/dE_orig|']
                        if abs(deriv) < 10e-5: deriv = 0.0
                        if abs(frac) < 10e-5: frac = 0.0
                        partial_deriv.append(deriv)
                        frac_E.append(frac)
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t3.add_columns([deriv])
            t4.add_columns([frac])
            print(t3, '\n', t4)

            #get partial derivs for lines that change linearly
            partial_deriv, frac_E = [], []
            for x in t5:
                for y in new_table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        deriv, frac = y['dE/dR'], y['|dE/dE_orig|']
                        if abs(deriv) < 10e-5: deriv = 0.0
                        if abs(frac) < 10e-5: frac = 0.0
                        partial_deriv.append(deriv)
                        frac_E.append(frac)
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t5.add_columns([deriv])
            t6.add_columns([frac])


            #vary A value
            extras.update({'process':A, 'transition':transition})
            inputs, values, transition = set_up(Z, z1, peak_Te, dens, extras=extras)
            new_inputs, new_values = vary_a(inputs, values, transition)
            q_max, q_min = new_values.get('q_max'), new_values.get('q_min')
            table, new_table, inputs, results = get_tables(new_inputs, new_values)


            partial_deriv, frac_E = [], []
            for x in t:
                for y in table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        epsilon_max = y['Epsilon_' + str(npnts)]
                        print("epsilon max =",epsilon_max, 'epsilon min=', y['Epsilon_1'],'orig', y['Epsilon_orig'])
                        deriv = (epsilon_max - y['Epsilon_1']) / (q_max - q_min)
                        if abs(deriv) < 0.0001: deriv = 0.0
                        partial_deriv.append(deriv)
                        frac = (epsilon_max - y['Epsilon_1'])/y['Epsilon_orig']
                        if abs(frac) < 0.0001: frac = 0.0
                        frac_E.append(frac)

            name_1 = 'dE/dR A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
            name_2 = 'dE/E A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t.add_columns([deriv])
            t2.add_columns([frac])

            #remove existing lines
            idx1 = []
            idx2 = []
            for x in trans_list:
                for i in range(len(new_table)):
                    if (table[i]['Upper'], table[i]['Lower']) == x:
                        idx1.append(i)
                for i in range(len(new_table)):
                    if (new_table[i]['Upper'], new_table[i]['Lower']) == x:
                        idx2.append(i)
            table.remove_rows(idx1)
            new_table.remove_rows(idx2)

            # now search through tables for other lines affected and add to tables
            for x in new_table:
                # interesting lines
                if 0.02 <= x['|dE/dE_orig|'] <= 0.98:
                    print("Found additional line impacted - ", x['|dE/dE_orig|'])
                    if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                        t3.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t3.colnames) - 4))
                        t4.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t3.colnames) - 4))
                        print("Adding this line, \n", t3, '\n', t4)
                # linear change lines
                elif 0.98 <= x['|dE/dE_orig|'] <= 1.02:
                    print("Found additional line impacted - ", x['|dE/dE_orig|'])
                    if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                        t5.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))
                        t6.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))

            # get partial derivs for interesting lines
            partial_deriv, frac_E = [], []
            for x in t3:
                for y in new_table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        deriv, frac = y['dE/dR'], y['|dE/dE_orig|']
                        if abs(deriv) < 10e-5: deriv = 0.0
                        if abs(frac) < 10e-5: frac = 0.0
                        partial_deriv.append(deriv)
                        frac_E.append(frac)
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t3.add_columns([deriv])
            t4.add_columns([frac])
            print(t3, '\n', t4)

            # get partial derivs for lines that change linearly
            partial_deriv, frac_E = [], []
            for x in t5:
                for y in new_table:
                    if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                        deriv, frac = y['dE/dR'], y['|dE/dE_orig|']
                        if abs(deriv) < 10e-5: deriv = 0.0
                        if abs(frac) < 10e-5: frac = 0.0
                        partial_deriv.append(deriv)
                        frac_E.append(frac)
            deriv = Column(name=name_1, data=partial_deriv)
            frac = Column(name=name_2, data=frac_E)
            t5.add_columns([deriv])
            t6.add_columns([frac])

            counter2+=0.1
        counter1+=1

    #write csv file
    element = pyatomdb.atomic.Ztoelsymb(Z)
    for number in range(1,10):
        file = pathlib.Path('varied ' + element + str(z1) + ' (' + str(number)+').csv')
        if file.exists():
            continue
        else:
            with open('varied ' + element + str(z1) + ' (' + str(number)+').csv', mode='w') as csv_file:
                fieldnames = t.colnames
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                for row in t:
                    data = {}
                    for col in fieldnames:
                        data.update({col: row[col]})
                    writer.writerow(data)

                emptyrow = {}
                for col in fieldnames:
                    emptyrow.update({col: ' '})
                writer.writerow(emptyrow)

                fieldnames2=t2.colnames
                print("Fieldnames1: ", fieldnames, '\nFieldnames2', fieldnames2)
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames2)
                writer.writeheader()

                for row in t2:
                    data = {}
                    for col in fieldnames2:
                        data.update({col: row[col]})
                    writer.writerow(data)

                emptyrow = {}
                for col in fieldnames2:
                    emptyrow.update({col: ' '})
                writer.writerow(emptyrow)
                writer.writerow({'Upper':'Additional lines impacted:'})

                #sort interesting lines table by emissivity
                t3.sort('Orig', reverse=True)
                t4.sort('Orig', reverse=True)
                t5.sort('Orig', reverse=True)
                t6.sort('Orig', reverse=True)

                fieldnames3=t3.colnames
                writer=csv.DictWriter(csv_file,fieldnames=fieldnames3)
                writer.writeheader()

                for row in t3:
                    data = {}
                    for col in fieldnames3:
                        data.update({col: row[col]})
                    writer.writerow(data)
                emptyrow = {}
                for col in fieldnames3:
                    emptyrow.update({col: ' '})
                writer.writerow(emptyrow)

                fieldnames4 = t4.colnames
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames4)
                writer.writeheader()

                for row in t4:
                    data = {}
                    for col in fieldnames4:
                        data.update({col: row[col]})
                    writer.writerow(data)

            with open('varied ' + element + str(z1) + ' linear lines (' + str(number) + ').csv', mode='w') as csv_file:
                fieldnames = t5.colnames
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                for row in t5:
                    data = {}
                    for col in fieldnames:
                        data.update({col: row[col]})
                    writer.writerow(data)

                emptyrow = {}
                for col in fieldnames:
                    emptyrow.update({col: ' '})
                writer.writerow(emptyrow)

                fieldnames2 = t6.colnames
                writer2 = csv.DictWriter(csv_file, fieldnames=fieldnames2)
                writer2.writeheader()

                for row in t6:
                    data = {}
                    for col in fieldnames2:
                        data.update({col: row[col]})
                    writer2.writerow(data)

            break

def wrapper_get_partial_deriv(list, errors, dens={}):
    for x in list:
        for error in errors:
            for ne in dens:
                if type(x[0]) == int:
                    Z, z1 = x[0], x[1]
                    get_partial_deriv(Z, z1, [error], ne)
                elif type(x[0]) == str:
                    Z, z1 = pyatomdb.atomic.elsymb_to_Z(x[0]), x[1]
                    get_partial_deriv(Z, z1, [error], ne)

def find_peak_Te(Z, z1, unit='K'):
    # find temperature where ion fraction peaks
    Te_array = numpy.logspace(4, 9, 51)
    ion_frac = []
    ionlist = numpy.zeros([len(Te_array), Z])
    reclist = numpy.zeros([len(Te_array), Z])

    for temp_z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Te_array, False, Z=Z, z1=temp_z1, extrap=True,
                                                         settings=False)

        ionlist[:, temp_z1 - 1] = iontmp
        reclist[:, temp_z1 - 1] = rectmp
    eqpopn = solve_ionrec(Te_array, ionlist, reclist, Z)
    ion_frac = eqpopn[:, z1-1]

    index = numpy.argmax(ion_frac)
    peak_Te = Te_array[index]

    if unit == 'keV':
        peak_Te = peak_Te/11604525.0061657

    return peak_Te