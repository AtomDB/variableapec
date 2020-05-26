<<<<<<< HEAD
"""
This module contains methods for looking at emission and line ratios
as a function of varying atomic data from the AtomDB files. Requires
PyAtomdB and Python 3."""

# Keri Heuer
# Version 2.0, March 25, 2020

import matplotlib.pyplot as plt, matplotlib.ticker as mtick, scipy.stats as stats, \
pyatomdb, numpy, pickle, pathlib, csv, os, errno, hashlib, requests, urllib.request, \
urllib.parse, urllib.error, time, subprocess, shutil, wget, glob, datetime, ftplib, \
pathlib, collections, operator, requests, matplotlib.pylab as pylab, glob, math
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

def ionize(Z, z1, Te, dens, in_range, pop_fraction, datacache={}, newsettings={}, newsettings2={}):
    z1_drv, d = z1-1, datacache
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1_drv, Te, dens, do_la= True, \
                            do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)
    
    lvdat = pyatomdb.atomdb.get_data(Z, z1_drv, 'LV', datacache=d)
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

    if newsettings == {}:
        ion_levpop = pyatomdb.apec.calc_ioniz_popn(drv_lev_pop*pop_fraction[z1_drv-1], Z, z1, z1_drv, Te, dens, datacache=d)
    elif newsettings != {}:
        which_transition, x, old_rate = [newsettings.get(k) for k in newsettings]
        ion_levpop = calc_ioniz_popn(which_transition, x, old_rate, drv_lev_pop*pop_fraction[z1_drv-1], Z, z1, z1_drv, Te, dens, datacache=d)

    ion_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    ion_linelist['lambda'] = in_range['WAVELEN']
    ion_linelist['epsilon'] = [x['EINSTEIN_A'] * ion_levpop[x['UPPER_LEV'] - 1] for x in in_range]
    return ion_linelist
            
def recombine(Z, z1, Te, dens, in_range, pop_fraction, datacache={}, newsettings={}, newsettings2={}):
    z1_drv, d = z1+1, datacache
    if z1 < Z:
        init, final, rates = pyatomdb.apec.gather_rates(Z, z1_drv, Te, dens, do_la=True, \
                                                        do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=d)

        lvdat = pyatomdb.atomdb.get_data(Z, z1_drv, 'LV', datacache=d)
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

    if newsettings == {}:
        recomb_levpop = pyatomdb.apec.calc_recomb_popn(drv_lev_pop * pop_fraction[z1_drv - 1], Z, z1, z1_drv, Te, dens,
                                                       drlevrates=0, rrlevrates=0,datacache=d)
    elif newsettings != {}:
        which_transition, x, old_rate = [newsettings.get(k) for k in newsettings]
        recomb_levpop = calc_recomb_popn(which_transition, x, old_rate, drv_lev_pop*pop_fraction[z1_drv-1], Z, z1,
                                        z1_drv, Te, dens,drlevrates=0, rrlevrates=0,datacache=d)

    recomb_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    recomb_linelist['lambda'] = in_range['WAVELEN']
    recomb_linelist['epsilon'] = [x['EINSTEIN_A']*recomb_levpop[x['UPPER_LEV']-1] for x in in_range]
    return recomb_linelist

def set_up(Z, z1, Te, dens, extras={}):
    
    """ Uses inputs from check_sensitivity(), see routine for variable definitions as they are same throughout.
=======
#Version 0.0, Oct 30, 2019
#Keri Heuer


import matplotlib.pyplot as plt
import pyatomdb, numpy, pickle
from astropy.table import Table, Column
from astropy.io import ascii
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from decimal import Decimal

def set_up(Z, z1, Te, dens, process, delta_r, transition, transition_2, all_levels, \
           nlev, npnts, show, wavelen, corrthresh, e_signif):
    
    """ Uses inputs from check_sensitivity(), see that subroutine for variable definitions.
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    Set_up() gets original rates, radiative collision matrix, level populations, line lists, 
    and intensities and returns these in a dictionary called values. Creates sensitivity tables.
    
    Returns inputs in dictionary called inputs, returns the variable transition."""
<<<<<<< HEAD

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
=======
    
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= True, \
                            do_ec=True, do_ir=True, do_pc=True, do_ai=True)
        
    if all_levels == False:
        print("Only solving for levels 1 to", nlev)
        index = numpy.where((init < nlev) & (final < nlev))[0]
        init = init[index]
        final = final[index]
        rates = rates[index]    #sum all rates (the q's)
    else:
        lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV')
        lvdat = lvdat[1].data
        nlev = len(lvdat)
                   
    matrix = numpy.zeros((nlev,nlev))
    B = numpy.zeros(nlev)
    #populate matrix by summing rates for all processes
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    for i in range(len(init)):
        x = final[i]
        y = init[i]
        matrix[x][y] += rates[i]
<<<<<<< HEAD

    # find fraction of each ion in plasma to multiply epsilons by
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

    #set up and solve CR matrix for level populations
    matrix[0][:] = 1.0
    B[0]=1.0
    lev_pop = numpy.linalg.solve(matrix,B)
    lev_pop *= pop_fraction[z1-1]

    #convert level populations into line lists & intensities for excitation only
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
    in_range = ladat[1].data
        
    linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    linelist['lambda'] = in_range['WAVELEN']
    linelist['epsilon'] = [x['EINSTEIN_A']*lev_pop[x['UPPER_LEV']-1] for x in in_range]

    #set up complete line list (emiss only due to excitation at this point)
    full_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
    full_linelist['lambda'] = linelist['lambda']
    full_linelist['epsilon'] = linelist['epsilon']

    #now add emissivity from ionization and recombination to excitation linelist (depending on z1)
    if z1 == 1: #skip ionization
        recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
        full_linelist['epsilon'] += recomb_emiss['epsilon']
    elif z1 == Z+1: #skip recombination
        ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
        full_linelist['epsilon'] += ion_emiss['epsilon']
    else: #do both
        recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
        ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
        full_linelist['epsilon'] += recomb_emiss['epsilon']
        full_linelist['epsilon'] += ion_emiss['epsilon']

    #set up sensitivity & partial derivatives tables
    table = Table([full_linelist['lambda'], in_range['UPPER_LEV'], in_range['LOWER_LEV'], full_linelist['epsilon']], \
        names=('Lambda', 'Upper', 'Lower', 'Epsilon_orig'))
    new_table = Table([full_linelist['lambda'], in_range['UPPER_LEV'], in_range['LOWER_LEV']], names=('Lambda', 'Upper', 'Lower'))

    #save variables
    values = {'matrix': matrix, 'B': B, 'in_range': in_range, 'linelist': linelist, 'table': table, 'new_table': new_table}
    if extras == {}:
        inputs = {'Z': Z, 'z1': z1, 'Te': Te, 'dens': dens}
        return inputs, values
    else:
        process, delta_r, transition, transition_2, npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [extras.get(k) for k in extras]
        inputs = {'Z': Z, 'z1': z1, 'Te': Te, 'dens': dens, 'process': process, 'delta_r': delta_r,
        'transition': transition, 'transition_2': transition_2, 'npnts': npnts, 'wavelen': wavelen,
        'Te_range': Te_range, 'dens_range': dens_range, 'corrthresh': corrthresh, 'e_signif': e_signif}
        return inputs, values, transition

def vary_a(inputs, values, which_transition):
=======
        
    #set up and solve matrix for level populations
    matrix[0][:] = 1.0
    B[0]=1.0
    lev_pop = numpy.linalg.solve(matrix,B)
    
    #convert level populations into line lists & intensities
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA')
    ladat = ladat[1].data
    
    if all_levels == False:
        in_range = ladat[ladat['UPPER_LEV'] < nlev+1] 
    else:
        in_range = ladat
        
    linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype')) 
    
    upper_level= numpy.zeros(len(in_range))
    lower_level = numpy.zeros(len(in_range))
    
    for i in range(len(in_range)):
        linelist['lambda'][i] = in_range['WAVELEN'][i]
        pop_level = in_range['UPPER_LEV'][i]
        linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*lev_pop[pop_level-1]
        upper_level[i] = in_range['UPPER_LEV'][i]
        lower_level[i] = in_range['LOWER_LEV'][i] 
    
    #set up sensitivity & partial derivatives tables
    table = Table([linelist['lambda'], upper_level, lower_level, linelist['epsilon']], \
        names=('Lambda', 'Upper', 'Lower', 'Epsilon_orig'))
    new_table = Table([linelist['lambda'], upper_level, lower_level], names=('Lambda', 'Upper', 'Lower'))

    #save variables 
    inputs = {}
    inputs.update( {'Z': Z, 'z1': z1, 'Te': Te, 'dens': dens, 'process': process, 'delta_r': delta_r, \
                    'transition': transition, 'transition_2': transition_2, \
                    'all_levels': all_levels, 'nlev': nlev, 'npnts': npnts, \
                    'show': show, 'wavelen': wavelen, 'corrthresh': corrthresh, 'e_signif': e_signif} )
    values = {}
    values.update( {'matrix': matrix, 'B': B, 'in_range': in_range, 'linelist': linelist, 'table': table, 'new_table': new_table} )
    return inputs, values, transition

def vary_a(inputs, values, transition):
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the A value for by the specified delta_r
    in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
<<<<<<< HEAD
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
=======
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("****************************************************")
<<<<<<< HEAD
    print("in vary_a for " + str(which_transition) + ", calculating rate for Te=%e and dens=%e" %(Te, dens))
    print("****************************************************")
    
    initial_lev, final_lev = which_transition[0], which_transition[1]
    
    old_a = in_range['EINSTEIN_A'][(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)][0]
    a_index = numpy.where([(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)])[1][0]
    table['Epsilon_orig'].unit='A'

    # find fraction of each ion in plasma and multiple level pops by z1 frac
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

    if old_a == 0:old_a = 1e-40
    
    #vary A values
    min_a, max_a = 1-delta_r, 1+delta_r
    new_a = [min_a*old_a, max_a*old_a]
    print("old A value =", old_a)
    print("new A values =", new_a)
    q_max, q_min = new_a[-1], new_a[0]
=======
    print("in vary_a, calculating rate for Te=%e and dens=%e"%(Te, dens))
    print("****************************************************")
    
    initial_lev = transition[0]
    final_lev = transition[1]
    
    old_a = in_range['EINSTEIN_A'][(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)][0]
    a_index = numpy.where([(in_range['UPPER_LEV'] == initial_lev) & (in_range['LOWER_LEV'] == final_lev)])[1][0]
    
    table['Epsilon_orig'].unit='A'
    
    #vary A values
    min_a = 1-delta_r
    max_a = 1+delta_r
    new_a = numpy.linspace(min_a*old_a, max_a*old_a, npnts)
    print("old A value =", old_a)
    print("new A values =", new_a)
    q_max = new_a[-1]
    q_min = new_a[0]
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
        
    index=1
    for x in new_a:
        #update LA data for specified transition
        in_range['EINSTEIN_A'][a_index] = x
<<<<<<< HEAD

        #get new CR matrix and resolve level pops with new A
        frac = str(round(x/old_a,2))
        new_matrix = matrix.copy()
        new_matrix[final_lev-1, initial_lev-1] += (x-old_a)   #off diagonal term
        new_matrix[initial_lev-1, initial_lev-1] -= (x-old_a)   #diagonal term
        new_matrix[0][:] = 1.0

        #find new level populations and get new epsilon values from excitation
        new_lev_pop = numpy.linalg.solve(new_matrix,B)
        new_lev_pop *= pop_fraction[z1-1]
        new_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
        new_linelist['lambda'] = in_range['WAVELEN']
        new_linelist['epsilon'] = [y['EINSTEIN_A'] * new_lev_pop[y['UPPER_LEV']-1] for y in in_range]

        # now add emissivity from ionization and recombination to full linelist (depending on z1)
        nset = {'which_transition': which_transition, 'x': x, 'old_rate': old_a}
        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            new_linelist['epsilon'] += ion_emiss['epsilon']
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
            new_linelist['epsilon'] += ion_emiss['epsilon']
=======
        
        #get new matrix and resolve level pops with new A
        frac = str(round(x/old_a,2))
        new_matrix = matrix.copy()
        new_linelist=linelist.copy()
        new_matrix[final_lev-1, initial_lev-1] += (x-old_a)   #off diagonal term
        new_matrix[initial_lev-1, initial_lev-1] -= (x-old_a)   #diagonal term
        new_matrix[0][:] = 1.0
        
        #find new level populations and get new epsilon values
        new_lev_pop = numpy.linalg.solve(new_matrix,B)
            
        for i in range(len(in_range)):
            new_linelist['lambda'][i] = in_range['WAVELEN'][i] #wavelength will be same
            pop_level = in_range['UPPER_LEV'][i]
            new_linelist['epsilon'][i] = in_range['EINSTEIN_A'][i]*new_lev_pop[pop_level-1]
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
        
        #update sensitivity table 
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' A')
        table.add_columns([new_col])
<<<<<<< HEAD
        index+=1
        
    values = {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, 'q_max': q_max, 'q_min': q_min}
    return inputs, values
    
def vary_exc(inputs, values, which_transition):
=======
            
        index+=1
        
    values = {}
    values.update( {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, \
                    'q_max': q_max, 'q_min': q_min} )
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    percent=str(delta_r*100)
    print("\nFor", element, ion + ", changed inputs:", "A value for level", \
          str(initial_lev)+"->"+str(final_lev))
    
    return inputs, values
    
def vary_exc(inputs, values, transition):
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the excitation rate for by the specified
    delta_r in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
<<<<<<< HEAD
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("***********************************************************************************")
    print("in vary_exc for " + str(which_transition) + ", calculating rate for Te=%e and dens=%e" %(Te, dens))
    print("***********************************************************************************")
    
    exc_init, exc_final, exc_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= False, \
                        do_ec=True, do_ir=False, do_pc=False, do_ai=False, datacache=d)
    
    initial_lev, final_lev = which_transition[0], which_transition[1]

    try:
        for a,b,c in zip(exc_init, exc_final, exc_rates):
            if (a,b) == (which_transition[0]-1, which_transition[1]-1):
                old_rate = c
    except:
        which_transition = which_transition[::-1]
        for a,b,c in zip(exc_init, exc_final, exc_rates):
            if (a,b) == (which_transition[0]-1, which_transition[1]-1): old_rate = c
    try:
        if old_rate == 0: old_rate = 1e-40
    except UnboundLocalError:
        print("Could not find transition", which_transition, " - please check input transition levels")

    table['Epsilon_orig'].unit='orig rate'

    # find fraction of each ion in plasma and multiple level pops by z1 frac
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
        
    #vary rate
    min_rate, max_rate = 1-delta_r, 1+delta_r
    new_rate = [min_rate*old_rate, max_rate*old_rate]
    print("old exc rate =", old_rate)
    print("new exc rates =", new_rate)
    q_max, q_min = new_rate[-1], new_rate[0]

=======
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("****************************************************")
    print("in vary_exc, calculating rate for Te=%e and dens=%e"%(Te, dens))
    print("****************************************************")
    
    exc_init, exc_final, exc_rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la= False, \
                        do_ec=True, do_ir=False, do_pc=False, do_ai=False)
    
    initial_lev=transition[0]
    final_lev=transition[1]
    
    old_rate_index = numpy.where((exc_final == final_lev-1) & (exc_init == initial_lev-1))[0]
    old_rate = exc_rates[old_rate_index][0]
    
    table['Epsilon_orig'].unit='orig rate'
        
    #vary rate
    min_rate = 1-delta_r
    max_rate = 1+delta_r
    new_rate = numpy.linspace(min_rate*old_rate, max_rate*old_rate, npnts)
    print("old rate =", old_rate)
    print("new rates =", new_rate)
    q_max = new_rate[-1]
    q_min = new_rate[0]
        
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    index=1
    for x in new_rate:
        #loop through varied rates, get new matrix and resolve level pops
        frac = str(round(x/old_rate,2))
        new_matrix = matrix.copy()
<<<<<<< HEAD
        new_matrix[final_lev-1, initial_lev-1] += (x-old_rate)   #off diagonal term
        new_matrix[initial_lev-1, initial_lev-1] -= (x-old_rate)   #diagonal term
        new_matrix[0][:] = 1.0

        # find new level populations and get new epsilon values from excitation
        new_lev_pop = numpy.linalg.solve(new_matrix, B)
        new_lev_pop *= pop_fraction[z1 - 1]
        new_linelist = numpy.zeros(len(in_range), dtype=pyatomdb.apec.generate_datatypes('linetype'))
        new_linelist['lambda'] = in_range['WAVELEN']
        new_linelist['epsilon'] = [y['EINSTEIN_A'] * new_lev_pop[y['UPPER_LEV'] - 1] for y in in_range]

        # now add emissivity from ionization and recombination to full linelist (depending on z1)
        nset={'which_transition':which_transition, 'x':x, 'old_rate':old_rate}
        if z1 == 1:  # skip ionization
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction,newsettings=nset)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
        elif z1 == Z + 1:  # skip recombination
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            new_linelist['epsilon'] += ion_emiss['epsilon']
        else:  # do both
            recomb_emiss = recombine(Z, z1, Te, dens, in_range, pop_fraction, newsettings=nset)
            ion_emiss = ionize(Z, z1, Te, dens, in_range, pop_fraction,newsettings=nset)
            new_linelist['epsilon'] += recomb_emiss['epsilon']
            new_linelist['epsilon'] += ion_emiss['epsilon']

        #update sensitivity table
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' rate')
        table.add_columns([new_col])
        index +=1
    
    values = {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, 'q_max': q_max, 'q_min': q_min}
=======
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
        
        #update sensitivity table 
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' rate')
        table.add_columns([new_col])
        
        index+=1
    
    values = {}
    values.update( {'table': table, 'new_table': new_table, 'new_linelist': new_linelist, \
                    'q_max': q_max, 'q_min': q_min} )

    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    percent=str(delta_r*100)
    print("\nFor", element, ion + ", changed inputs:", "excitation rate from level", \
          str(initial_lev)+"->"+str(final_lev)) 
    
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    return inputs, values

def get_tables(inputs, values):
    
    """ Inputs and values are dictionaries outputted by either vary_exc() or vary_a().
    
    Calculates dE/dR as the difference in min and max emissivities divided by the difference in changed rates,
    dE/dE_orig as the difference in min and max emissivities divided by the original emissivity, and sorts
<<<<<<< HEAD
    the sensitivity table by dE/E.
=======
    the sensitivity table by dE/dE_orig. 
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    If the variables corrthresh and e_signif were specified by the user, it will filter the
    sensitivity table for values with dE/dE_orig greater than corrthresh and/or for lines
    with a dE/dR greater than e_signif.
    
    Prints the sensitivity table.
    
    Returns table and new_table and the dictionaries inputs and results."""
    
<<<<<<< HEAD
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    npnts, wavelen, Te_range, dens_range, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    table, new_table, new_linelist, q_max, q_min = [values.get(k) for k in values]
    min, max = table['Epsilon_1'], table['Epsilon_'+str(npnts)]
    
    #ratio_col = Column(name = 'Epsilon_max/Epsilon_min', data = (epsilon_max/table['Epsilon_1']))
    #table.add_columns([ratio_col])
    # table['Epsilon_max/Epsilon_min'] = max/min
    # table['Epsilon_max/Epsilon_min'].unit = None
        
    #add partial derivative, dE/dR
    # delta_epsilon = Column(name = 'dE/dR', data = (epsilon_max-table['Epsilon_1'])/(q_max-q_min))
    # new_table.add_columns([delta_epsilon])
    new_table['dE/dR'] = (max-min)/(q_max-q_min)
    new_table['dE/dR'].unit = None
    
    # orig_eps = Column(name = 'Epsilon_orig', data = table['Epsilon_orig'])
    # new_table.add_columns([orig_eps])
    new_table['Epsilon_orig'] = table['Epsilon_orig']
    
    #add "correlation factor" dE/dE_orig
    # epsilon_corr = Column(name = 'dE/E', data = (epsilon_max-table['Epsilon_1'])/table['Epsilon_orig'])
    # abs_epsilon_corr = Column(name = '|dE/E|', data = [abs(val) for val in epsilon_corr])
    # new_table.add_columns([abs_epsilon_corr])
    new_table['|dE/E|'] = [abs(val) for val in max-min/table['Epsilon_orig']]
    try:
      new_table.sort('|dE/E|', reverse=True)
    except TypeError:
      new_table.sort('|dE/E|')
      new_table = new_table[::-1]
    
    #apply filters
    if corrthresh != 0.0:     #only show lines whose "epsilon correlation" >= than specified value
        new_table = new_table[new_table['|dE/E|'] >= corrthresh]
    elif e_signif != 0.0:    #only show lines with partial epsilon/partial rate derivative is >= specified value
        new_table = new_table[new_table['dE/dR'] >= eps_der]

    results = {'inputs': inputs, 'wavelength': numpy.array(new_table['Lambda']),
                    'upper': numpy.array(new_table['Upper']), 'lower': numpy.array(new_table['Lower']), \
                    'dE/dR': numpy.array(new_table['dE/dR']), 'epsilon_orig': numpy.array(new_table['Epsilon_orig']),\
                    '|dE/E|': numpy.array(new_table['|dE/E|']), 'min_eps': numpy.array(min),'max_eps': max}
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
=======
    Z, z1, Te, dens, process, delta_r, transition, transition_2, all_levels, \
    nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
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
    new_table.sort('|dE/dE_orig|', reverse=True)
    
    #apply filters
    if corrthresh != 0.0:     #only show lines whose "epsilon correlation" >= than specified value
        new_table = new_table[new_table['|dE/dE_orig|'] >= corrthresh]
    elif e_signif != 0.0:    #only show lines with partial epsilon/partial rate derivative is >= specified value
        new_table = new_table[new_table['dE/dR'] >= eps_der]

    #print("Lines affected at Te="+str(Te) + ", dens="+str(dens))
    #print(table)
    #print(new_table)
    
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
                    'max_eps': epsilon_max })
    
    return table, new_table, inputs, results
    
def plot_sensitivity(inputs, new_table):
    
    """ Inputs is a dictionary outputted by get_tables() and new_table is the sensitivity
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    table outputted by get_tables().
    
    Plots the lines affected by the parameter(s) changed, including the original transition.
    Currently plotting for wavelengths between 10-20 Angstroms (as set by mask1 and mask2 below).""" 
    
<<<<<<< HEAD
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
    cutoff_data = new_table[new_table['|dE/E|'] >= corrthresh]
    if wavelen != {}:
        cutoff_data = cutoff_data[(wavelen[0] < cutoff_data['Lambda']) & (cutoff_data['Lambda'] < wavelen[1])]
    
    #plot wavelength vs. % emissivity change
    axs[0,1].semilogy(cutoff_data['Lambda'], cutoff_data['|dE/E|'], linestyle='', marker='o', label=transition)
=======
    Z, z1, Te, dens, process, delta_r, transition, transition_2, all_levels, \
    nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(which_transition[0])+'-'+str(which_transition[1])
    name=element+'_'+ion+'_'+process+'_'+line+'_'

    #set up plot of "sensitive epsilons"
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('Lambda ($\AA$)', fontsize=18) 
    ax1.set_ylabel('% Emissivity Change', fontsize=18)
    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    temp = '%.E' % Decimal(Te) +'K'
    percentage = str(delta_r*100)+'%'
    if process == 'A':
        text='A value varied by $\Delta$ $\pm$' + percentage + ', ' + temp
    elif process == 'exc':
        text='Direct excitation rate varied by $\Delta$ $\pm$' + percentage + ', ' + temp  
    anchored_text = AnchoredText(text, loc='center', frameon=False)
    ax1.add_artist(anchored_text)
    
    #filter data for significance 
    cutoff_data = new_table[new_table['|dE/dE_orig|'] >= corrthresh]
    mask1 = wavelen[0] < cutoff_data['Lambda'] 
    cutoff_data = cutoff_data[mask1]
    mask2 = cutoff_data['Lambda'] < wavelen[1]
    cutoff_data = cutoff_data[mask2]
    
    #plot wavelength vs. % emissivity change
    ax1.semilogy(cutoff_data['Lambda'], cutoff_data['|dE/dE_orig|'], linestyle='', marker='o', label=transition)
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    #label each point w/ transition
    transition_labels=[]
    for x in cutoff_data:
        if process == 'A':
            transition_name = '{0:g}'.format(x['Upper'])+'->'+'{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name = '{0:g}'.format(x['Lower'])+'->'+'{0:g}'.format(x['Upper'])
        transition_labels.append(transition_name)
<<<<<<< HEAD
    for (x,y,label) in zip(numpy.array(cutoff_data['Lambda']),numpy.array(cutoff_data['|dE/E|']),transition_labels):
        axs[0,1].annotate(label, xy=(x,y))

    axs[0,1].legend(fontsize='x-small')

def run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range={}, dens_range={}, num={}, plot=False):
    
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

    print("Running", type, "diagnostics")
    start = timer()

    extras = {'process': vary, 'delta_r': delta_r, 'transition': (up, lo), 'transition_2': [], 'npnts': 2,
              'wavelen': (10, 20), 'Te_range': Te_range, 'dens_range': dens_range, 'corrthresh': 10e-5, 'e_signif': 0.0}

    if temperature == True:
        temp_bins = list(numpy.geomspace(Te_range[0], Te_range[1], num=temp_num))  # 20
        Te_eps_orig, Te_eps_min, Te_eps_max =[], [], []
        counter=0
        for temp_Te in temp_bins:
            if vary == 'A':
                Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, extras=extras)
                Te_new_inputs, Te_new_values = vary_a(Te_inputs, Te_values, (up, lo))
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, extras=extras)
                Te_new_inputs, Te_new_values = vary_exc(Te_inputs, Te_values, (lo, up))
            Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
            for x in Te_table:
                if (x['Upper'], x['Lower']) == (up, lo):
                    if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                    if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                    if x[-1] == 0.0: x[-1] = 1e-40
                    Te_eps_orig.append(x['Epsilon_orig'])
                    Te_eps_min.append(x['Epsilon_1'])
                    Te_eps_max.append(x[-1])
            counter += 1
            print(str(temp_num - counter), 'temperatures left\n')

    if density == True:
        dens_bins = list(numpy.geomspace(dens_range[0], dens_range[1], num=dens_num))
        dens_eps_orig, dens_eps_min, dens_eps_max =[],[],[]
        counter=0
        for temp_dens in dens_bins:
            if vary == 'A':
                dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, extras=extras)
                dens_new_inputs, dens_new_values = vary_a(dens_inputs, dens_values, (up, lo))
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, extras=extras)
                dens_new_inputs, dens_new_values = vary_exc(dens_inputs, dens_values, (lo, up))
            dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)
            for x in dens_table:
                if (x['Upper'], x['Lower']) == (up, lo):
                    for i in [x['Epsilon_orig'], x['Epsilon_1'], x[-1]]:
                        if i==0: i=1e-40
                    # if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                    # if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                    # if x[-1] == 0.0: x[-1] = 1e-40
                    dens_eps_orig.append(x['Epsilon_orig'])
                    dens_eps_min.append(x['Epsilon_1'])
                    dens_eps_max.append(x[-1])
            counter += 1
            print(str(dens_num-counter), 'densities left\n')
    print("Took", (timer() - start)/60, "minutes")

    if type == 'temp':
        line_diagnostics = {'type': 'temp', 'temps': list(temp_bins),'orig': list(Te_eps_orig),
                                 'min': list(Te_eps_min),'max': list(Te_eps_max)}
    elif type == 'dens':
        line_diagnostics = {'type': 'dens', 'dens': list(dens_bins),'orig': list(dens_eps_orig),
                                 'min': list(dens_eps_min),'max': list(dens_eps_max)}
    elif type == 'both':
        line_diagnostics = {'type': 'both', 'temps': list(temp_bins), 'dens': list(dens_bins),
                                 'Te_orig': list(Te_eps_orig), 'Te_min': list(Te_eps_min),
                                 'Te_max': list(Te_eps_max), 'dens_orig': list(dens_eps_orig),
                                 'dens_min': list(dens_eps_min), 'dens_max': list(dens_eps_max)}
    if plot == True:
        plot_line_diagnostics((up, lo), vary, delta_r, line_diagnostics)
    return line_diagnostics
    
def plot_line_diagnostics(transition, vary, delta_r, line_diagnostics):
=======
    for (x,y,label) in zip(numpy.array(cutoff_data['Lambda']),numpy.array(cutoff_data['|dE/dE_orig|']),transition_labels):
        ax1.annotate(label, xy=(x,y))
    
    if process == 'A':
        plt.savefig(name+'epsilon_sensitivity.pdf')
    elif process == 'exc':
        plt.savefig(name+'epsilon_sensitivity.pdf')

    plt.tight_layout()
    plt.legend()
    plt.show()
    
def wrapper_plot_sensitivity(inputs, new_table1, new_table2): #needs to be finished
    
    """ Results is a dictionary.
    
    Plots lines affected by the changed parameter(s) including the original transition,
    for the multiple transitions specified, in a different color for each transition.
    
    By default, will plot for wavelengths between 10-20 Angstroms, unless a wavelength
    range is specified by wavelen=()."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, all_levels, \
    nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    
    print("******")
    print("wavelen is ", wavelen)
    print("corrthres is", corrthresh)
    print("***")
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    temp = '%.E' % Decimal(Te) +'K'
    percentage = '{0:g}'.format(delta_r*100)+'%'
    density = 'dens='+str(dens)
    line1=str(transition[0])+'-'+str(transition[1])
    line2=str(transition_2[0])+'-'+str(transition_2[1])
    name=element+'_'+ion+'_'+process+'_'+line1+'_'+line2+'_'

    #set up plot of "sensitive epsilons"
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('Lambda ($\AA$)', fontsize=18) 
    ax1.set_ylabel('% Emissivity Change', fontsize=18)
    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    temp = '%.E' % Decimal(Te) +'K'
    percentage = str(delta_r*100)+'%'
    if process == 'A':
        text='A value varied by $\Delta$ $\pm$' + percentage + ', ' + temp
    elif process == 'exc':
        text='Direct excitation rate varied by $\Delta$ $\pm$' + percentage + ', ' + temp  
    anchored_text = AnchoredText(text, loc='center', frameon=False)
    ax1.add_artist(anchored_text)
    
    #filter data for significance 
    cutoff_data1 = new_table1[new_table1['|dE/dE_orig|'] >= corrthresh]
    mask1 = wavelen[0] < cutoff_data1['Lambda'] 
    cutoff_data1 = cutoff_data1[mask1]
    mask2 = cutoff_data1['Lambda'] < wavelen[1]
    cutoff_data1 = cutoff_data1[mask2]
    
    cutoff_data2 = new_table2[new_table2['|dE/dE_orig|'] >= corrthresh]
    mask1 = wavelen[0] < cutoff_data2['Lambda'] 
    cutoff_data2 = cutoff_data2[mask1]
    mask2 = cutoff_data2['Lambda'] < wavelen[1]
    cutoff_data2 = cutoff_data2[mask2]
    
    #plot wavelength vs. % emissivity change for both transitions
    ax1.semilogy(cutoff_data1['Lambda'], cutoff_data1['|dE/dE_orig|'], linestyle='', marker='o', label=transition)
    print("done")
    ax1.semilogy(cutoff_data2['Lambda'], cutoff_data2['|dE/dE_orig|'], linestyle='', marker='o', label=transition_2)
    print("done with 2")
    #label each point with transition
    transition_labels_1=[]
    transition_labels_2=[]
    for x in cutoff_data1:
        if process == 'A':
            transition_name_1 = '{0:g}'.format(x['Upper'])+'->'+'{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name_1 = '{0:g}'.format(x['Lower'])+'->'+'{0:g}'.format(x['Upper'])
        transition_labels_1.append(transition_name_1)
    for y in cutoff_data2:
        if process == 'A':
            transition_name_2 = '{0:g}'.format(y['Upper'])+'->'+'{0:g}'.format(y['Lower'])
        if process == 'exc':
            transition_name_2 = '{0:g}'.format(y['Lower'])+'->'+'{0:g}'.format(y['Upper'])
        transition_labels_2.append(transition_name_2)
    for (x,y,label) in zip(numpy.array(cutoff_data1['Lambda']),numpy.array(cutoff_data1['|dE/dE_orig|']),transition_labels_1):
        ax1.annotate(label, xy=(x,y))
    for (x,y,label) in zip(numpy.array(cutoff_data2['Lambda']),numpy.array(cutoff_data2['|dE/dE_orig|']),transition_labels_2):
        ax1.annotate(label, xy=(x,y))
    
    plt.savefig(name+'epsilon_sensitivity.pdf')

    plt.tight_layout()
    plt.legend()
    plt.show()
        
def run_line_diagnostics(table, inputs, values, which_transition):
    
    """ Table is table of transitions and emissivity values. Inputs and values are dictionaries and
    which_transition is a tuple specifying which transition to run line diagnostics for.
    
    Varies temperature and density separately and recalculates emissivitiy values. 
    
    Returns the dictionary line_diagnostics of arrays of the temperature and density bins used,
    the original, min, and max emissivities, the name of element and ion, label of what process varied
    and by how much, and which_transition."""
    
    #read in inputs and values from set_up() and table from get_table()
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]

    print("////////////////////////////////////////////////////////////////")
    print("Now running line diagnostics for transition:", which_transition)
    print("////////////////////////////////////////////////////////////////")
    
    table = values.get('table')
    new_table = values.get('new_table')
    new_linelist = values.get('new_linelist')
    q_max = values.get('new_linelist')
    q_min = values.get('q_min')
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line=str(which_transition[0])+'-'+str(which_transition[1])
    name=element+'_'+ion+'_'+process+'_'+line+'_'
    percentage=str(delta_r*100)+'%'
    
    line_diagnostics = {}
    
    #vary temperature and recalculate emissivities 
    print("\nChanging temperature now.\n")
    if isinstance(Te, list) != True:
        t_bins = numpy.geomspace(Te/10, Te*10, num=20)
    elif isinstance(Te, list) == True:
        t_bins = numpy.geomspace(Te[0], Te[2], num=20)

    temp_bins=[]
    for x in t_bins:
        temp_bins.append(x)
    
    Te_eps_orig=[]
    Te_eps_min=[]
    Te_eps_max=[]
    for temp_Te in temp_bins:
        print("Temperature is:", temp_Te)
        Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, process, delta_r, which_transition, transition_2, all_levels, \
           nlev, npnts, show, wavelen, corrthresh, e_signif)
        if process == 'A':
            Te_new_inputs, Te_new_values = vary_a(Te_inputs, Te_values, which_transition)
            Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
            line_index = numpy.where((Te_table['Upper'] == transition[0]) & (Te_table['Lower'] == transition[1]))[0]
        elif process == 'exc':
            Te_new_inputs, Te_new_values = vary_exc(Te_inputs, Te_values, which_transition)
            Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
            line_index = numpy.where((Te_table['Upper'] == transition[1]) & (Te_table['Lower'] == transition[0]))[0]
              
        Te_table = Te_table[line_index]
                       
        Te_eps_orig.append(Te_table['Epsilon_orig'][0])
        Te_eps_min.append(Te_table['Epsilon_1'][0])
        Te_eps_max.append(Te_table['Epsilon_'+str(npnts)][0])
    
    #vary density and recalculate emissivities
    print("\nChanging density now.\n")
    if isinstance(dens, int) == True:
        bins = numpy.geomspace(10e0, 10e16, num=8)
    elif isinstance(dens, list) == True:
        bins = numpy.geomspace(dens[0], dens[2], num=8)

    dens_bins=[]
    for x in bins:
        dens_bins.append(x)
        
    
    dens_eps_orig=[]
    dens_eps_min=[]
    dens_eps_max=[]
    for temp_dens in dens_bins:
        print("Density is:", temp_dens)
        dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, process, delta_r, which_transition, transition_2, all_levels, \
           nlev, npnts, show, wavelen, corrthresh, e_signif)
        if process == 'A':
            dens_new_inputs, dens_new_values = vary_a(dens_inputs, dens_values, which_transition)
            dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)
            line_index = numpy.where((dens_table['Upper'] == transition[0]) & (dens_table['Lower'] == transition[1]))[0]
        elif process == 'exc':
            dens_new_inputs, dens_new_values = vary_exc(dens_inputs, dens_values, which_transition)
            dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)
            line_index = numpy.where((dens_table['Upper'] == transition[1]) & (dens_table['Lower'] == transition[0]))[0]

        dens_table = dens_table[line_index]
                       
        dens_eps_orig.append(dens_table['Epsilon_orig'][0])
        dens_eps_min.append(dens_table['Epsilon_1'][0])
        dens_eps_max.append(dens_table['Epsilon_'+str(npnts)][0])
    
    if process == 'A':
        label = 'Range of emissivities due to $\Delta $\pm$' + percentage + ' in A value'
    elif process == 'exc':
        label = 'Range of emissivities due to $\Delta $\pm$' + percentage + ' in direct excitation rate'
        
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

    file = open('line_diagnostics', 'wb')
    pickle.dump(line_diagnostics, file)
    file.close()

    return line_diagnostics
    
def plot_line_diagnostics(inputs, line_diagnostics):
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    """ Line_diagnostics is a dictionary outputted by run_line_diagnostics. It holds arrays of
    the temperature and density bins, the original, minimum, and maximum emissivities from varying
    temperature and then density, the name of the element and ion, label of which process varied
    and by how much, and the transition.
    
    Plots emissivity as a function of temperature and density for the specified single transition."""

<<<<<<< HEAD
    print("Plotting line diagnostics now.")
    type = line_diagnostics.get('type')
    text = '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1])

    if type == 'temp':  # plot emissivity versus temperature
        type, temps, Te_orig, Te_min, Te_max = [line_diagnostics.get(k) for k in line_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax.set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        #ax.set_xlabel('Temperature in K')
        ax.semilogx(temps, Te_orig, label='Original', color='b')
        ax.fill_between(temps, Te_min, Te_max, alpha=0.5, color='b', label="Range")
        min, max = [x / y for x, y in zip(Te_min, Te_orig)], [x / y for x, y in zip(Te_max, Te_orig)]
        ax2.axhline(y=1, color='k')
        ax2.fill_between(temps, min, max, label='Range', color='b', alpha=0.5)
        ax2.set_xlabel('Temperature in K')
        ax2.set_ylabel('New/Original Emissivity')
        ax.legend(fontsize='xx-small')
        ax2.legend(fontsize='xx-small')
    elif type == 'dens':
        type, dens, dens_orig, dens_min, dens_max = [line_diagnostics.get(k) for k in line_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax.set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        #ax2.set_xlabel('Density in cm$^{-3}$')
        ax.semilogx(dens, dens_orig, label='Original', color='b')
        ax.fill_between(dens, dens_min, dens_max, alpha=0.5, color='b', label="Range")
        min, max = [x / y for x, y in zip(dens_min, dens_orig)], [x / y for x, y in zip(dens_max, dens_orig)]
        ax2.axhline(y=1, color='k')
        ax2.fill_between(dens, min, max, label='Range', color='b', alpha=0.5)
        ax2.set_xlabel('Density in cm$^{-3}$')
        ax2.set_ylabel('New/Original Emissivity')
        ax.legend(fontsize='xx-small')
        ax2.legend(fontsize='xx-small')
    elif type == 'both':
        type, temps, dens, Te_orig, Te_min, Te_max, dens_orig, dens_min, dens_max = [line_diagnostics.get(k) for k in line_diagnostics]
        fig, ax = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col')
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax[0, 0].set_ylabel('Emissivity in \n$\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$')
        ax[1,0].set_xlabel('Temperature in K')
        ax[1, 1].set_xlabel('Density in cm$^{-3}$')
        ax[1, 0].set_ylabel('New/Original Emissivity')
        ax[0, 1].semilogx(dens, dens_orig, label='Original', color='b')
        ax[0, 0].semilogx(temps, Te_orig, label='Original', color='b')
        ax[0, 1].fill_between(dens, dens_min, dens_max, alpha=0.5, color='b', label="Range")
        ax[0, 0].fill_between(temps, Te_min, Te_max, color='b', alpha=0.5, label="Range")
        ax[1, 0].axhline(y=1, color='k')
        ax[1, 1].axhline(y=1, color='k')
        min, max = [x / y for x, y in zip(dens_min, dens_orig)], [x / y for x, y in zip(dens_max, dens_orig)]
        ax[1, 1].fill_between(dens, min, max, label='Range', color='b', alpha=0.5)
        min, max = [x / y for x, y in zip(Te_min, Te_orig)], [x / y for x, y in zip(Te_max, Te_orig)]
        ax[1, 0].fill_between(temps, min, max, label='Range', color='b', alpha=0.5)
        for ax in [ax[0,0], ax[0,1], ax[1,1], ax[1,0]]:
            ax.legend(fontsize='xx-small')

    plt.tight_layout()
    plt.subplots_adjust(top=0.86, hspace=0)
    plt.savefig('test line diagnostics.pdf')
    plt.show()
    plt.close('all')

def run_line_ratio_diagnostics(transition1, transition2, line_diagnostics_1, line_diagnostics_2, vary, delta_r, plot=False):
=======
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    file = open('line_diagnostics', 'rb')
    line_diagnostics = pickle.load(file)
    file.close()

    temp_bins, dens_bins, Te_eps_orig, Te_eps_min, Te_eps_max, dens_eps_orig, \
            dens_eps_min, dens_eps_max, name, label, transition = [line_diagnostics.get(k) for k in line_diagnostics]

    temp = '%.E' % Decimal(Te) +'K'
    percentage = '{0:g}'.format(delta_r*100)+'%'
    density = 'dens='+str(dens)
    if process == 'A':
        orig_text='{0:g}'.format(transition[0])+'->'+'{0:g}'.format(transition[1]) + \
        ' A value $\Delta$ $\pm$' + percentage 
    elif process == 'exc':
        orig_text='{0:g}'.format(transition[0])+'->'+'{0:g}'.format(transition[1]) + \
        ' Direct excitation rate $\Delta$ $\pm$' + percentage   
    
    
    #plot emissivity versus temperature
    text = orig_text + ', ' + density
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    ax1.set_xlabel('Temperature in K', fontsize=18)
    ax1.set_ylabel('Emissivity in $\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=18)
    anchored_text = AnchoredText(text, loc='lower right', frameon=False)
    ax1.add_artist(anchored_text)
    ax1.loglog(temp_bins, Te_eps_orig, label='Original emissivity')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_temp.pdf')
    ax1.fill_between(temp_bins, Te_eps_min, Te_eps_max, alpha=0.5, color='g', \
                     label="Range of emissivities")
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_temp_range.pdf')
    plt.show()
    
    #plot emissivity versus density
    text = orig_text + ', ' + temp
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.tick_params(axis='x', labelsize=18)
    ax2.tick_params(axis='y', labelsize=18)
    ax2.set_xlabel('Density in cm$^{-3}$', fontsize=18)
    ax2.set_ylabel('Emissivity in $\mathit{ph}$ $cm^3$ $s^{-1}$ $bin^{-1}$', fontsize=18)
    anchored_text = AnchoredText(text, loc='lower right', frameon=False)
    ax2.add_artist(anchored_text)
    ax2.loglog(dens_bins, dens_eps_orig, label='Original emissivity')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_dens.pdf')
    ax2.fill_between(dens_bins, dens_eps_min, dens_eps_max, alpha=0.5, color='g', label='Range of emissivities')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_dens_range.pdf')
    plt.show()
    
def wrapper_run_line_diagnostics(table1, inputs1, values1, table2, inputs2, values2):
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    
    """ Table1 and table2 are tables from individually run get_tables() on the two transitions
    specified by the user. Inputs1, values1, inputs2, values2, are  dictionaries holding the
    inputs and the sensitivity table/emissivity values for each transition respectively.
    
    Varies temperature and density separately for each transition and recalculates emissivity.
    
    Returns dictionary line_ratio_diagnostics containing arrays of the temperature and density bins
    for each transition, as well as the original, minimum, and maximum line ratios calculated
    from varying temperature and density independently."""
<<<<<<< HEAD
    type = line_diagnostics_1.get('type')

    if type == 'temp':
        type1, temp_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1 = [line_diagnostics_1.get(k) for k in line_diagnostics_1]
        type2, temp_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2 = [line_diagnostics_2.get(k) for k in line_diagnostics_2]
        Te_line_ratios = [x / y for x, y in zip(Te_eps_orig1, Te_eps_orig2)]
        Te_line_ratios_min = [x / y for x, y in zip(Te_eps_min1, Te_eps_max2)]
        Te_line_ratios_max = [x / y for x, y in zip(Te_eps_max1, Te_eps_min2)]

        line_ratio_diagnostics = {'type': 'temp', 'temps': numpy.asarray(temp_bins1), 'orig': numpy.asarray(Te_line_ratios),
                                       'min': numpy.asarray(Te_line_ratios_min), 'max': numpy.asarray(Te_line_ratios_max)}

    elif type == 'dens':
        type1, dens_bins1, dens_eps_orig1, dens_eps_min1, dens_eps_max1 = [line_diagnostics_1.get(k) for k in line_diagnostics_1]
        type2, dens_bins2, dens_eps_orig2, dens_eps_min2, dens_eps_max2 = [line_diagnostics_2.get(k) for k in line_diagnostics_2]
        dens_line_ratios = [x / y for x, y in zip(dens_eps_orig1, dens_eps_orig2)]
        dens_line_ratios_min = [x / y for x, y in zip(dens_eps_min1, dens_eps_max2)]
        dens_line_ratios_max = [x / y for x, y in zip(dens_eps_max1, dens_eps_min2)]

        line_ratio_diagnostics = {'type': 'dens', 'dens': numpy.asarray(dens_bins1), 'orig': numpy.asarray(dens_line_ratios),
                                       'min': numpy.asarray(dens_line_ratios_min), 'max': numpy.asarray(dens_line_ratios_max)}
    elif type == 'both':
        type1, temp_bins1, dens_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1, dens_eps_orig1, \
            dens_eps_min1, dens_eps_max1 = [line_diagnostics_1.get(k) for k in line_diagnostics_1]
        type2, temp_bins2, dens_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2, dens_eps_orig2, \
            dens_eps_min2, dens_eps_max2 = [line_diagnostics_2.get(k) for k in line_diagnostics_2]

        Te_line_ratios = [x/y for x,y in zip(Te_eps_orig1, Te_eps_orig2)]
        Te_line_ratios_min = [x/y for x,y in zip(Te_eps_min1, Te_eps_max2)]
        Te_line_ratios_max = [x/y for x,y in zip(Te_eps_max1, Te_eps_min2)]
        dens_line_ratios = [x/y for x,y in zip(dens_eps_orig1, dens_eps_orig2)]
        dens_line_ratios_min = [x/y for x,y in zip(dens_eps_min1, dens_eps_max2)]
        dens_line_ratios_max = [x/y for x,y in zip(dens_eps_max1, dens_eps_min2)]

        line_ratio_diagnostics = {'type': 'both', 'temps': numpy.asarray(temp_bins1), 'dens': numpy.asarray(dens_bins1),
                'Te_orig': numpy.asarray(Te_line_ratios), 'Te_min': numpy.asarray(Te_line_ratios_min),
                'Te_max': numpy.asarray(Te_line_ratios_max), 'dens_orig': numpy.asarray(dens_line_ratios),
                'dens_min': numpy.asarray(dens_line_ratios_min), 'dens_max': numpy.asarray(dens_line_ratios_max)}

    if plot == True:
        plot_line_ratio_diagnostics(transition1, transition2, vary, delta_r, line_ratio_diagnostics)

    return line_ratio_diagnostics

def plot_line_ratio_diagnostics(transition, transition2, vary, delta_r, line_ratio_diagnostics):
    
    """ Line_ratio_diagnostics is a dictionary from run_line_ratio_diagnostics() containing
=======
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs1.get(k) for k in inputs1]

    line_ratio_diagnostics = {}
    
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    line1=str(transition[0])+'-'+str(transition[1])
    line2=str(transition_2[0])+'-'+str(transition_2[1])
    line_ratio = line1+'/'+line2
    name=element+'_'+ion+'_'+process+'_'+line_ratio
    temp = '%.E' % Decimal(Te) +'K'
    percentage = '{0:g}'.format(delta_r*100)+'%'
    density = 'dens='+str(dens)
    if process == 'A':
        text= line_ratio + ' A value $\Delta$ $\pm$' + percentage 
    elif process == 'exc':
        text= line_ratio + ' Direct excitation rate $\Delta$ $\pm$' + percentage
    
    print("////////////////////////////////////////////////////////////////")
    print("Running diagnostics for line ratio", line_ratio)
    print("////////////////////////////////////////////////////////////////")
    
    line_diagnostics1 = run_line_diagnostics(table1, inputs1, values1, transition)
    line_diagnostics2 = run_line_diagnostics(table2, inputs2, values2, transition_2)
    
    temp_bins1, dens_bins1, Te_eps_orig1, Te_eps_min1, Te_eps_max1, dens_eps_orig1, \
        dens_eps_min1, dens_eps_max1, name1, label1, transition1 = [line_diagnostics1.get(k) for k in line_diagnostics1]
    temp_bins2, dens_bins2, Te_eps_orig2, Te_eps_min2, Te_eps_max2, dens_eps_orig2, \
        dens_eps_min2, dens_eps_max2, name2, label2, transition2 = [line_diagnostics2.get(k) for k in line_diagnostics2]
    
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
                                   'ratio': line_ratio, 'name': name, 'label': label1})
        
        
    file = open('line_ratio_diagnostics', 'wb')
    pickle.dump(line_ratio_diagnostics, file)
    file.close()
    
    return line_ratio_diagnostics
    
def plot_line_ratio_diagnostics(inputs, line_ratio_diagnostics):
    
    """ Line_ratio_diagnostics is a dictionary from wrapper_run_line_diagnostics() containing
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
    arrays of temperature and density bins, and line ratios (original, min, max) calculated
    from varying temperature and density.
    
    Plots the line ratios of emissivities for specified two transitions as a function of
    temperature and density."""

<<<<<<< HEAD
    print("Plotting line ratio diagnostics.")
    type = line_ratio_diagnostics.get('type')

    text = '{0:g}'.format(transition[0]) + '->' + '{0:g}'.format(transition[1]) + \
        ' / ' + '{0:g}'.format(transition2[0]) + '->' + '{0:g}'.format(transition2[1])

    if type == 'temp':  # plot emissivity versus temperature
        type, temps, Te_orig, Te_min, Te_max = [line_ratio_diagnostics.get(k) for k in line_ratio_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax.set_ylabel('Line Ratio')
        ax.semilogx(temps, Te_orig, label='Original', color='b')
        ax.fill_between(temps, Te_min, Te_max, alpha=0.5, color='b', label="Range")
        min, max = [x / y for x, y in zip(Te_min, Te_orig)], [x / y for x, y in zip(Te_max, Te_orig)]
        ax2.axhline(y=1, color='k')
        ax2.fill_between(temps, min, max, label='Range', color='b', alpha=0.5)
        ax2.set_xlabel('Temperature in K')
        ax2.set_ylabel('New/Original Ratio')
    elif type == 'dens':
        type, dens, dens_orig, dens_min, dens_max = [line_ratio_diagnostics.get(k) for k in line_ratio_diagnostics]
        fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax.set_ylabel('Line Ratio')
        ax.set_xlabel('Density in cm$^{-3}$')
        ax.semilogx(dens, dens_orig, label='Original', color='b')
        ax.fill_between(dens, dens_min, dens_max, alpha=0.5, color='b', label="Range")
        min, max = [x / y for x, y in zip(dens_min, dens_orig)], [x / y for x, y in zip(dens_max,dens_orig)]
        ax2.axhline(y=1, color='k')
        ax2.fill_between(dens, min, max, label='Range', color='b', alpha=0.5)
        ax2.set_xlabel('Density in cm$^{-3}$')
        ax2.set_ylabel('New/Original Ratio')
    elif type == 'both':
        type, temps, dens, Te_orig, Te_min, Te_max, dens_orig, dens_min, dens_max = \
            [line_ratio_diagnostics.get(k) for k in line_ratio_diagnostics]
        fig, ax = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col')
        fig.suptitle(text + ' ' + vary + ' ' + '$\pm$' + str(delta_r*100)+'%')
        ax[0, 0].set_ylabel('Line Ratio')
        for ax in [ax[0, 0], ax[1, 0]]:
            ax.set_xlabel('Temperature in K')
        for ax in [ax[0, 1], ax[1, 1]]:
            ax.set_xlabel('Density in cm$^{-3}$')
        ax[1, 0].set_ylabel('New/Original Ratio')
        ax[0, 1].semilogx(dens, dens_orig, label='Original', color='b')
        ax[0, 0].semilogx(temps, Te_orig, label='Original', color='b')
        ax[0, 1].fill_between(dens, dens_min, dens_max, alpha=0.5, color='b',label="Range")
        ax[0, 0].fill_between(temps, Te_min, Te_max, color='b', alpha=0.5, label="Range")
        ax[1, 0].axhline(y=1, color='k')
        ax[1, 1].axhline(y=1, color='k')
        min, max = [x / y for x, y in zip(dens_min, dens_orig)], [x / y for x, y in zip(dens_max, dens_orig)]
        ax[1, 1].fill_between(dens, min, max, label='Range', color='b', alpha=0.5)
        min, max = [x / y for x, y in zip(Te_min, Te_orig)], [x / y for x, y in zip(Te_max, Te_orig)]
        ax[1, 0].fill_between(temps, min, max, label='Range', color='b', alpha=0.5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0, top=0.86)
    plt.legend(fontsize='xx-small')
    plt.savefig('test plot ratio diagnostics.pdf')
    plt.show()
    plt.close('all')

def plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, vary, lines, wavelen={}, corrthresh={}, e_signif={}):

    """ Plots sensitive epsilons for multiple transitions all on one plot.
    Lines is dict of {'name': (up, lo)} or just list of [(up, lo), (up2, lo2)]
    Can specify wavelength range as well as correlation threshold and minimum
    emissivity (e_signif) for data plotted. Corrthresh refers to fractional change
    in emissivity, |dE/E|. Vary is either 'A' value or 'exc' rate."""

    if vary == 'A': processes = ['A']
    elif vary == 'exc': processes = ['exc']
    elif vary == 'both': processes = ['A', 'exc']
    if corrthresh == {}: corrthresh = 0.02
    if isinstance(lines,list):
        set = {}
        for (up, lo) in lines:
            set.update({str(up)+'->'+str(lo): (up,lo)})
    elif isinstance(lines, dict):
        set = lines

    for process in processes:
        plt.figure()
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('% Emissivity Change')
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
        plt.title(text, fontsize=16)

        for k in set:
            transition = set.get(k)
            if process == 'A':
                extras = {'process': process, 'delta_r': delta_r, 'transition': transition, 'transition_2': None,
                          'npnts': 2, 'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': corrthresh, 'e_signif': 0.0}
                inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_a(inputs, values, transition)
            elif process == 'exc':
                extras = {'process': process, 'delta_r': delta_r, 'transition': transition[::-1], 'transition_2': None,
                          'npnts': 2,'wavelen': {}, 'Te_range': {}, 'dens_range': {}, 'corrthresh': corrthresh, 'e_signif': 0.0}
                inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_exc(inputs, values, transition)
            table, new_table, inputs, results = get_tables(new_inputs, new_values)

            # filter data for significance, only want other lines changed
            for i in range(len(new_table)+1):
                if (new_table['Upper'][i],new_table['Lower'][i]) == transition or (new_table['Upper'][i],new_table['Lower'][i]) == transition[::-1]:
                    new_table.remove_rows([i])
                    break
            cutoff_data = new_table[new_table['|dE/E|'] >= corrthresh]
            if wavelen != {}:
                cutoff_data = cutoff_data[wavelen[0] < cutoff_data['Lambda'] < wavelen[1]]
            if e_signif != {}:
                cutoff_data = cutoff_data[cutoff_data['Epsilon_orig'] > e_signif]

            # plot wavelength vs. % emissivity change
            plt.semilogy(cutoff_data['Lambda'], cutoff_data['|dE/E|'], linestyle='', marker='o', label=k)

            # label each point w/ transition
            # transition_labels = []
            # for x in cutoff_data:
            #     transition_name = '{0:g}'.format(x['Upper']) + '->' + '{0:g}'.format(x['Lower'])
            #     transition_labels.append(transition_name)
            # for (x, y, label) in zip(numpy.array(cutoff_data['Lambda']), numpy.array(cutoff_data['|dE/E|']),
            #                          transition_labels):
            #     plt.annotate(label, xy=(x, y))

            plt.tight_layout()
            if vary == 'exc':
                plt.legend(fontsize='xx-small', loc='lower right')
            elif vary == 'A':
                plt.legend(fontsize='xx-small', loc='lower right')

    file_name = 'sensitive lines '+vary+' '+element + str(z1) + '_'
    for number in range(1, 20, 1):
        file = pathlib.Path(file_name + str(number) + '.pdf')
        if file.exists():
            continue
        else:
            plt.savefig(file_name+ str(number) + '.pdf')
            break

    plt.show()

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

            diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type='temp', num=num, plot=False)
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

            diagnostics = variableapec.run_line_diagnostics(table, inputs, values, transition, type='dens', num=num, plot=False)

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
  #factor = delta_r
  Telist = numpy.logspace(4, 9, 1251)
  element = pyatomdb.atomic.Ztoelsymb(Z)
  z1_test = z1

  # Change ionization
  varyir = 'i'

  eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, varyir, delta_r)

  # find peak temperatures
  peak = numpy.max(eqpopn[:, z1_test - 1])
  index = numpy.argmax(eqpopn[:, z1_test-1])
  next_peak = numpy.max(eqpopn[:, z1_test])
  next_index = numpy.argmax(eqpopn[:, z1_test])

  print("Original abundance at", frac, "*peak is", peak, "or in -Log10T =", -numpy.log10(frac*peak))
  print("New range of", frac, "*peak is", numpy.max(negpopn[:, z1_test-1]), "to", numpy.max(pospopn[:, z1_test-1]))
  print("in -Log10T this is", -numpy.log10(numpy.max(negpopn[:, z1_test-1])), "to", -numpy.log10(numpy.max(pospopn[:, z1_test-1])))

  #interpolate
  min_temp = numpy.interp(frac * peak, negpopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_min_temp = numpy.interp(frac * next_peak, negpopn[:next_index, z1_test], Telist[:next_index])
  max_temp = numpy.interp(frac * peak, pospopn[index:, z1_test - 1][::-1], Telist[index:][::-1])
  next_max_temp = numpy.interp(frac * next_peak, pospopn[:next_index, z1_test], Telist[:next_index])

  orig_dt = max_temp - min_temp
  orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
  next_dt = next_max_temp - next_min_temp
  next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)

  print("Ionize - comparing temps: orig", Telist[index], "min:", min_temp, "max:", max_temp)


  dex_changes = {'Uncertainty': '+/-' + str(delta_r * 100) + '%', 'Ion frac': str(frac) + ' peak',
                 'Unit': 'dex', element + ' ' + str(z1_test - 1) + '+ (ionize)': str(orig_dex_change),
                 element + ' ' + str(z1_test) + '+ (ionize)': str(next_dex_change)}
  temp_changes = {'Uncertainty': ' ', 'Ion frac': ' ',
                  'Unit': 'Kelvin', element + ' ' + str(z1_test - 1) + '+ (ionize)': str(orig_dt),
                  element + ' ' + str(z1_test) + '+ (ionize)': str(next_dt)}

  actual_temps = {'Ion': element + ' ' + str(z1_test-1) + '+', 'frac peak': str(frac) + '*peak',
                  'Error': str(delta_r*100)+'%', 'orig': numpy.log10(Telist[index]),
                  'ionize min': numpy.log10(min_temp), 'ionize max': numpy.log10(max_temp)}

  # Change recombination
  varyir = 'r'
  eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, varyir, delta_r)

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

  orig_dt = max_temp - min_temp
  orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
  next_dt = next_max_temp - next_min_temp
  next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)

  print("Recomb - comparing temps: orig", Telist[index], "min:", min_temp, "max:", max_temp)
  #update with recomb data
  dex_changes.update({element + ' ' + str(z1_test - 1) + '+ (recomb)': str(orig_dex_change),
                      element + ' ' + str(z1_test) + '+ (recomb)': str(next_dex_change)})
  temp_changes.update({element + ' ' + str(z1_test - 1) + '+ (recomb)': str(orig_dt),
                      element + ' ' + str(z1_test) + '+ (recomb)': str(next_dt)})

  actual_temps.update({'recomb min': numpy.log10(min_temp), 'recomb max': numpy.log10(max_temp)})
  names = ['Ion', 'frac peak', 'Error', 'orig', 'ionize min', 'ionize max', 'recomb min', 'recomb max']

  # write to csv file
  fieldnames = ['Uncertainty', 'Ion frac', 'Unit', element + ' ' + str(z1_test - 1) + '+ (ionize)',
                element + ' ' + str(z1_test) + '+ (ionize)',
                element + ' ' + str(z1_test - 1) + '+ (recomb)', element + ' ' + str(z1_test) + '+ (recomb)']

  file = pathlib.Path('temp change ' + element + str(z1_test) + '.csv')
  if file.exists():
    with open('temp change ' + element + str(z1_test) + '.csv', mode='a+') as write_obj:
      dict_writer = csv.DictWriter(write_obj, fieldnames=fieldnames)
      dict_writer.writerow(dex_changes)
      dict_writer.writerow(temp_changes)
    with open('actual dex temps' + element + str(z1_test) + '.csv', mode='a+') as write_obj:
        dict_writer = csv.DictWriter(write_obj, fieldnames=names)
        dict_writer.writerow(actual_temps)
  else:
    with open('temp change ' + element + str(z1_test) + '.csv', mode='w') as csv_file:
      dict_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
      dict_writer.writeheader()
      dict_writer.writerow(dex_changes)
      dict_writer.writerow(temp_changes)
    with open('actual dex temps ' + element +' '+ str(z1_test) + '.csv', mode='w') as obj:
        dict_writer = csv.DictWriter(obj, fieldnames=names)
        dict_writer.writeheader()
        dict_writer.writerow(actual_temps)

def wrapper_find_temp_change(list, frac_list, errors):
  for (Z, z1) in list:
    for frac in frac_list:
      for delta_r in errors:
        print("Finding temp shift at", frac, '* peak for', (Z, z1), 'and', str(delta_r*100), '% uncertainty')
        find_temp_change(Z, z1, frac, delta_r)

def vary_csd(Z, z1, varyir, delta_r):
    """ Varies CSD by changing z1 'i' or 'r' rate only.
    Delta_r is fractional change, i.e. 0.10 for 10%.
    Plots new CSD for Z and changes in slope.
    Returns dictionary of orig min, max populations."""

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

    eqpopn, popspopn, negpopn = get_new_popns(Telist, Z, varyir, delta_r)
    ret = {'orig': eqpopn, 'max': pospopn, 'min': negpopn}

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

    return ret

def monte_carlo_csd(Z, max_ionize, max_recomb, runs=100, Telist=numpy.logspace(4,9,1251), makefiles=False, plot=False):
    """ Varies CSD 100 Monte Carlo calculations with a max error on
    ionization and recombination rates. If need varied eigenfiles,
    set makefiles=True."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    clist = get_cmap(Z+2)

    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp

    mc_popn = numpy.zeros([runs, Z+1, len(Telist)])
    random_ion = numpy.zeros([len(Telist), Z])
    random_rec = numpy.zeros([len(Telist), Z])
    for run in range(runs):
        #get random ionization errors
        lower, upper = -2*max_ionize, 2*max_ionize
        mu, sigma = 0, max_ionize
        random1 = stats.truncnorm.rvs((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=Z)
        #get random recombination errors
        lower, upper = -2 * max_recomb, 2 * max_recomb
        mu, sigma = 0, max_recomb
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

        for z1 in range(1, Z + 2):
            mc_popn[run, z1 - 1, :] = newpopn[:, z1 - 1]

    # find mean and standard dev from all runs for each ion at each temp
    median = numpy.zeros([len(Telist), Z + 1])
    min = numpy.zeros([len(Telist), Z + 1])
    max = numpy.zeros([len(Telist), Z + 1])
    for z1 in range(1, Z + 2):
        for i in range(len(Telist)):
            pop = mc_popn[:, z1 - 1, i]
            pop_list = numpy.sort(pop)
            median[i, z1-1] = numpy.median(pop_list)
            min_index = int(0.16 * runs - 1)
            max_index = int(0.84 * runs - 1)
            min[i, z1-1] = pop_list[min_index]
            max[i, z1-1] = pop_list[max_index]

    if plot == True:
        # set up plot
        fig, ax = plt.subplots()
        ax.set_xlabel('Temperature in K', fontsize=12)
        ax.set_ylabel('Ion Fraction', fontsize=12)
        for z1 in range(1, Z+2):
            ax.semilogx(Telist, median[:, z1-1], color=clist(z1-1), linestyle='-')
            ax.fill_between(Telist, min[:, z1-1], max[:, z1-1], color = clist(z1-1), alpha=0.4)
        fig.suptitle('Ionize error between +/-'+str(max_ionize*100)+'% \nRecomb error between' +
                        '+/-' +str(max_recomb*100)+'%')
        plt.savefig(element+' Monte Carlo CSD.pdf')
        plt.show()
        plt.close('all')

    return median, min, max

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

def gaussian_varied_rates(Z, max_error, Telist):
    """Applies random error from Gaussian distribution
    of sigma = max_error for all rate coefficients
    in an element over Telist range.
    Returns arrays of varied ion and rec rates."""

    #get original rates
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

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

def get_partial_deriv(Z, z1, vary, errors, Te={}, dens={}, lines={}):
    """Creates .csv file with partial log derivatives
    dE/dR and dE/E for specified error(s) for key lines
    and other affected lines.

    dens: int or str 'critical'
        If no dens supplied, uses dens=1
        If 'critical', will use critical dens
    errors: list []
    vary: str of either 'exc', 'A', or 'both'
    Te: float
        If no Te supplied, uses Te where ion peaks in plasma
    lines: list []
        If key lines not specified,
        will vary H-like, He-like, or Fe line complexes."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    print("Calculating partial derivatives for", element, str(z1))

    if Te == {}:
        #find temperature where ion fraction peaks
        Te = find_peak_Te(Z, z1, unit='K')
        print("Using peak Te =", Te)

    if lines == {}:
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
    elif lines != {}: trans_list = lines

    if dens=={}: dens=1
    elif dens=='critical':
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

    uppers, lowers = [x[0] for x in trans_list], [x[1] for x in trans_list]

    #get orig emissivities so use any random transition for set up
    inputs, values = set_up(Z, z1, Te, dens)
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
    t2 = Table([uppers, lowers, wavelen, orig], names=('Upper', 'Lower', 'Lambda', 'Orig'))
    #t3 = Table(names=('Upper', 'Lower', 'Lambda', 'Orig'))
    #t4 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))
    #t5 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))
    #t6 = Table(names=('Upper', 'Lower', 'Lambda','Orig'))
    t3, t4, t5, t6 = t2.copy(), t2.copy(), t2.copy(), t2.copy()

    for delta_r in errors:
        for transition in trans_list:
            print("Varying the transition", transition, 'by', str(delta_r*100)+'%')
            if (vary == 'exc') or (vary == 'both'): #vary exc rate
                extras={'process':'exc', 'delta_r': delta_r, 'transition':transition[::-1], 'transition_2':[],
                    'npnts':2, 'wavelen':(10,20), 'Te_range':(Te/10, Te*10), 'dens_range':(1,10e16),
                        'corrthresh':0.0, 'e_signif':0.0}
                npnts=2
                inputs, values, exc_transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_exc(inputs, values, exc_transition)
                table, new_table, inputs, results = get_tables(new_inputs, new_values)
                partial_deriv, frac_E = [], []
                for x in t:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if frac < 0.0001: frac = 0.0
                            frac_E.append(frac)
                            partial_deriv.append(deriv)

                name_1 = 'dE/dR exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                name_2 = 'dE/E exc ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                t[name_1], t2[name_2] = partial_deriv, frac_E

                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t.add_columns([deriv])
                # t2.add_columns([frac])

                #remove rows of lines already checked separately because tables are ordered differently
                idx1, idx2 = [],[]
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
                    if 0.02 <= x['|dE/E|'] <= 0.98:
                        print("Found interesting line with dE/dR = ", x['dE/dR'])
                        if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                            row = [x['Upper'], x['Lower'], x['Lambda'], x['Epsilon_orig']]+['0']*(len(t3.colnames)-4)
                            t3.add_row(row)
                            t4.add_row(row)
                    #linear change lines
                    elif 0.98 <= x['|dE/E|'] <= 1.02:
                        print("Found linear change line with dE/dR = ", x['dE/dR'])
                        if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                            t5.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']]+['0']*(len(t5.colnames)-4))
                            t6.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))

                #get partial derivs for interesting lines
                partial_deriv, frac_E = [], []
                for x in t3:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if abs(frac) < 0.0001: frac = 0.0
                            partial_deriv.append(deriv)
                            frac_E.append(frac)
                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t3.add_columns([deriv])
                # t4.add_columns([frac])
                t3[name_1], t4[name_2] = partial_deriv, frac_E

                #get partial derivs for lines that change linearly
                partial_deriv, frac_E = [], []
                for x in t5:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if abs(frac) < 0.0001: frac = 0.0
                            partial_deriv.append(deriv)
                            frac_E.append(frac)
                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t5.add_columns([deriv])
                # t6.add_columns([frac])
                t5[name_1], t6[name_2] = partial_deriv, frac_E

            if (vary == 'A') or (vary == 'both'):
                #vary A value
                extras.update({'process':'A', 'transition':transition})
                inputs, values, transition = set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = vary_a(inputs, values, transition)
                table, new_table, inputs, results = get_tables(new_inputs, new_values)
                partial_deriv, frac_E = [], []
                for x in t:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if frac < 0.0001: deriv = 0.0
                            partial_deriv.append(deriv)
                            frac_E.append(frac)

                name_1 = 'dE/dR A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                name_2 = 'dE/E A ' + str(delta_r * 100) + '% ' + str(transition[0]) + '->' + str(transition[1])
                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t.add_columns([deriv])
                # t2.add_columns([frac])
                t[name_1], t2[name_2] = partial_deriv, frac_E

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
                    if 0.02 <= x['|dE/E|'] <= 0.98: #"interesting" lines
                        print("Found additional line impacted - ", x['|dE/E|'])
                        if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                            t3.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t3.colnames) - 4))
                            t4.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t3.colnames) - 4))

                    elif 0.98 <= x['|dE/E|'] <= 1.02: #linear change lines
                        print("Found additional line impacted - ", x['|dE/E|'])
                        if (x['Upper'], x['Lower']) not in zip(t3['Upper'], t3['Lower']):
                            t5.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))
                            t6.add_row([x['Upper'], x['Lower'], x['Lambda'],x['Epsilon_orig']] + ['0'] * (len(t5.colnames) - 4))

                # get partial derivs for interesting lines
                partial_deriv, frac_E = [], []
                for x in t3:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if abs(frac) < 0.0001: frac = 0.0
                            partial_deriv.append(deriv)
                            frac_E.append(frac)
                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t3.add_columns([deriv])
                # t4.add_columns([frac])
                t3[name_1], t4[name_2] = partial_deriv, frac_E

                # get partial derivs for lines that change linearly
                partial_deriv, frac_E = [], []
                for x in t5:
                    for y in new_table:
                        if (y['Upper'], y['Lower']) == (x['Upper'], x['Lower']):
                            deriv, frac = y['dE/dR'], y['|dE/E|']
                            if abs(deriv) < 0.0001: deriv = 0.0
                            if abs(frac) < 0.0001: frac = 0.0
                            partial_deriv.append(deriv)
                            frac_E.append(frac)
                # deriv = Column(name=name_1, data=partial_deriv)
                # frac = Column(name=name_2, data=frac_E)
                # t5.add_columns([deriv])
                # t6.add_columns([frac])
                t5[name_1], t6[name_2] = partial_deriv, frac_E

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
    Telist, ion_frac = numpy.logspace(4, 9, 51), []
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for temp_z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=temp_z1, extrap=True,settings=False)
        ionlist[:, temp_z1 - 1] = iontmp
        reclist[:, temp_z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)
    idx = numpy.argmax(eqpopn[:, z1-1])
    peak_Te = Telist[idx]

    if unit == 'keV':
        peak_Te = peak_Te/11604525.0061657

    return peak_Te

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calc_direct_exc_emiss(Z, z1, up, lo, Te, dens, unit='K', datacache={}):
    if unit == 'keV': Te = Te/11604525.0061657
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=True, \
                                                        do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
    lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV')
    nlev = len(lvdat[1].data)

    matrix, B = numpy.zeros((nlev, nlev)), numpy.zeros(nlev)
    # populate full CR matrix by summing rates for all processes
    for i in range(len(init)):
        x, y = final[i], init[i]
        matrix[x][y] += rates[i]

    # find fraction of each ion in plasma to multiply epsilons by
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=datacache)

    # set up and solve CR matrix for level populations
    matrix[0][:], B[0] = 1.0, 1.0
    lev_pop = numpy.linalg.solve(matrix, B)
    lev_pop *= pop_fraction[z1 - 1]

    exc_rate = pyatomdb.atomdb.get_maxwell_rate(Te, Z=Z, z1=1, dtype='EC', finallev=up, initlev=lo, datacache=datacache, exconly=True)
    emiss = exc_rate * lev_pop[lo - 1]

    return emiss[0]

def calc_lev_pop(Z, z1, Te, dens, Teunit='K', datacache={}):
    if Teunit == 'keV': Te = Te / 11604525.0061657
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
    lvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV')
    nlev = len(lvdat[1].data)

    matrix, B = numpy.zeros((nlev, nlev)), numpy.zeros(nlev)
    # populate full CR matrix by summing rates for all processes
    for i in range(len(init)):
        x, y = final[i], init[i]
        matrix[x][y] += rates[i]

    # find fraction of each ion in plasma to multiply epsilons by
    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=datacache)

    # set up and solve CR matrix for level populations
    matrix[0][:], B[0] = 1.0, 1.0
    lev_pop = numpy.linalg.solve(matrix, B)
    lev_pop *= pop_fraction[z1 - 1]

    return lev_pop

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

def get_all_new_emiss(Z, z1, up, lo, Te, dens, vary, delta_r):
    """ Varies Z, z1 'exc' or 'A' by delta_r at specified Te and dens
    and returns table of all lines with original epsilon, dE/dR and dE/E,
    sorted by greatest change in emissivity dE/E to smalelst.
    Writes table to csv file as well."""

    element = pyatomdb.atomic.Ztoelsymb(Z)
    extras = {'process':vary, 'delta_r':delta_r,'transition':(up, lo), 'transition_2': [], 'npnts':2,
            'wavelen': {}, 'Te_range':{},'dens_range': {},'corrthresh':10e-5, 'e_signif':0.0}
    if vary == 'exc':
        extras.update({'transition': (lo, up)})
        inputs, values = variableapec.set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_exc(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)
    elif vary == 'A':
        inputs, values = variableapec.set_up(Z, z1, Te, dens, extras=extras)
        new_inputs, new_values = vary_a(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)

    try:
      new_table.sort('|dE/E|', reverse=True)
    except TypeError:
      new_table.sort('|dE/E|')
      new_table = new_table[::-1]

    for number in range(1, 20):
        fname = 'new epsilons for ' + element + str(z1) + ' ' + str(number) + '.csv'
        file = pathlib.Path(fname)
        if file.exists():
            continue
        else:
            with open(fname, mode='w') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=new_table.colnames)
                writer.writeheader()
                for row in new_table:
                    data = {}
                    for col in new_table.colnames:
                        data.update({col: row[col]})
                    writer.writerow(data)
            break

    print(new_table)
    return new_table

def get_line_emiss(Z, z1, up, lo, vary, ionize_delta_r, recomb_delta_r={}, Te={}, dens=1):
    #don't like this function
    """ Te can be int or list. Returns line emiss values in the following order:
    min, orig, max. If Te input is a list, returns array of length Te.
    Vary can be 'exc' or 'A' to change transition specific rate,
    or 'ir' to change all rate coefficients for Z via Monte Carlo.
    Default is 51 temperatures from 10e4 to 10e9 and dens = 1.
    If only want to specify general max error (i.e. not specific to rate),
    leave recomb_delta_r empty. """

    tau, npnts, Telist = 1e13, 2, numpy.logspace(4, 9, 51)
    if Te == {}: Te = Telist
    if recomb_delta_r == {}: delta_r = ionize_delta_r

    #get line emissivity at particular temp
    if isinstance(Te, int):
        extras = {'process':vary, 'delta_r':delta_r,'transition':(up, lo), 'transition_2': [], 'npnts':2,
            'wavelen':(10,20), 'Te_range':(Te/10, Te*10),'dens_range':(1, 10e16),'corrthresh':10e-5, 'e_signif':0.0}
        inputs, values = variableapec.set_up(Z, z1, Te, dens, extras=extras)
        if vary == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, transition)
            table, new_table, inputs, results = get_tables(new_inputs, new_values)
            for a,b,c,d,e in zip(table['Upper'], table['Lower'], table['Epsilon_orig'], table['Epsilon_1'], table['Epsilon_2']):
                if (a,b) == (up, lo):
                    ret = {'temps': Te, 'orig': c, 'min': d, 'max':e}
        elif vary == 'A':
            new_inputs, new_values = vary_a(inputs, values, transition)
            table, new_table, inputs, results = get_tables(new_inputs, new_values)
            for a,b,c,d,e in zip(table['Upper'], table['Lower'], table['Epsilon_orig'], table['Epsilon_1'], table['Epsilon_2']):
                if (a,b) == (up, lo):
                    ret = {'temps': Te, 'orig': c, 'min': d, 'max': e}
        elif vary == 'ir':
            avg, low, high = monte_carlo_csd(Z, max_ionize, max_recomb)
        return ret

    #get line emissivity as a function of temp:
    elif (isinstance(Te, list)) or (isinstance(Te, numpy.ndarray)):
        extras = {'process': vary, 'delta_r': delta_r,'transition': (up, lo), 'transition_2': [], 'npnts': 2,
                'wavelen': (10, 20), 'Te_range': Te, 'dens_range': (1, 10e16), 'corrthresh': 10e-5, 'e_signif': 0.0}
        inputs, values = variableapec.set_up(Z, z1, Te[0], dens, extras=extras)


        ret = {'temps': Telist, 'orig': orig, 'min': min, 'max': max}

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

def find_max_error_csd(Z, z1, delta_ionize, delta_recomb):     ### Not working
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

def find_new_temp(Z, frac, delta_ionize, delta_recomb, vary='ir', z1={}, unit='K'):
    """ Returns dict of orig/median, min, and max temp values at specified
    fractional abundance (as a fraction of peak). To vary 'i' or 'r'
    rate, need to specify z1. To vary all rate coefficients, set vary
    to 'ir.' Default is vary 'ir' with Monte Carlo CSD and
    return dict of new value arrays indexed by z1-1. Frac can be int
    (single frac) or list. Dict returned is of arrays with frac
    number of rows, z1 columns (i.e. if z1 specified, returns array
    of values for only z1, otherwise no z1 specified, gives new temps
    for all ions).
    """
    if isinstance(frac, float) or isinstance(frac, int): frac = [frac]
    orig_temps = numpy.zeros([len(frac), Z + 1])
    min_temps, max_temps = numpy.zeros([len(frac), Z + 1]), numpy.zeros([len(frac), Z + 1])
    Telist = numpy.logspace(4,9,51)
    if (vary == 'i') or (vary == 'r'):
        eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, varyir, delta_r)
    else: eqpopn, pospopn, negpopn = monte_carlo_csd(Z, delta_ionize, delta_recomb)

    for i in range(len(frac)):
        for z1 in range(1, Z+1):
            peak, index = numpy.max(eqpopn[:, z1 - 1]), numpy.argmax(eqpopn[:, z1 - 1])
            orig_temps[i, z1-1] = numpy.interp(frac[i] * peak, eqpopn[:, z1-1][::-1], Telist[::-1])
            min_temp[i, z1-1] = numpy.interp(frac[i] * peak, negpopn[:, z1 - 1][::-1], Telist[::-1])
            max_temp[i, z1-1] = numpy.interp(frac[i] * peak, pospopn[:, z1 - 1][::-1], Telist[::-1])

    if z1 == {}: ret = {'orig': orig_temps, 'min': min, 'max': max}
    else: ret = {'orig': orig_temps[:, z1-1], 'min': min_temps[:, z1-1], 'max': max_temps[:, z1-1]}
    return ret

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
            pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te[i], teunit='K', datacache=d)
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
        avg, low, high = monte_carlo_csd(Z, delta_ionize, delta_recomb)
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

def temps_of_interest(Z, z1, Telist):
    """Returns dict of temperatures in K and labels of where ion is at
    1%, 10%, peak fractional abundance, and back down to 10% and 1%."""

    z1_test = z1
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])

    # get orig rates
    for z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True)

        ionlist[:, z1 - 1] = iontmp
        reclist[:, z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)

    fracs = [1e-2, 0.1]
    peak = numpy.argmax(eqpopn[:, z1_test - 1])
    increasing_temps = numpy.interp(fracs, eqpopn[:peak, z1_test - 1], Telist[:peak])
    decreasing_temps = numpy.interp(fracs, eqpopn[peak:, z1_test - 1][::-1], Telist[peak:][::-1])
    temp_bins = numpy.array([increasing_temps[0], increasing_temps[1], Telist[peak], decreasing_temps[1], decreasing_temps[0]])

    #ret = {'temps': temp_bins, 'labels': ['1%', '10%', 'peak', '10%', '1%']}
    return temp_bins


def line_sensitivity(Z, z1, up, lo, vary, errors, trans_list, temps={}, dens=1):
    """ Calculates fractional change in emissivity dE/E for
    specified Z (element), z1 (ion charge+1), up, lo line.
    Vary either 'exc' or 'A' for each transition in trans_list
    by each error in errors. Default density is 1 and default
    temps are temps where Z, z1 ionic fraction is 1%, 10%, peak,
    back down to 10% and 1%. Temps can be 'peak' for just the peak
    temperature or a single value.
    trans_list : list of transitions to vary
    vary : str
    errors : list of delta_r's
    """

    if temps == 'peak': temps = [variableapec.find_peak_Te]
    if isinstace(temps, int) or isinstance(temps, float): temps = [temps]
    elif temps == {}: temps = temps_of_interest(Z, z1, numpy.logspace(4,9,1251))

    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)

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
    plt.figure()
    plt.xlabel('Fractional error')
    plt.ylabel('Fractional change in emissivity')
    plt.title(element + ' ' + ion + ' ' + l + ' sensitivity')

    matrix, legend_labels = numpy.zeros([len(trans_list), len(temps), len(errors)]), []
    clist = variableapec.get_cmap(len(temps) + 1)    #color is for # of temps
    markers = ['o', 'v', 's', 'P', '^', '2']    #marker is for each transition
    temp_str = ['%.1E' % Decimal(x) for x in temps]

    for a in range(len(trans_list)):
        transition = trans_list[a]
        legend_labels.append(Line2D([0], [0], marker=markers[a], label=str(transition[0]) + '->' + str(transition[1])))
        for b in range(len(temps)):
            Te, changes = temps[b], []
            legend_labels.append(Line2D([0], [0], color=clist(b), label=temp_str[b]))
            for c in range(len(errors)):
                delta_r, npnts = errors[c], 2
                extras = {'process': vary, 'delta_r': delta_r, 'transition': transition, 'transition_2': [],
                      'npnts': 2, 'wavelen': (10, 20), 'Te_range': {}, 'dens_range': {}, 'corrthresh': 0.0, 'e_signif': 0.0}
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                if vary == 'exc':
                    new_inputs, new_values = variableapec.vary_exc(inputs, values, transition[::-1])
                elif vary == 'A':
                    new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
                table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)
                for y in new_table:
                    if (y['Upper'], y['Lower']) == transition or (y['Upper'], y['Lower']) == transition[::-1]:
                        print("fractional change is", y['|dE/E|'])
                        if y['|dE/E|'] < 0.0001: matrix[a,b,c] = 0.0
                        else: matrix[a,b,c] = y['|dE/E|']
                        changes.append(y['|dE/E|'])
            # for each transition, we'll plot errors vs frac E
            #plt.plot(errors, matrix[a, b, :], color=clist(b), linestyle='--', marker=markers[a])
            plt.plot(errors, changes, color=clist(b), linestyle='--', marker=markers[a])

    plt.legend(handles=legend_labels, fontsize='xx-small', loc='upper left')
    plt.show()

def line_ratio_diagnostics(Z, z1, up, lo, up2, lo2, Te, dens, vary, delta_r, Te_range={}, dens_range={}, num=10,
                           plot=False):
    """ Varies 'exc' or 'A' value of each line individually and recalculates
    line ratio with the various errors."""
    # set defaults
    if Te_range == -1: Te_range = (Te / 10, Te * 10)
    if dens_range == -1: dens_range = (10e0, 10e16)

    # set type
    if (Te_range != {}) & (dens_range != {}):
        type = 'both'
    elif (Te_range == {}) & (dens_range != {}):
        type = 'dens'
    elif (Te_range != {}) & (dens_range == {}):
        type = 'temp'

    line1 = variableapec.run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range=Te_range,
                                              dens_range=dens_range, num=num)
    line2 = variableapec.run_line_diagnostics(Z, z1, up2, lo2, Te, dens, vary, delta_r, Te_range=Te_range,
                                              dens_range=dens_range, num=num)
    line_ratio = variableapec.run_line_ratio_diagnostics((up, lo), (up2, lo2), line1, line2, vary, delta_r, plot=plot)
    return line_ratio

def get_orig_popn(Z, Telist):
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
    for temp_z1 in range(1, Z + 1):
        iontmp, rectmp = pyatomdb.atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=temp_z1, extrap=True, settings=False)
        ionlist[:, temp_z1 - 1] = iontmp
        reclist[:, temp_z1 - 1] = rectmp
    eqpopn = solve_ionrec(Telist, ionlist, reclist, Z)
    return eqpopn

def blended_line_ratio(Z, z1, Te, dens, vary, delta_r, transition_list, denom, type={}, num={}, Te_range={}, dens_range={}, plot=False):
    #specify equation of blended line ratio
    #i.e. denom=1, then blended line ratio = [line 1 + line 2] / line 3
    #denm=2, then blended line ratio = line 1 / [line 2 + line 3]

    #input checking
    if Te_range == {} or Te_range == -1: Te_range = (Te/10, Te*10)
    if dens_range == {} or dens_range == -1: dens_range = (1, 1e16)
    if num == {}: num = 10
    if type == {}: type = 'both'

    emiss1, emiss2, emiss3 = {}, {}, {}
    emiss_list = (1,2,3)
    extras = {'process': vary, 'delta_r': delta_r, 'transition': [], 'transition_2': [], 'npnts': 2,
              'wavelen': (10, 20), 'Te_range': Te_range, 'dens_range': dens_range, 'corrthresh': 10e-5, 'e_signif': 0.0}
    element = pyatomdb.atomic.Ztoelsymb(Z)

    if type == 'temp':
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating temperature diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A':
                extras.update({'transition': transition})
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
            elif vary == 'exc':
                extras.update({'transition': transition[::-1]})
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)
            table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

            lines=[]
            for upper, lower, wavelen in zip(table['Upper'], table['Lower'], table['Lambda']):
                if (upper, lower) == (up, lo):
                    lines.append(wavelen)

            diagnostics = variableapec.run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r,
                            Te_range=Te_range, dens_range=dens_range, num=num, plot=plot)
            type, temp_bins, Te_eps_orig, Te_eps_min, Te_eps_max = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'type': type, 'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}
            if emiss == 2:
                emiss2 = {'type': type, 'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}
            if emiss == 3:
                emiss3 = {'type': type, 'temp_bins': temp_bins, 'Te_min': Te_eps_min, 'Te_max': Te_eps_max, 'Te_orig': Te_eps_orig}

        type, temp_bins1, Te_eps_min, Te_eps_max, Te_eps_orig = [emiss1.get(k) for k in emiss1]
        type, temp_bins2, Te_eps_min2, Te_eps_max2, Te_eps_orig2 = [emiss2.get(k) for k in emiss2]
        type, temp_bins3, Te_eps_min3, Te_eps_max3, Te_eps_orig3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            Te_total_min = Te_eps_min / (Te_eps_min2 + Te_eps_min3)
            Te_total_max = Te_eps_max / (Te_eps_max2 + Te_eps_max3)
            Te_total_orig = Te_eps_orig / (Te_eps_orig2 + Te_eps_orig3)
        elif denom==1:
            Te_total_min = (Te_eps_min + Te_eps_min2) / Te_eps_min3
            Te_total_max = (Te_eps_max + Te_eps_max2) / Te_eps_max3
            Te_total_orig = (Te_eps_orig + Te_eps_orig2) / Te_eps_orig3

        t1 = Table([temp_bins1, Te_total_min, Te_total_max, Te_total_orig],names=('temp', 'Te_min', 'Te_max', 'Te_orig'))
        for number in range(1, 20, 1):
            file = pathlib.Path(element+ '_'+str(z1) + '_blended_ratio_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t1.write(element+ '_'+str(z1) + '_blended_ratio_Te_' + str(number) + '.fits', format='fits')
                break

        if plot == True:
            blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
            fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig.suptitle(blended_ratio)
            ax2.set_xlabel('Temperature in K')
            ax2.set_ylabel('New Ratio/Original')
            ax1.ylabel('Blended Ratio')
            plt.semilogx(temp_bins1, Te_total_orig, label='Original', color='b')
            plt.fill_between(temp_bins1, Te_total_min, Te_total_max, label='Range', color='b', alpha=0.5)
            min = [x / y for x, y in zip(Te_total_min, Te_total_orig)]
            max = [x / y for x, y in zip(Te_total_max, Te_total_orig)]
            ax2.axhline(y=1, color='b')
            ax2.fill_between(temp_bins1, min, max, color='b', alpha=0.5, label='Range')
            plt.tight_layout()
            fig.subplots_adjust(hspace=0, top=0.86)
            ax1.legend(fontsize='xx-small')
            ax2.legend(fontsize='xx-small')
            for number in range(1, 10):
                outfile = pathlib.Path(element + ' ' + str(z1) + ' blended ratio Te' + str(number) + '.pdf')
                if outfile.exists():
                    continue
                else:
                    plt.savefig(element + ' ' + str(z1) + ' blended ratio Te' + str(number) + '.pdf')
                    break
            plt.show()
            plt.close('all')

    elif type == 'dens':
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating density diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A':
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)
            table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

            lines = []
            for upper, lower, wavelen in zip(table['Upper'], table['Lower'], table['Lambda']):
                if (upper, lower) == (up, lo): lines.append(wavelen)

            diagnostics = variableapec.run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r,
                                               Te_range=Te_range, dens_range=dens_range, num=num, plot=plot)

            type, dens_bins, dens_eps_orig, dens_eps_min, dens_eps_max = [diagnostics.get(k) for k in diagnostics]

            if emiss == 1:
                emiss1 = {'dens_bins': dens_bins, 'dens_min': dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
            if emiss == 2:
                emiss2 = {'dens_bins': dens_bins, 'dens_min':dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}
            if emiss == 3:
                emiss3 = {'dens_bins': dens_bins, 'dens_min':dens_eps_min, 'dens_max': dens_eps_max, 'dens_orig': dens_eps_orig}

        type1, dens_bins, dens_eps_min, dens_eps_max, dens_eps_orig = [emiss1.get(k) for k in emiss1]
        type2, dens_bins2, dens_eps_min2, dens_eps_max2, dens_eps_orig2 = [emiss2.get(k) for k in emiss2]
        type3, dens_bins3, dens_eps_min3, dens_eps_max3, dens_eps_orig3 = [emiss3.get(k) for k in emiss3]

        if denom==2:
            dens_total_min = dens_eps_min / (dens_eps_min2 + dens_eps_min3)
            dens_total_max = dens_eps_max / (dens_eps_max2 + dens_eps_max3)
            dens_total_orig = dens_eps_orig / (dens_eps_orig2 + dens_eps_orig3)
        elif denom==1:
            dens_total_min = (dens_eps_min + dens_eps_min2) / dens_eps_min3
            dens_total_max = (dens_eps_max + dens_eps_max2) / dens_eps_max3
            dens_total_orig = (dens_eps_orig + dens_eps_orig2) / dens_eps_orig3

        t2 = Table([dens_bins, dens_total_min, dens_total_max, dens_total_orig], names=('dens', 'dens_min', 'dens_max', 'dens_orig'))
        for number in range(1, 20, 1):
            file = pathlib.Path(element+ '_'+str(z1) + '_blended_ratio_dens_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t2.write(element+ '_'+str(z1) + '_blended_ratio_dens_' + str(number) + '.fits', format='fits')
                break

        if plot == True:
            blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
            fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig.suptitle(blended_ratio)
            ax2.set_xlabel('Density in cm$^{-3}$')
            ax2.set_ylabel('New Ratio/Original')
            ax1.ylabel('Blended Ratio')
            ax1.semilogx(dens_bins, dens_total_orig, label='Original', color='b')
            ax1.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='b', alpha=0.5)
            ax2.axhline(y=1, color='b')
            min = [x / y for x, y in zip(dens_total_min, dens_total_orig)]
            max = [x / y for x, y in zip(dens_total_max, dens_total_orig)]
            ax2.fill_between(dens_bins, min, max, color='b', alpha=0.5, label='Range')
            plt.tight_layout()
            fig.subplots_adjust(hspace=0, top=0.86)
            ax1.legend(fontsize='xx-small')
            ax2.legend(fontsize='xx-small')
            for number in range(1, 10):
                outfile = pathlib.Path(element + ' ' + str(z1) + ' blended ratio dens ' + str(number) + '.pdf')
                if outfile.exists():
                    continue
                else:
                    plt.savefig(element + ' ' + str(z1) + ' blended ratio dens ' + str(number) + '.pdf')
                    break
            plt.show()
            plt.close('all')

    elif type == 'both':
        for transition, emiss in zip(transition_list, emiss_list):
            print('Calculating temperature and density diagnostics for', transition)
            up, lo = transition[0], transition[1]
            if vary == 'A':
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens,extras=extras)
                new_inputs, new_values = variableapec.vary_a(inputs, values, transition)
            elif vary == 'exc':
                extras.update({'transition': (lo, up)})
                inputs, values, transition = variableapec.set_up(Z, z1, Te, dens, extras=extras)
                new_inputs, new_values = variableapec.vary_exc(inputs, values, transition)
            table, new_table, inputs, results = variableapec.get_tables(new_inputs, new_values)

            lines = []
            for upper, lower, wavelen in zip(table['Upper'], table['Lower'], table['Lambda']):
                if (upper, lower) == (up, lo): lines.append(wavelen)

            diagnostics = variableapec.run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r,
                                               Te_range=Te_range, dens_range=dens_range, num=num, plot=plot)

            type, temps, dens, Te_orig, Te_min, Te_max, dens_orig, dens_min, dens_max = [diagnostics.get(k) for k in diagnostics]

            if emiss == emiss1:
                emiss1 = {'type':type, 'temp_bins': temps, 'dens_bins': dens, 'Te_min': Te_min, 'Te_max': Te_max,
                          'Te_orig': Te_orig, 'dens_min': dens_min, 'dens_max': dens_max, 'dens_orig': dens_orig}
            if emiss == emiss2:
                emiss2 = {'type':type, 'temp_bins': temps, 'dens_bins': dens, 'Te_min': Te_min, 'Te_max': Te_max,
                          'Te_orig': Te_orig, 'dens_min': dens_min, 'dens_max': dens_max, 'dens_orig': dens_orig}
            if emiss == emiss3:
                emiss3 = {'type':type, 'temp_bins': temps, 'dens_bins': dens, 'Te_min': Te_min, 'Te_max': Te_max,
                          'Te_orig': Te_orig, 'dens_min': dens_min, 'dens_max': dens_max, 'dens_orig': dens_orig}

        type1, temps1, dens1, Te_min, Te_max, Te_orig, dens_min, dens_max, dens_orig = [emiss1.get(k) for k in emiss1]
        type2, temps2, dens2, Te_min2, Te_max2, Te_orig2, dens_min2, dens_max2, dens_orig2 = [emiss2.get(k) for k in emiss2]
        type3, temps3, dens3, Te_min3, Te_max3, Te_orig3, dens_min3, dens_max3, dens_orig3 = [emiss3.get(k) for k in emiss3]

        if denom == 1:
            Te_total_min = (Te_min + Te_min2) / Te_min3
            Te_total_max = (Te_max + Te_max2) / Te_max3
            Te_total_orig = (Te_orig + Te_orig2) / Te_orig3
            dens_total_min = (dens_min + dens_min2) / dens_min3
            dens_total_max = (dens_max + dens_max2) / dens_max3
            dens_total_orig = (dens_orig + dens_orig2) / dens_orig3
        elif denom == 2:
            Te_total_min = Te_min / (Te_min2 + Te_min3)
            Te_total_max = Te_max / (Te_max2 + Te_max3)
            Te_total_orig = Te_orig / (Te_orig2 + Te_orig3)
            dens_total_min = dens_min / (dens_min2 + dens_min3)
            dens_total_max = dens_max / (dens_max2 + dens_max3)
            dens_total_orig = dens_orig / (dens_orig2 + dens_orig3)

        t1 = Table([temps1, Te_total_min, Te_total_max, Te_total_orig], names=('temp', 'Te_min', 'Te_max', 'Te_orig'))
        t2 = Table([dens1, dens_total_min, dens_total_max, dens_total_orig], names=('dens', 'dens_min', 'dens_max', 'dens_orig'))
        for number in range(1, 20, 1):
            file = pathlib.Path(element+ '_' +str(z1) + '_blended_ratio_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                t1.write(element+ '_' +str(z1) +'_blended_ratio_Te_' + str(number) + '.fits', format='fits')
                t2.write(element+ '_' +str(z1) +'_blended_ratio_dens_' + str(number) + '.fits', format='fits')
                break

        if plot == True:
            blended_ratio = str(lines[0]) + str(lines[1]) + '$\AA$/' + str(lines[2]) + '$\AA$'
            fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig.suptitle(blended_ratio)
            ax1.set_ylabel('Blended ratio') #blended_ratio, fontsize=12)
            ax2.set_xlabel('Temperature in K', fontsize=12)
            ax2.set_ylabel('New Ratio/Original', fontsize=12)
            ax1.semilogx(temp_bins, Te_total_orig, label='Original', color='b')
            ax1.fill_between(temp_bins, Te_total_min, Te_total_max, label='Range', color='b', alpha=0.5)
            min = [x / y for x, y in zip(Te_total_min, Te_total_orig)]
            max = [x / y for x, y in zip(Te_total_max, Te_total_orig)]
            ax2.axhline(y=1, color='b')
            ax2.fill_between(temp_bins, min, max, color='b', alpha=0.5, label='Range')
            plt.tight_layout()
            fig.subplots_adjust(hspace=0, top=0.86)
            ax1.legend(fontsize='xx-small')
            ax2.legend(fontsize='xx-small')

            fig2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
            fig2.suptitle(blended_ratio)
            ax2.set_xlabel('Density in cm$^{-3}$')
            ax2.set_ylabel('New Ratio/Original')
            ax1.ylabel('Blended Ratio')
            ax1.semilogx(dens_bins, dens_total_orig, label='Original', color='b')
            ax1.fill_between(dens_bins, dens_total_min, dens_total_max, label='Range', color='b', alpha=0.5)
            ax2.axhline(y=1, color='b')
            min = [x / y for x, y in zip(dens_total_min, dens_total_orig)]
            max = [x / y for x, y in zip(dens_total_max, dens_total_orig)]
            ax2.fill_between(dens_bins, min, max, color='b', alpha=0.5, label='Range')
            plt.tight_layout()
            fig2.subplots_adjust(hspace=0, top=0.86)
            ax1.legend(fontsize='xx-small')
            ax2.legend(fontsize='xx-small')

            for number in range(1,10):
                outfile = pathlib.Path(element+' '+str(z1)+' blended ratio Te '+str(number)+'.pdf')
                if outfile.exists():
                    continue
                else:
                    fig1.savefig(element + ' ' + str(z1) + ' blended ratio Te ' + str(number) + '.pdf')
                    fig2.savefig(element + ' ' + str(z1) + ' blended ratio dens ' + str(number) + '.pdf')
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

def four_line_diagnostic(Z, z1, Te, dens, vary, delta_r, Te_range={}, dens_range={}, type={}, num={}):
    element = pyatomdb.atomic.Ztoelsymb(Z)

    if type == {}:
        type = 'both'
        if Te_range == -1 or Te_range == {}: Te_range = (Te / 10, Te * 10)
        if (dens_range == -1) or dens_range == {}: dens_range = (1, 1e16)
    elif type == 'temp':
        if Te_range == -1 or Te_range == {}: Te_range = (Te/10, Te*10)
    elif type == 'dens':
        if (dens_range == -1) or dens_range == {}: dens_range = (1, 1e16)
    if num == {}: num = 8

    if (Z % 2) == 0: list = {'r': (7, 1), 'f': (2, 1), 'i': (6, 1), 'i2': (5, 1)}
    else: list = {'r': (7, 1), 'f': (2, 1), 'i': (4, 1), 'i2': (5, 1)}

    for line, transition in list.items():
        up, lo = transition[0], transition[1]
        diagnostics = variableapec.run_line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, Te_range=Te_range,
                                                              dens_range=dens_range, num=num)
        print("For line", line, "diagnostics are:", diagnostics)
        if line == 'r': r_diagnostics = diagnostics
        elif line == 'f': f_diagnostics = diagnostics
        elif line == 'i': i_diagnostics = diagnostics
        elif line == 'i2': i2_diagnostics = diagnostics

    # write diagnostics
    if type == 'temp':
        type, temp_bins1, Te_r_orig, Te_r_min, Te_r_max = [r_diagnostics.get(k) for k in r_diagnostics]
        type2, temp_bins2, Te_f_orig, Te_f_min, Te_f_max = [f_diagnostics.get(k) for k in f_diagnostics]
        type3, temp_bins3, Te_i_orig, Te_i_min, Te_i_max = [i_diagnostics.get(k) for k in i_diagnostics]
        type4, temp_bins4, Te_i_orig2, Te_i_min2, Te_i_max2 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temp_bins1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table.write(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits', format='fits')
                break
    elif type == 'dens':
        type, dens_bins1, dens_r_orig, dens_r_min, dens_r_max = [r_diagnostics.get(k) for k in r_diagnostics]
        type2, dens_bins2, dens_f_orig, dens_f_min, dens_f_max = [f_diagnostics.get(k) for k in f_diagnostics]
        type3, dens_bins3, dens_i_orig, dens_i_min, dens_i_max = [i_diagnostics.get(k) for k in i_diagnostics]
        type4, dens_bins4, dens_i_orig2, dens_i_min2, dens_i_max2 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([dens_bins1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2], names=
                       ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table.write(element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits', format='fits')
                break
    elif type == 'both':
        t1, temps1, dens1, Te_r_orig, Te_r_min, Te_r_max, dens_r_orig, dens_r_min, dens_r_max = [r_diagnostics.get(k) for k in r_diagnostics]
        t2, temps2, dens2, Te_f_orig, Te_f_min, Te_f_max, dens_f_orig, dens_f_min, dens_f_max = [f_diagnostics.get(k) for k in f_diagnostics]
        t3, temps3, dens3, Te_i_orig, Te_i_min, Te_i_max, dens_i_orig, dens_i_min, dens_i_max = [i_diagnostics.get(k) for k in i_diagnostics]
        t4, temps4, dens4, Te_i_orig2, Te_i_min2, Te_i_max2, dens_i_orig2, dens_i_min2, dens_i_max2, name4 = [i2_diagnostics.get(k) for k in i2_diagnostics]

        table = Table([temps1, Te_r_orig, Te_r_min, Te_r_max, Te_f_orig, Te_f_min, Te_f_max,
                       Te_i_orig, Te_i_min, Te_i_max, Te_i_orig2, Te_i_min2, Te_i_max2], names=('temp',
                       'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
                        'i2 orig', 'i2 min', 'i2 max'))
        table2 = Table([dens1, dens_r_orig, dens_r_min, dens_r_max, dens_f_orig, dens_f_min, dens_f_max,
                        dens_i_orig, dens_i_min, dens_i_max, dens_i_orig2, dens_i_min2, dens_i_max2], names=
            ('dens', 'r orig', 'r min', 'r max', 'f orig', 'f min', 'f max', 'i orig', 'i min', 'i max',
             'i2 orig', 'i2 min', 'i2 max'))

        for number in range(1, 20, 1):
            file = pathlib.Path(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits')
            if file.exists():
                continue
            else:
                table2.write(element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits', format='fits')
                table.write(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits', format='fits')
                break

def g_ratio(Z, z1, Te, dens, vary, delta_r, Te_range={}, num={}, need_data=True, plot=True):

    if (Te_range == {}) or (Te_range == -1): Te_range = (Te/10, Te*10)
    if num == {}: num = 8
    if need_data == True:
        four_line_diagnostics(Z, z1, Te, dens, vary, delta_r, Te_range=Te_range, type='temp', num=num)

    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)

    # retrieve diagnostics
    for number in range(20, 0, -1):
        file = pathlib.Path(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits')
        if file.exists():
            hdul = fits.open(element + '_' + str(z1) + '_four_line_Te_' + str(number) + '.fits')
            data = hdul[1].data
            temp_bins, r_orig, r_min, r_max = data['temp'], data['r orig'], data['r min'], data['r max']
            f_orig, f_min, f_max = data['f orig'], data['f min'], data['f max']
            i_orig, i_min, i_max = data['i orig'], data['i min'], data['i max']
            i2_orig, i2_min, i2_max = data['i2 orig'], data['i2 min'], data['i2 max']

            hdul.close()
            break
        else:
            continue

    # do math
    g_min, g_orig, g_max = [], [], []
    for rmin, r, rmax, fmin, f, fmax, imin, i, imax, i2min, i2, i2max in zip(r_min, r_orig, r_max, f_min, f_orig, f_max, \
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

    if plot == True:
        # temp dependent g ratio
        fig, (ax_1, ax_2) = plt.subplots(nrows=2, sharex=True)
        name = element + ' ' + ion + ' G ratio'
        ax_1.set_xlabel('T (K)')
        ax_1.set_ylabel(name)
        ax_1.semilogx(temp_bins, g_orig, label='Original', color='b')
        ax_1.fill_between(temp_bins, g_min, g_max, alpha=0.5, color='b',label="Range")
        min = [a / b for a, b in zip(g_min, g_orig)]
        max = [a / b for a, b in zip(g_max, g_orig)]
        ax_2.axhline(y=1, color='b')
        ax_2.fill_between(temp_bins, min, max, color='b', alpha=0.5, label='Range')
        ax_2.set_ylabel('New Ratio/Original')
        ax_1.legend(fontsize='xx-small')
        ax_2.legend(fontsize='xx-small')
        fig.tight_layout()
        fig.subplots_adjust(hspace=0)
        fig.savefig(element + '_' + ion + '_' + 'G ratio.pdf')

        plt.show()

    return {'temps': temp_bins, 'orig': g_orig, 'min': g_min, 'max':g_max}

def r_ratio(Z, z1, Te, dens, vary, delta_r, dens_range={}, num={}, need_data=True, plot=True):

    if (dens_range == {}) or (dens_range == -1): dens_range = (1, 1e16)
    if num == {}: num = 8
    if need_data == True:
        four_line_diagnostics(Z, z1, Te, dens, vary, delta_r, dens_range=dens_range, type='dens', num=num)

    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)

    # retrieve diagnostics
    for number in range(20, 0, -1):
        file = pathlib.Path(element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits')
        if file.exists():
            hdul = fits.open(element + '_' + str(z1) + '_four_line_dens_' + str(number) + '.fits')
            data = hdul[1].data
            dens_bins, r_orig, r_min, r_max = data['dens'], data['r orig'], data['r min'], data['r max']
            f_orig, f_min, f_max = data['f orig'], data['f min'], data['f max']
            i_orig, i_min, i_max = data['i orig'], data['i min'], data['i max']
            i2_orig, i2_min, i2_max = data['i2 orig'], data['i2 min'], data['i2 max']
            hdul.close()
            break
        else:
            continue

    # do math
    R_orig, R_min, R_max = [], [], []
    for rmin, r, rmax, fmin, f, fmax, imin, i, imax, i2min, i2, i2max in zip(r_min, r_orig, r_max, f_min,
                   f_orig, f_max, i_min, i_orig, i_max,i2_min, i2_orig, i2_max):
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

    if plot == True:
        # density dependent R ratio
        fig2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        fig2.subplots_adjust(hspace=0)
        name = element + ' ' + ion + ' R ratio'
        ax1.set_ylabel(name)
        ax1.set_xlabel('Density in cm$^{-3}$')
        ax1.semilogx(dens_bins, R_orig, label='Original', color='b')
        ax1.fill_between(dens_bins, R_min, R_max, alpha=0.5, color='b', label="Range")
        ax1.legend(fontsize='xx-small')
        min = [a / b for a, b in zip(R_min, R_orig)]
        max = [a / b for a, b in zip(R_max, R_orig)]
        ax2.axhline(y=1, color='b')
        ax2.fill_between(dens_bins, min, max, color='b', alpha=0.5, label='Range')
        ax2.set_ylabel('New Ratio/Original')
        ax2.set_xlabel('Density in cm$^{-3}$')
        ax2.legend(fontsize='xx-small')
        fig2.savefig(element + '_' + ion + '_' + 'R ratio.pdf')
        plt.show()

    return {'dens': dens_bins, 'orig': R_orig, 'min': R_min, 'max': R_max}

def find_CSD_change(Z, z1, delta_ionize, delta_recomb={}, Te={}, frac={}, varyir={}, datacache={}, printout=True):
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

    Telist = numpy.logspace(4, 9, 1251)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    d = {}

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

            eqpopn, pospopn, negpopn = get_new_popns(Telist, Z, z1, varyir, delta_r)

            #if frac abundance specified, find temp change in new CSD at that fraction
            if frac != {}:
                peak = numpy.max(eqpopn[:, z1 - 1])
                index = numpy.argmax(eqpopn[:, z1 - 1])
                next_peak = numpy.max(eqpopn[:, z1])
                next_index = numpy.argmax(eqpopn[:, z1])

                # interpolate
                min_temp = numpy.interp(frac * peak, negpopn[index:, z1-1][::-1], Telist[index:][::-1])
                next_min_temp = numpy.interp(frac * next_peak, negpopn[:next_index, z1], Telist[:next_index])
                max_temp = numpy.interp(frac * peak, pospopn[index:, z1-1][::-1], Telist[index:][::-1])
                next_max_temp = numpy.interp(frac * next_peak, pospopn[:next_index, z1], Telist[:next_index])

                orig_dt = max_temp - min_temp
                orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
                next_dt = next_max_temp - next_min_temp
                next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)

                if printout == True:
                    print(element, str(z1-1), "+ original temperature at", frac, "*peak is", Telist[index], "(Log(T) = ", numpy.log10(Telist[index]))
                    print("\n New temperature range in K is", min_temp, "to", max_temp, "i.e. dT =", orig_dt)
                    print("\n Log(T) range is", numpy.log10(min_temp), "to", numpy.log10(max_temp), "i.e. dLog(T) =", orig_dex_change)
                    # print("\n\nThis also affects", element, str(z1), "+:")
                    # print("Original temperature in K at", frac, "*peak was", Telist[next_index], "(Log(T) = ", numpy.log10(Telist[next_index]))
                    # print("\n New temperature range in K is", next_min_temp, "to", next_max_temp, "i.e. dT =", next_dt)
                    # print("\n Log(T) new range is", numpy.log10(next_min_temp), "to", numpy.log10(next_max_temp), "i.e. dLog(T) =", next_dex_change)
                    print("\n")
                percent_change_Te.append(abs(orig_dt/Telist[index]))

            #if temperature given, find abundance change in new CSD at that temperature
            if Te != {}:
                #original abundance
                pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
                orig = pop_fraction[z1-1]
                next_orig = pop_fraction[z1]

                min = numpy.interp(Te, Telist, negpopn[:, z1 - 1])
                max = numpy.interp(Te, Telist, pospopn[:, z1 - 1])
                next_min = numpy.interp(Te, Telist, negpopn[:, z1])
                next_max = numpy.interp(Te, Telist, pospopn[:, z1])

                #find change
                dA = max - min
                log_dA = -numpy.log10(max) - (-numpy.log10(min))
                next_dA = next_max - next_min
                next_log_dA = -numpy.log10(next_max) - (-numpy.log10(next_min))

                if printout == True:
                    print(element, str(z1 - 1), "+ original abundance at", "%.2E" % Decimal(Te), "K =", orig, "(-Log10(X) =", -numpy.log10(orig), ')\n')
                    print("New abundance range is", min, "to", max, "i.e. d_abund =", dA)
                    print("\n(-Log10(X) new range is", -numpy.log10(min), "to", -numpy.log10(max), '), i.e. d_log_abund =', log_dA)
                    # print("\n\n This also affects neighbor ion", element, str(z1), "+:")
                    # print("\nOriginal abundance at", "%.2E" % Decimal(Te), "K =", next_orig, "(-Log10(X) =", -numpy.log10(next_orig), ')\n')
                    # print("New abundance range is", next_min, "to", next_max, "i.e. d_abund =", next_dA)
                    # print("\n(-Log10(X) new range is", -numpy.log10(next_min), "to", -numpy.log10(next_max), '), i.e. d_log_abund =', next_log_dA)
                    print("\n")
                percent_change_abund.append(abs(dA/orig))

    #if varyir not specified, monte carlo rates
    else:
        print("\n Varying all rate coefficients for", element, '\n')
        eqpopn, pospopn, negpopn = monte_carlo_csd(Z, delta_ionize, delta_recomb)

        for z1 in ions:
            # if frac abundance specified, find temp change in new CSD at that fraction
            if frac != {}:
                peak = numpy.max(eqpopn[:, z1 - 1])
                index = numpy.argmax(eqpopn[:, z1 - 1])
                next_peak = numpy.max(eqpopn[:, z1])
                next_index = numpy.argmax(eqpopn[:, z1])

                # interpolate
                min_temp = numpy.interp(frac * peak, negpopn[index:, z1 - 1][::-1], Telist[index:][::-1])
                next_min_temp = numpy.interp(frac * next_peak, negpopn[:next_index, z1], Telist[:next_index])
                max_temp = numpy.interp(frac * peak, pospopn[index:, z1 - 1][::-1], Telist[index:][::-1])
                next_max_temp = numpy.interp(frac * next_peak, pospopn[:next_index, z1], Telist[:next_index])

                orig_dt = max_temp - min_temp
                orig_dex_change = numpy.log10(max_temp) - numpy.log10(min_temp)
                next_dt = next_max_temp - next_min_temp
                next_dex_change = numpy.log10(next_max_temp) - numpy.log10(next_min_temp)

                if printout == True:
                    print(element, str(z1 - 1), "+ original temperature at", frac, "*peak is", Telist[index], "(Log(T) = ",
                    numpy.log10(Telist[index]))
                    print("\n New temperature range in K is", min_temp, "to", max_temp, "i.e. dT =", orig_dt)
                    print("\n Log(T) range is", numpy.log10(min_temp), "to", numpy.log10(max_temp), "i.e. dLog(T) =",
                          orig_dex_change)
                    # print("\n\nThis also affects", element, str(z1), "+:")
                    # print("Original temperature in K at", frac, "*peak was", Telist[next_index], "(Log(T) = ",
                    #       numpy.log10(Telist[next_index]))
                    # print("\n New temperature range in K is", next_min_temp, "to", next_max_temp, "i.e. dT =", next_dt)
                    # print("\n Log(T) new range is", numpy.log10(next_min_temp), "to", numpy.log10(next_max_temp),
                    #       "i.e. dLog(T) =", next_dex_change)
                    print("\n")
                percent_change_Te.append(abs(orig_dt / Telist[index]))

            # if temperature given, find abundance change in new CSD at that temperature
            if Te != {}:
                # original abundance
                # pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)
                # orig = pop_fraction[z1 - 1]
                # next_orig = pop_fraction[z1]
                orig = numpy.interp(Te, Telist, eqpopn[:, z1-1])
                min = numpy.interp(Te, Telist, negpopn[:, z1 - 1])
                max = numpy.interp(Te, Telist, pospopn[:, z1 - 1])
                print("orig:", orig, "min:", min, "max:", max)
                next_min = numpy.interp(Te, Telist, negpopn[:, z1])
                next_max = numpy.interp(Te, Telist, pospopn[:, z1])

                # find change
                dA = max - min
                log_dA = -numpy.log10(max) - (-numpy.log10(min))
                next_dA = next_max - next_min
                next_log_dA = -numpy.log10(next_max) - (-numpy.log10(next_min))

                if printout == True:
                    print('\n', element, str(z1 - 1), "+ original abundance at", "%.2E" % Decimal(Te), "K =", orig,
                          "(-Log10(X) =", -numpy.log10(orig), ')\n')
                    print("New abundance range is", min, "to", max, "i.e. d_abund =", dA)
                    print("\n(-Log10(X) new range is", -numpy.log10(min), "to", -numpy.log10(max), '), i.e. d_log_abund =',
                          log_dA)
                    # print("\n\n This also affects neighbor ion", element, str(z1), "+:")
                    # print("\nOriginal abundance at", "%.2E" % Decimal(Te), "K =", next_orig, "(-Log10(X) =",
                    #       -numpy.log10(next_orig) + ')\n')
                    # print("New abundance range is", next_min, "to", next_max, "i.e. d_abund =", next_dA)
                    # print("\n(-Log10(X) new range is", -numpy.log10(next_min), "to",
                    #       -numpy.log10(next_max) + '), i.e. d_log_abund =', next_log_dA)
                    print("\n")
                percent_change_abund.append(abs(dA / orig))

    if (Te != {}) & (frac == {}):
        print("Fractional change in abundances:", percent_change_abund)
        return percent_change_abund
    elif (Te == {}) & (frac != {}):
        print("Fractional change in temperatures:", percent_change_Te)
        return percent_change_Te
    else:
        print("Fractional change in abundances:", percent_change_abund, "\nFractional change in temperatures:", percent_change_Te)
        return percent_change_Te, percent_change_abund

def get_new_popns(Telist, Z, z1_test, varyir, delta_r):
    """ Returns equilibrium popn, min, and max from varying single z1
    rate coefficient."""
    factor = delta_r
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

    return eqpopn, pospopn, negpopn

def error_analysis(Z, z1, up, lo, Te, dens, errors, linewidth=2.5, filename={}):
    """ Generates error analysis PDF of CSD and emissivity errors
    for Z, z1 and transition (up, lo) at Te in K. Linewidth is in eV.
    Errors is a list [] of fractional error, i.e. 0.1 for 10%, errors
    should contain only one error if error is not sensitive to rate type
    (ionization or recombination), otherwise errors = [ionize, recomb].
    If filename for output pdf has spaces, they will be replaced
    with underscores (_), i.e. 'O 7 analysis' makes O_7_analysis.pdf
    Default file generated is ErrorAnalysis.pdf"""

    print("Linewidth is", linewidth, "eV")
    if len(errors) == 1: delta_ionize = delta_recomb = errors[0]
    elif len(errors) == 2: delta_ionize, delta_recomb = errors[0], errors[1]
    if filename == {}: filename = 'ErrorAnalysis'
    if '' in filename:
        filename = filename.split()
        filename = '_'.join(filename)
    element = pyatomdb.atomic.Ztoelsymb(Z)
    ion = pyatomdb.atomic.int_to_roman(z1)
    name = element + ' ' + ion
    d, dir = {}, os.environ['ATOMDB'] + '/APED/' + element.lower() + '/' + element.lower() + '_' + str(z1) + '/'

    # get LV data
    lvdata = pyatomdb.atomdb.get_data(Z, z1, 'LV', datacache=d)
    lvdata = lvdata[1].data
    up_config, lo_config = lvdata['ELEC_CONFIG'][up - 1], lvdata['ELEC_CONFIG'][lo - 1]

    # get emissivity and multiply by elemental abundance
    kT = Te / 11604525.0061657  # temp in keV
    ab = pyatomdb.atomdb.get_abundance()[Z]  # abundance of element Z
    s = pyatomdb.spectrum.NEISession(elements=[Z])
    tau = 1e13  # leave fixed
    ret = s.return_line_emissivity(kT, tau, Z, z1, up, lo, init_pop=kT)
    emiss = ret['epsilon']  # already abundance corrected

    # get A values
    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', datacache=d)
    in_range = ladat[1].data
    cascades, cascade_emiss = {}, 0
    lev_pop = calc_lev_pop(Z, z1, Te, dens, Teunit='K', datacache=d)
    other_cascades = 0
    for x in in_range:
        if (x['UPPER_LEV'], x['LOWER_LEV']) == (up, lo):
            lambda_obs, lambda_theory = x['WAVE_OBS'], x['WAVELEN']
            print("Lambda obs", lambda_obs, "lambda theory", lambda_theory)
            ref_obs, ref_theory = x['WV_OBS_REF'], x['WAVE_REF']
            A_val, A_ref = x['EINSTEIN_A'], x['EIN_A_REF']
        # get any cascade data
        if x['LOWER_LEV'] == up:
            epsilon = (x['EINSTEIN_A'] * lev_pop[x['UPPER_LEV'] - 1]) * ab
            percent = round((epsilon / emiss) * 100, 2)
            cascade_emiss += epsilon
            if percent > 1:
                cascades.update(
                    {str(x['UPPER_LEV']) + '$\\rightarrow$' + str(x['LOWER_LEV']): str(percent) + '\\%'})
            else:
                other_cascades += epsilon

    # check if references are bibcodes or just text
    citations, refs = {'obs': '', 'theory': '', 'A': ''}, [ref_obs, ref_theory, A_ref]
    for x, y in zip(refs, list(citations.keys())):
        if " " not in x:
            citations[y] = '\\href{https://ui.adsabs.harvard.edu/abs/' + x + '/abstract}{' + x + '})'
        else:
            citations[y] = x

    pop_fraction = pyatomdb.apec.solve_ionbal_eigen(Z, Te, teunit='K', datacache=d)

    # compare emiss to peak line emissivity
    kTlist = numpy.linspace(kT / 10, kT * 10, 51)  # list of temperatures to look at
    ldata = s.return_line_emissivity(kTlist, tau, Z, z1, up, lo)
    peak_emiss = 0
    for i in range(len(kTlist)):
        while ldata['epsilon'][i] > peak_emiss:
            peak_emiss = ldata['epsilon'][i]

    if emiss < peak_emiss: compare_emiss = 'below'
    elif emiss > peak_emiss: compare_emiss = 'above'
    elif emiss == peak_emiss: compare_emiss = 'at'

    # get emissivities from excitation, ionization, recombination
    init, final, rates = pyatomdb.apec.gather_rates(Z, z1, Te, dens, do_la=False, \
                                                    do_ec=True, do_ir=False, do_pc=False, do_ai=False, datacache=d)
    exc_index = ((final == up - 1) & (init != up - 1))
    init, final, rates = init[exc_index], final[exc_index], rates[exc_index]
    rates *= lev_pop[init] * ab * pop_fraction[z1 - 1]
    igood = numpy.where(rates > 0.01 * sum(rates))[0]
    print("List of significant sources of excitation. Probably only ground state in most cases. \n")
    for iigood in igood:
        print(
            'EXC %i -> %i : %e %f%%' % (
            init[iigood], final[iigood], rates[iigood], rates[iigood] * 100. / sum(rates)))
    exc_emiss = sum(rates)

    i_linelist = ionize(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
    r_linelist = recombine(Z, z1, Te, dens, in_range, pop_fraction, datacache=d)
    for x, y in zip(i_linelist, r_linelist):
        if (x['lambda'] == lambda_theory) or (x['lambda'] == lambda_obs):
            ionize_emiss = x['epsilon'] * ab
        if (y['lambda'] == lambda_theory) or (y['lambda'] == lambda_obs):
            recomb_emiss = y['epsilon'] * ab

    # CSD uncertainty
    ions=[]
    print(pop_fraction)
    for i in range(1, Z+2):
        if pop_fraction[i-1] >= 0.001:     #get ions with abund > 0.1%
            print("For z1=", i, "pop frac", pop_fraction[i-1])
            ions.append(i)

    str_ions = [element + '$^{+' + str(x - 1) + '}$' for x in ions]
    #list_ions = '[' + str_ions[0] + ', ' + str_ions[1] + ', ' + str_ions[2] + ']'
    CSD = [round(pop_fraction[x-1]*100,2) for x in ions]
    CSD_errors = find_CSD_change(Z, ions, delta_ionize, delta_recomb, Te=Te, datacache=d, printout=False)  # MC
    CSD_errors = [int(round(x*100)) for x in CSD_errors]
    #flux_per_ion_stage = [round(ionize_emiss/emiss, 2), round((exc_emiss+cascade_emiss)/emiss,2), round(recomb_emiss/emiss,2)]

    #print CSD
    print("\nIons over 1% abundance are", ions, "compared to", [z1-1, z1, z1+1])
    print("\nCSD are", CSD)
    print("\nCSD errors are", CSD_errors)

    # get energies and flux uncertainty
    ionization_energy = pyatomdb.atomdb.get_ionpot(Z, z1, datacache=d)
    ionization_energy /= 1000  # in keV

    if numpy.isnan(lambda_obs):
        excitation_energy = (12.398 / lambda_theory) * 1000
        print("\nCalculating excitation energy with theoretical wavelength.")
    else:
        excitation_energy = (12.398 / lambda_obs) * 1000  # in eV
        print("\nCalculating excitation energy with observational wavelength.")

    print("kT", kT, "and excitation_energy", excitation_energy)
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
    a['epsilon'] *= ab  # should also multiply by ion frac here, don't know
    a = a[a['epsilon'] > 1e-20]
    sat_blend_flux = sum(a['epsilon'])
    new = sum(a['epsilon']) * pop_fraction[z1 - 1]
    print("Comparing sum", sat_blend_flux, "and sum*pop fraction", new)
    DR_list.update({str(len(a)): element + ' ' + pyatomdb.atomic.int_to_roman(z1)})

    other_lines, total_blend_flux = {}, sat_blend_flux
    sess = pyatomdb.spectrum.CIESession()
    #load response file if given:
    # if (arf != {}) & (rmf != {}):
    #     sess.set_response(rmf, arf = arf)
    #     lines = sess.return_linelist(kT, [wvmin, wvmax], apply_aeff=True)
    # elif (rmf != {}):
    #     sess.set_response(rmf)
    #     lines = sess.return_linelist(kT, [wvmin, wvmax], apply_aeff=True)
    # else:
    print("calling lines on kT", kT, "wvmin", wvmin, "wvmax", wvmax)
    lines = sess.return_linelist(kT,[wvmin,wvmax])

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

    # error propagation (sum CSD errors in quadrature)
    #final_error = math.sqrt(CSD_errors[z1-2] ** 2 + CSD_errors[z1-1] ** 2 + CSD_errors[z1] ** 2)
    #print("\nFinal error is", final_error)
    final_error = 1000000

    # write to tex file
    title = 'Uncertainty Analysis for ' + name + ' ' + str(up) + '$\\rightarrow$' + str(lo) + '\n' \
            + 'at ' + '%.2E' % Decimal(Te) + ' K and ' + '%.1E' % Decimal(dens) + ' cm$^{-3}$ \\'
    with open(filename+'.tex', 'w') as file:
        file.write('\\documentclass[11pt]{article}\n\n')
        file.write('\\usepackage{hyperref}\n\n')
        file.write('\\begin{document}\n\n')
        file.write('\\begin{center}\n{\Large ' + title + '\n\\large ')
        file.write('AtomDB v3.0.9; pyAtomDB error analysis v0.1}\n')
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
            file.write(ion + ' & ' + str(frac) + '\\% $\\pm$ ' + str(error/2) + '\\% \\\ \n')
        file.write('\\end{tabular}\n')
        file.write('\\vskip 0.1in \n')

        file.write('\\noindent Detailed breakdown: \\\ \n')
        file.write('\\noindent \\begin{tabular}{lll} \n')
        file.write('\\hline Direct Excitation & & ' + str(round((exc_emiss / emiss) * 100, 1)) + '\\% \\\ \n')

        for key, val in cascades.items():
            file.write('Cascade & ' + key + ' & ' + val + '\\\ \n')

        file.write('Other cascades & & ' + str(round((other_cascades/emiss)*100, 1)) + '\\% \\\ \n')
        file.write('Ionization & & ' + str(round((ionize_emiss / emiss) * 100, 1)) + '\\% \\\ \n')
        file.write('Recombination & & ' + str(round((recomb_emiss / emiss) * 100, 1)) + '\\% \\\ \\hline \n')
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

        if len(DR_list) > 0:
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

        if len(DR_list) > 0:
            file.write('\\indent Note: DR lines have highly uncertain fluxes in general; 300\\% in some cases. ' +
                       '$\lambda$\ also uncertain.  Need to compare to Badnell for totals. \\\ \n')

        file.write('\\end{document}')

    os.system("pdflatex " + filename+'.tex')
    #remove extra files
    for ext in ['.log', '.out', '.aux', '.tex']:
        os.remove(filename+ext)

def find_max_error_csd(Z, z1, delta_ionize, delta_recomb={}):
    """ Finds temperature where CSD range from +/- delta_r
    is the largest. Returns dictionary of 'avg abund', 'temp'
    'CSD error' and 'frac error', where avg abund is the ionic
     fraction at the temperature, CSD error is the max change
     in CSD from the Monte Carlo uncertainty perturbation,
     and frac error is the CSD error as a fractino of avg abund.'"""

    if delta_recomb == {}: delta_recomb = delta_ionize

    median, min, max = variableapec.monte_carlo_csd(Z, delta_ionize, delta_recomb)
    Telist = numpy.logspace(4,9,1251)
    spread, i = 0, 0
    for x,y in zip(min[:, z1-1], max[:, z1-1]):
        if y-x >= spread:
            spread = y-x
            i+=1
        else:
            break

    print("max change in CSD is", spread, "and frac error is", spread/median[i, z1-1], "at temperature", Telist[i])
    ret = {'avg abund': median[i, z1-1], 'temp': Telist[i], 'CSD error': spread, 'frac error': spread/median[i, z1-1]}
    print(ret)
    return ret

def peak_frac_sensitivity(Z, errors, z1_list={}, makefiles=True, plot=True):
    """ Finds fractional change in ion peak abundances from
    varying CSD with Gaussian selected random errors using each delta_r
    from errors as sigma. If z1_list empty, calculates for all ions."""

    Telist, element,  = numpy.logspace(4, 9, 1251), pyatomdb.atomic.Ztoelsymb(Z)
    ions=[element+' +']+[element+' '+str(z1-1)+'+' for z1 in range(2, Z+2)]
    if z1_list == {}:
        z1_list = [z1 for z1 in range(1,Z+2)]
    clist = get_cmap(Z+2)
    t = values = Table()
    t['Ion'] = values['Ion'] = ions
    temp_fracs = numpy.zeros([len(errors), len(Telist), Z+1])

    counter = 0
    for max_error in errors:
        median, min, max = variableapec.monte_carlo_csd(Z, max_error, 100, makefiles=False, plot=False)
        for z1 in range(1, Z+2):
            temp_fracs[counter, :, z1-1] = abs((max[:, z1-1]-min[:, z1-1])/median[:, z1-1])
        percent, numbers_only = [], []
        for z1 in range(1, Z+2):
            #find ion abundances at peak abundance
            peak = numpy.max(median[:, z1 - 1])
            idx = numpy.argmax(median[:, z1 - 1])
            min_p = min[idx, z1 - 1]
            max_p = max[idx, z1 - 1]
            val = abs((max_p - min_p) / peak)
            numbers_only.append(val)
            val = numpy.format_float_positional(val*100, precision=2)
            percent.append(val+' %')

        t['+/- '+str(max_error*100)+'%'] = percent
        values['+/- '+str(max_error*100)+'%'] = numbers_only
        counter+=1

    if plot == True:
        plt.xlabel('Fractional error')
        plt.ylabel('Fractional change in peak abundance')
        plt.title('Ion peak abundance sensitivity to error')
        for z1 in z1_list:
            v = [values[col][z1-1] for col in values.colnames[1:]]
            plt.plot(errors, v, label=ions[z1-1], color=clist(z1-1))
        plt.legend(fontsize='xx-small', loc='upper left')
        plt.savefig(element+' peak abund sensitivity.pdf')
        plt.show()
        plt.close('all')

    if makefiles == True:
        for number in range(1, 30):
            fname = element + ' peak frac changes ' + str(number) + '.csv'
            file = pathlib.Path(fname)
            if file.exists():
                continue
            else:
                with open(fname, mode='w') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames = values.colnames)
                    writer.writeheader()
                    for row in values:
                        data = {}
                        for col in values.colnames:
                            data.update({col: row[col]})
                        writer.writerow(data)
                    writer.writerow({values.colnames[0]:' '})
                break
=======
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    file = open('line_diagnostics', 'rb')
    line_diagnostics = pickle.load(file)
    file.close()
    
    temp_bins1, temp_bins2, dens_bins1, dens_bins2, Te_line_ratios, Te_line_ratios_min, \
               Te_line_ratios_max, dens_line_ratios, dens_line_ratios_min, \
               dens_line_ratios_max, ratio, name, label = \
               [line_ratio_diagnostics.get(k) for k in line_ratio_diagnostics]
        
    temp = '%.E' % Decimal(Te) +'K'
    percentage = '{0:g}'.format(delta_r*100)+'%'
    density = 'dens='+str(dens)
    if process == 'A':
        orig_text= ratio + ' A value $\Delta$ $\pm$' + percentage 
    elif process == 'exc':
        orig_text= ratio + ' Direct excitation rate $\Delta$ $\pm$' + percentage

    #plot emissivity versus temperature
    text = orig_text + ', ' + density
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    ax1.set_xlabel('Temperature in K', fontsize=18)
    ax1.set_ylabel('Line Ratio of Emissivities', fontsize=18)
    anchored_text = AnchoredText(text, loc='lower right', frameon=False)
    ax1.add_artist(anchored_text)
    ax1.semilogx(temp_bins1, Te_line_ratios, label='Original line ratio')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_temp.pdf')
    ax1.fill_between(temp_bins1, Te_line_ratios_min, Te_line_ratios_max, alpha=0.5, color='g', \
                     label="Range of line ratios")
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_temp_range.pdf')
    plt.show()
    
    #plot emissivity versus density
    text = orig_text + ', ' + temp
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.tick_params(axis='x', labelsize=18)
    ax2.tick_params(axis='y', labelsize=18)
    ax2.set_xlabel('Density in cm$^{-3}$', fontsize=18)
    ax1.set_ylabel('Line Ratio of Emissivities', fontsize=18)
    anchored_text = AnchoredText(text, loc='lower right', frameon=False)
    ax2.add_artist(anchored_text)
    ax2.semilogx(dens_bins1, dens_line_ratios, label='Original emissivity')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_dens.pdf')
    ax2.fill_between(dens_bins1, dens_line_ratios_min, dens_line_ratios_max, alpha=0.5, color='g', label='Range of line ratios')
    plt.tight_layout()
    plt.legend()
    plt.savefig(name+'line_diagnostic_dens_range.pdf')
    plt.show()

#wrapper_check_sensitivity has yet to be tested
def wrapper_check_sensitivity(Z, z1, Te, dens, process, delta_r, \
            transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif):
    transition_list = transition
    for transition in transition_list:
        check_sensitivity(Z, z1, Te, dens, process, delta_r, \
            transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif)
        
        
def check_sensitivity(Z, z1, Te, dens, process, delta_r, transition, transition_2=None, \
            all_levels=None, nlev=None, npnts=None, show=None, wavelen=(10,20), corrthresh=10e-5, e_signif=0.0):

    """
    Check emissivity sensitivity for specified element, ion, and transition.
    
    Parameters
    ----------
    Z: int
    The nuclear charge of the element
    
    z1 : int
    ion charge +1
    
    te : float
    temperture (Kelvin)
    
    dens : float
    electron density (cm^-3)
    
    process : str
    specify rate to vary, i.e. 'A' or 'exc'
    
    transition : tuple
    (upper, lower) transition to vary
        i.e. to vary 5->4 A value, transition = (5,4)
    
    delta_r : float
    delta of rate to vary, i.e. delta_r = 0.1 varies rate by +-0.1
    
    npts : int
    number of points to calculate emissivity for
    
    corrthresh : float
    the minimum desired correlation threshold for epsilon, dE/dE_orig
    
    e_signif : float
    the minimum value of the partial derivative of epsilon to rate
        i.e. significance of dE/dR
  
    Returns
    -------
    dict of lambda, transition info, correlation factor (if deemed prominent change, dE/dR, and the orig E)
    
    """
    
    print("Z="+str(Z), "z1="+str(z1), "Te="+str(Te), "dens="+str(dens), "process="+str(process), \
          "delta_r="+str(delta_r), "transition="+str(transition), "transition2="+str(transition_2), \
          "all_levels="+str(all_levels), "nlev="+str(nlev), "npnts="+str(npnts), "show="+str(show),
            "wavelength range="+str(wavelen), "correlation threshold="+str(corrthresh), "epsilon significance="+str(e_signif))
    
    if npnts is None:
        npnts = 3
    elif transition_2 is None:
        transition_2 = []
    elif all_levels is None:
        all_levels = True
    elif show is None:
        show = True
    elif corrthresh is None:
        corrthresh = 10e-5
    elif e_signif is None:
        e_signif = 0.0
    elif wavelen is None:
        wavelen = (10, 20)
        
    if transition_2 is None:    #check sensitivity for a single transition
        inputs, values, transition = set_up(Z, z1, Te, dens, process, delta_r, \
                transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif)
        if process == 'A':
            new_inputs, new_values = vary_a(inputs, values, transition)
        elif process == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, transition)
        table, new_table, inputs, results = get_tables(new_inputs, new_values)
        
        plot_sensitivity(inputs, new_table)

        line_diagnostics = run_line_diagnostics(table, inputs, values, transition)

        plot_line_diagnostics(inputs, line_diagnostics)

        file = open('results_'+process, 'wb')
        pickle.dump([results, line_diagnostics], file)
        file.close()

        plt.show()

    elif transition_2 != None:  #calculate line ratio diagnostics 
        inputs, values, transition = set_up(Z, z1, Te, dens, process, delta_r, \
                transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif)

        inputs_2, values_2, transition_2 = set_up(Z, z1, Te, dens, process, delta_r, \
                transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif)
        
        if process == 'A':
            new_inputs, new_values = vary_a(inputs, values, transition)
            new_inputs_2, new_values_2, = vary_a(inputs_2, values_2, transition_2)
        elif process == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, transition)
            new_inputs_2, new_values_2, = vary_exc(inputs_2, values_2, transition_2)

        table1, new_table1, inputs1, results1 = get_tables(new_inputs, new_values)
        table2, new_table2, inputs2, results2 = get_tables(new_inputs_2, new_values_2)
        
        wrapper_plot_sensitivity(inputs, new_table1, new_table2)

        #line_ratio_diagnostics = wrapper_run_line_diagnostics(table1, inputs, values, table2, inputs_2, values_2)
        #plot_line_ratio_diagnostics(inputs, line_ratio_diagnostics)
        
        
        file = open('results_'+process, 'wb')
        pickle.dump([results1, results2, line_ratio_diagnostics], file)
        file.close()


    elif isinstance(transition, list) == True: #run check_sensitivity on multiple transitions separately \
        #where transition=[] is list of transitions
        wrapper_check_sensitivity(Z, z1, Te, dens, process, delta_r, \
                    transition, transition_2, all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif)
        
    elif show == True:
        plt.show()
             
>>>>>>> d62ccbc3c5cec6ddf9ffad10769e49c748c5abdf
