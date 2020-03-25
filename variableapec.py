"""
This module contains methods for looking at emission and line ratios
as a function of varying atomic data from the AtomDB files. Requires
PyAtomdB and Python 3."""

# Keri Heuer
# Version 1.5, March 16, 2020

import matplotlib.pyplot as plt
import pyatomdb, numpy, pickle, pathlib, csv
from astropy.table import Table, Column
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from decimal import Decimal

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

def set_up(Z, z1, Te, dens, process, delta_r, transition, transition_2, \
           npnts, wavelen, Te_range, dens_range, corrthresh, e_signif):
    
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
    inputs.update( {'Z': Z, 'z1': z1, 'Te': Te, 'dens': dens, 'process': process, 'delta_r': delta_r, \
                    'transition': transition, 'transition_2': transition_2, \
                    'npnts': npnts, 'wavelen': wavelen, 'Te_range': Te_range, 'dens_range': dens_range,\
                    'corrthresh': corrthresh, 'e_signif': e_signif} )
    values = {}
    values.update( {'matrix': matrix, 'B': B, 'in_range': in_range, 'linelist': linelist, 'table': table, 'new_table': new_table} )
    return inputs, values, transition

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

    # idx1 = numpy.where((0 < new_table['|dE/dE_orig|']) & (new_table['|dE/dE_orig|'] < 0.02))[0]
    # idx2 = numpy.where(0.98 < new_table['|dE/dE_orig|'])[0] #< 1.02
    # print(idx1, idx2)
    # print(type(idx1), type(idx2))
    # idxs=[]
    # idxs.append(idx1)
    # idxs.append(idx2)
    # observable = new_table[idxs]
    # for number in range(1,10):
    #     file = pathlib.Path('observable partial deriv' + element + str(z1) + process + str(number)+'.csv')
    #     if file.exists():
    #         continue
    #     else:
    #         with open('observable partial deriv' + element + str(z1) + process + str(number)+'.csv', mode='w') as csv_file:
    #             fieldnames = observable.colnames
    #             writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    #             writer.writeheader()
    #             for row in observable:
    #                 data = {}
    #                 for col in fieldnames:
    #                     data.update({col: row[col]})
    #                 writer.writerow(data)

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
            #print("Temperature bins:", temp_bins)
            print("Current temperature is:", temp_Te)
            Te_inputs, Te_values, transition = set_up(Z, z1, temp_Te, dens, process, delta_r, which_transition, transition_2, \
               npnts, wavelen, Te_range, dens_range, corrthresh, e_signif)
            if process == 'A':
                Te_new_inputs, Te_new_values = vary_a(Te_inputs, Te_values, which_transition)
                Te_table, Te_new_table, Te_inputs, Te_results = get_tables(Te_new_inputs, Te_new_values)
                #line_index = numpy.where((Te_table['Upper'] == which_transition[0]) & (Te_table['Lower'] == which_transition[1]))[0]
                # for a, b, c, d, e in zip(Te_table['Upper'], Te_table['Lower'], Te_table['Epsilon_orig'], Te_table['Epsilon_1'],
                #                          Te_table['Epsilon_' + str(npnts)]):
                #     if (a, b) == which_transition:
                #         if c == 0.0: c = 1e-40
                #         if d == 0.0: d = 1e-40
                #         if e == 0.0: e = 1e-40
                #         Te_eps_orig.append(c)
                #         Te_eps_min.append(d)
                #         Te_eps_max.append(e)

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

                #line_index = numpy.where((Te_table['Upper'] == which_transition[1]) & (Te_table['Lower'] == which_transition[0]))[0]
                #Te_table = Te_table[(Te_table['Upper'] == which_transition[1]) & (Te_table['Lower'] == which_transition[0])]

                # Te_eps_orig.append(Te_table['Epsilon_orig'][0])
                # Te_eps_min.append(Te_table['Epsilon_1'][0])
                # Te_eps_max.append(Te_table['Epsilon_'+str(npnts)][0])

                # for a, b, c, d, e in zip(Te_table['Upper'], Te_table['Lower'], Te_table['Epsilon_orig'],
                #                          Te_table['Epsilon_1'], Te_table['Epsilon_' + str(npnts)]):
                #     if (b, a) == which_transition:
                #         if c == 0.0: c = 1e-40
                #         if d == 0.0: d = 1e-40
                #         if e == 0.0: e = 1e-40
                #         Te_eps_orig.append(c)
                #         Te_eps_min.append(d)
                #         Te_eps_max.append(e)

                for x in Te_table:
                    if (x['Lower'], x['Upper']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        Te_eps_orig.append(x['Epsilon_orig'])
                        Te_eps_min.append(x['Epsilon_1'])
                        Te_eps_max.append(x['Epsilon_'+str(npnts)])

                print(Te_table)
                #if Te_table['Epsilon_orig'] == 0.0: Te_table['Epsilon_orig'] = 1e-40

            print('\n the arrays should be\n', Te_eps_orig, '\n', Te_eps_min, '\n', Te_eps_max)

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
            #print("dens bins:", dens_bins)
            #print("Current density is:", temp_dens)
            dens_inputs, dens_values, transition = set_up(Z, z1, Te, temp_dens, process, delta_r, which_transition, transition_2, \
               npnts, wavelen, Te_range, dens_range, corrthresh, e_signif)
            if process == 'A':
                dens_new_inputs, dens_new_values = vary_a(dens_inputs, dens_values, which_transition)
                dens_table, dens_new_table, dens_inputs, dens_results = get_tables(dens_new_inputs, dens_new_values)
                #line_index = numpy.where((dens_table['Upper'] == transition[0]) & (dens_table['Lower'] == transition[1]))[0]
                # for a, b, c, d, e in zip(dens_table['Upper'], dens_table['Lower'], dens_table['Epsilon_orig'], dens_table['Epsilon_1'],
                #                          dens_table['Epsilon_' + str(npnts)]):
                #     if (a, b) == which_transition:
                #         if c == 0.0: c = 1e-40
                #         if d == 0.0: d = 1e-40
                #         if e == 0.0: e = 1e-40
                #         print(a,b,c,d,e)
                #         dens_eps_orig.append(c)
                #         dens_eps_min.append(d)
                #         dens_eps_max.append(e)

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
                #line_index = numpy.where((dens_table['Upper'] == transition[1]) & (dens_table['Lower'] == transition[0]))[0]
                # for a, b, c, d, e in zip(dens_table['Upper'], dens_table['Lower'], dens_table['Epsilon_orig'], dens_table['Epsilon_1'],
                #                          dens_table['Epsilon_' + str(npnts)]):
                #     if (b, a) == which_transition:
                #         if c == 0.0: c = 1e-40
                #         if d == 0.0: d = 1e-40
                #         if e == 0.0: e = 1e-40
                #

                for x in dens_table:
                    if (x['Lower'], x['Upper']) == which_transition:
                        if x['Epsilon_orig'] == 0.0: x['Epsilon_orig'] = 1e-40
                        if x['Epsilon_1'] == 0.0: x['Epsilon_1'] = 1e-40
                        if x['Epsilon_'+str(npnts)] == 0.0: x['Epsilon_'+str(npnts)] = 1e-40
                        dens_eps_orig.append(x['Epsilon_orig']/ temp_dens)
                        dens_eps_min.append(x['Epsilon_1']/temp_dens)
                        dens_eps_max.append(x['Epsilon_'+str(npnts)]/temp_dens)

            counter += 1
            print(str(dens_num-counter), 'densities left for', element, ion)
            #print('\nCurrent dens arrays:\n', dens_eps_orig, '\n', dens_eps_min, '\n', dens_eps_max, '\n')

        # dens_table = dens_table[line_index]
        # dens_eps_orig.append(dens_table['Epsilon_orig'][0]/temp_dens)
        # dens_eps_min.append(dens_table['Epsilon_1'][0]/temp_dens)
        # dens_eps_max.append(dens_table['Epsilon_'+str(npnts)][0]/temp_dens)
    
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

            inputs, values, transition = set_up(Z, z1, Te, dens, process, delta_r,
             transition, transition_2=None, npnts=2, wavelen=wavelen, Te_range=Te_range,
                    dens_range=dens_range,corrthresh=corrthresh, e_signif=0.0)
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

    print("Z=" + str(Z), "z1=" + str(z1), "Te=" + str(Te), "dens=" + str(dens), "process=" + str(process), \
          "delta_r=" + str(delta_r), "transition=" + str(transition), "transition2=" + str(transition_2), \
          "npnts=" + str(npnts), "wavelength range=" + str(wavelen), "temp range=" + str(Te_range), \
          "dens range=" + str(dens_range), "correlation threshold=" + str(corrthresh),
          "epsilon significance=" + str(e_signif))

    if transition_2 =={}:    #check sensitivity for a single transition
        inputs, values, transition = set_up(Z, z1, Te, dens, process, delta_r, \
                transition, transition_2, npnts, wavelen, Te_range, dens_range, corrthresh, e_signif)
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
        which_transition = transition
        print("Doing calculations for", which_transition)
        inputs, values, transition = set_up(Z, z1, Te, dens, process, delta_r, \
                which_transition, transition_2, npnts, wavelen, Te_range, dens_range, corrthresh, e_signif)
        if process == 'A':
            new_inputs, new_values = vary_a(inputs, values, which_transition)
        elif process == 'exc':
            new_inputs, new_values = vary_exc(inputs, values, which_transition)

        which_transition = transition_2
        transition_2 = None
        print("Doing calculations for", which_transition, "and transition 2 is", transition_2)
        inputs_2, values_2, transition_2 = set_up(Z, z1, Te, dens, process, delta_r, \
                which_transition, transition_2, npnts, wavelen, Te_range, dens_range, corrthresh, e_signif)
        if process == 'A':
            new_inputs_2, new_values_2, = vary_a(inputs_2, values_2, which_transition)
        elif process == 'exc':
            new_inputs_2, new_values_2, = vary_exc(inputs_2, values_2, which_transition)

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
        
        # file = open('results_'+process, 'wb')
        # pickle.dump([results1, results2, line_ratio_diagnostics], file)
        # file.close()

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