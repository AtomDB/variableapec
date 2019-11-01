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
    
    Set_up() gets original rates, radiative collision matrix, level populations, line lists, 
    and intensities and returns these in a dictionary called values. Creates sensitivity tables.
    
    Returns inputs in dictionary called inputs, returns the variable transition."""
    
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
    for i in range(len(init)):
        x = final[i]
        y = init[i]
        matrix[x][y] += rates[i]
        
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
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the A value for by the specified delta_r
    in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
    all_levels, nlev, npnts, show, wavelen, corrthresh, e_signif = [inputs.get(k) for k in inputs]
    
    matrix, B, in_range, linelist, table, new_table = [values.get(k) for k in values]
    
    print("****************************************************")
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
        
    index=1
    for x in new_a:
        #update LA data for specified transition
        in_range['EINSTEIN_A'][a_index] = x
        
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
        
        #update sensitivity table 
        new_col = Column(name='Epsilon_'+str(index), data = new_linelist['epsilon'], unit = frac+' A')
        table.add_columns([new_col])
            
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
    
    """ Inputs and values are dictionaries outputted by set_up() subroutine. Transition is a tuple
    outputted by set_up() and dictates which transition to vary the excitation rate for by the specified
    delta_r in inputs.
    
    Recalculates the emissivities and updates sensitivity table with these values. Updates the values
    dictionary with tables.
    
    Returns the dictionaries input and values."""
    
    Z, z1, Te, dens, process, delta_r, transition, transition_2, \
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
    table outputted by get_tables().
    
    Plots the lines affected by the parameter(s) changed, including the original transition.
    Currently plotting for wavelengths between 10-20 Angstroms (as set by mask1 and mask2 below).""" 
    
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
    
    #label each point w/ transition
    transition_labels=[]
    for x in cutoff_data:
        if process == 'A':
            transition_name = '{0:g}'.format(x['Upper'])+'->'+'{0:g}'.format(x['Lower'])
        if process == 'exc':
            transition_name = '{0:g}'.format(x['Lower'])+'->'+'{0:g}'.format(x['Upper'])
        transition_labels.append(transition_name)
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
    
    """ Line_diagnostics is a dictionary outputted by run_line_diagnostics. It holds arrays of
    the temperature and density bins, the original, minimum, and maximum emissivities from varying
    temperature and then density, the name of the element and ion, label of which process varied
    and by how much, and the transition.
    
    Plots emissivity as a function of temperature and density for the specified single transition."""

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
    
    """ Table1 and table2 are tables from individually run get_tables() on the two transitions
    specified by the user. Inputs1, values1, inputs2, values2, are  dictionaries holding the
    inputs and the sensitivity table/emissivity values for each transition respectively.
    
    Varies temperature and density separately for each transition and recalculates emissivity.
    
    Returns dictionary line_ratio_diagnostics containing arrays of the temperature and density bins
    for each transition, as well as the original, minimum, and maximum line ratios calculated
    from varying temperature and density independently."""
    
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
    arrays of temperature and density bins, and line ratios (original, min, max) calculated
    from varying temperature and density.
    
    Plots the line ratios of emissivities for specified two transitions as a function of
    temperature and density."""

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
             
