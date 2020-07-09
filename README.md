# variableapec
For looking at emission and line ratios as function of atomic data uncertainties.

Varies any fundamental atomic data by a set amount and re-runs the emissivity calculations to identify which lines are sensitive to the parameter(s) changed. By making use of uncertainty estimates within the AtomDB project, the routine allows identification of which processes will be of most impact in analyzing existing and upcoming high resolution data. Routines included in package: calculating new emissivities for all lines of an ion after adding uncertainty to 1 line, calculating sensitivity of one line to multiple magnitudes of uncertainties in another line, varying all lines of a line ratios (classical line ratio, g and r ratios, blended line ratios) and recalculating over range of temperatures and densities, line diagnostics over a range of temperatures and densities, partial derivatives of all lines (dE/dR and dE/E) due to changes in one or more rates, finding lines affected by more than 2% due to rate variations. 

Installation:
============
Requires PyAtomDB package and Python 3.


Usage Examples:
==============
Check new line emissivities from a 10% change in 27->1 A value for Fe XVII at 3e6 K and dens = 1:

	Z, z1, up, lo, Te, dens, vary, delta_r = 26, 17, 27, 1, 3e6, 1, 'A', 0.1
	new_emiss = variableapec.get_all_new_emiss(Z, z1, up, lo, Te, dens, vary, delta_r)

Check new line emissivities from a 10% change in 23->1 direct exc rate Fe XVII at 1e6 K and dens = 1:

	Z, z1, up, lo, Te, dens, vary, delta_r = 26, 17, 23, 1, 3e6, 1, 'exc', 0.1
	new_emiss = variableapec.get_all_new_emiss(Z, z1, up, lo, Te, dens, vary, delta_r)

Check sensitivity of Fe XVII 17->6 line to 3 uncertainties in 17->1 direct exc rate over 4 temperatures at dens = 1e5 cm^-3
	
	Z, z1, up, lo, vary, errors, trans_list, temps, dens, temps = 26, 17, 17, 6, 'exc', [0.1, 0.2, 0.3], [(17,1)], 1e5, [1e6,2e6,3e6,4e6,5e6]
	variableapec.line_sensitivity(Z, z1, up, lo, vary, errors, trans_list, temps=temps, dens=dens)
	
Run line diagnostics for O VII 2->1 transition (no plotting) for 8 densities over default range with 10% uncertainty on A value. 
	
	Z, z1, up, lo, Te, dens, vary, delta_r = 8, 7, 2, 1, 3e6, 1, 'A', 0.1 
	line_ratios = variableapec.line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, dens_range={}, num=8, plot=False)

Run and plot line ratio diganostics for Fe XVII 3C/3D (27->1/23->1) ratio with 10% uncertainty on direct excitation rates for 10 temperatures between 1e4 and 1e9 K. ***Specify {} or -1 for default dens and Te ranges, only provide input for either Te_range or dens_range if you want temperature or density diagnostics (default is both).***

	Z, z1, up, lo, up2, lo2, Te, dens, vary, delta_r, Te_range = 26, 17, 27, 1, 23, 1, 3e6, 1, 'exc', 0.10, [1e4, 1e9]
	line_ratios = variableapec.line_ratio_diagnostics(Z, z1, up, lo, up2, lo2, Te, dens, vary, delta_r, Te_range=Te_range, num=10, plot=True)
	
Calculate and plot O VII G ratio with 10% uncertainty on direct excitation rate over 10 temperatures between 1e6-3e6K.
	
	Z, z1, Te, dens, vary, delta_r, Te_range, num = 8, 7, 1e6, 1, 'exc', 0.1, (1e6,3e6), 10
	g_ratio = variableapec.g_ratio(Z, z1, Te, dens, vary, delta_r, Te_range=Te_range, num=num, need_data=True, plot=True)
	
Calculate Fe XXV R ratio (no plot) with 10% uncertainty on A value over 5 densities between (1, 1e16).
	
	Z, z1, Te, dens, vary, delta_r, dens_range, num = 26, 25, 1e6, 1, 'A', 0.1, (1,1e16), 5
	r_ratio = variableapec.r_ratio(Z, z1, Te, dens, vary, delta_r, dens_range=dens_range, num=num, need_data=True, plot=False)
	
Calculate and plot Fe XIII blended line ratio with 15% error on direct excitation rates at Te = 200 eV and 15 densities between 1e8 and 1e13 cm^-3.
	
	Z, z1, Te, dens, delta_r = 26, 13, 2.3e6, 10e11, 0.10       
	lines = [(25,3), (24,3), (20,1)]         
	denom = 1             #i.e blended line ratio = (25,3)+(24,3)/(20,1)
	dens_range, num = (1e8, 1e13), 15
	variableapec.blended_line_ratio(Z, z1, Te, dens, 'exc', delta_r, lines, denom=denom, type='dens', num=num, dens_range=dens_range, plot=True)
	
Plot Fe XXV R ratio from 3 different temperatures using Matplotlib plasma colormap with supplied legend labels.
	
	fnames = ['O VII_four_line_dens_2.fits', 'O VII_four_line_dens_5.fits', 'O VII_four_line_dens_3.fits']
	labels = ["Te = 5e5 K", "Te = 3e6 K", "Te = 1e7 K"]
	variableapec.plot_ratio(fnames, ratio='r', labels=labels, cmap='plasma')

Plot Fe XXII density-sensitive line ratio from .fits file created by line_ratio_diagnostics()
	
	variableapec.plot_ratio("Fe XXII_21-1-22-2 ratio dens 1.fits")
	
Plot Fe XXII density-sensitive line ratio from 2 .fits files created by line_ratio_diagnostics() for 10% and 30% errors.

	variableapec.plot_ratio(["Fe XXII_21-1-22-2 ratio dens 1.fits", "Fe XXII_21-1-22-2 ratio dens 2.fits"], labels=['30%', '10%'])

Generate error analysis PDF for O VI 2->1 transition at 1e6K and dens=1 with max uncertainty of 15% and a linewidth of 5eV for blended lines, output file named 'O7.pdf'.

	Z, z1, up, lo, Te, dens, delta_r = 8, 7, 2, 1, 1e6, 1, 0.15
	variableapec.error_analysis(Z, z1, up, lo, Te, dens, delta_r, linewidth=5, filename='O7')
	
Plot new O CSD from varying O 6+ ionization rate by 15%.
	
	Z, z1, varyir, delta_r = 8, 7, 'i', 0.15
	median, min, max = variableapec.vary_csd(Z, z1, varyir, delta_r)
	
Plot new Fe CSD from varying all rate coefficients via Monte Carlo calculations using ion specific errors over temperatures 1e6 to 1e9 K.
	
	max_ionize = {26: 0.15, 25: 0.14, 24: 0.13, 23: 0.265, 22: 0.4, 21: 0.4, 20: 0.4, 19: 0.4, 18: 0.15, 17: 0.15, 16: 0.2, 15: 0.13, 14: 0.06, 13: 0.1, 			   12: 0.14, 10: 0.15, 9: 0.16, 8: 0.16, 7: 0.16, 6: 0.16, 5: 0.16, 4: 0.16, 3: 0.16, 2: 0.16, 1: 0.16}
    	max_recomb = {26: 0.2, 25: 0.2, 24: 0.2, 23: 0.2, 22: 0.2, 21: 0.2, 20: 0.2, 19: 0.2, 18: 0.2, 17: 0.2, 16: 0.2, 15: 0.18, 14: 0.18, 13: 0.18, 
		      12: 0.24, 11: 0.3, 0: 0.3, 9: 0.3, 8: 0.26, 7: 0.26, 6: 0.26, 5: 0.26, 4: 0.26, 3: 0.26, 2: 0.26, 1: 0.26}
	Z = 26
	median, min, max = variableapec.monte_carlo_csd(Z, max_ionize, max_recomb, Te_range=[1e6,1e9], plot=True)
	
Plot new O CSD from varying all rate coefficients via Monte Carlo calculations over temperatures 1e5 to 1e7 K and a single max ionization/recombination error of 15%.
	
	Z, max_ionize = 8, 0.15
	median, min, max = variableapec.monte_carlo_csd(Z, max_ionize, Te_range=[1e5,1e7], plot=True)
	
Plot ratio of O VIII/O VII after perturbing CSD via Monte Carlo calculations over temperatures between 1e4 and 1e8.

	Z, z1, z2 = 8, 8, 7
	max_ionize = {8: 0.1, 7: 0.2, 6: 0.05, 5: 0.05, 4: 0.05, 3: 0.05, 2: 0.05, 1: 0.05}
	max_recomb = {8: 0.18, 7: 0.2, 6: 0.16, 5: 0.3, 4: 0.3, 3: 0.3, 2: 0.3, 1: 0.3}
	ion_ratio(Z, z1, z2, max_ionize, max_recomb, Te_range=[1e4,1e8])
	
Plot fractional change in peak abundance for all O ions over a range of errors, write values to csv file.

	Z, errors = 8, [0.1, 0.2, 0.3, 0.4, 0.5]
	variableapec.peak_frac_sensitivity(Z, errors, z1_list={}, makefiles=True, plot=True)
	
Plot fractional change in peak abundance for select highly charged Fe ions over a range of errors, don't write values to csv file.
	
	Z, errors = 26, [0.1, 0.2, 0.3, 0.4, 0.5]
	variableapec.peak_frac_sensitivity(Z, errors, z1_list=[16,17,18,19,20], makefiles=False, plot=True)

Plot sensitive lines between 10-20 A to 30% perturbations in the A value of key Fe XVII emission lines at 3e6 K, with a minimum emissivity of 1e-20 ph/cm^-3 and minimum change in emissivity of 1e-4. 

	Z, z1, Te, dens, delta_r, vary = 26, 17, 3e6, 1, 0.30, 'A'
	lines = {'3C':(27,1), '3D':(23,1), '3E': (17,1), 'M2': (2,1), '3G':(3,1), '3F':(5,1)}
	variableapec.plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, vary, lines, wavelen=[10,20], corrthresh=1e-4, e_signif=1e-20)
	
	### can also specify lines as a list of (up, lo) transitions if there is no common name, x-axis tick labels will just be transition ####
	
	Z, z1, Te, dens, delta_r, vary = 26, 17, 3e6, 1, 0.30, 'A'
	lines = [(27,1), (23,1), (17,1), (2,1), (3,1), (5,1)]
	variableapec.plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, vary, lines, wavelen=[10,20], corrthresh=1e-4, e_signif=1e-20)
	
Plot fractional change in abundance of O 5+ ion (full CSD curve) over temperatures 1e4 and 1e9 K and multiple errors. 

	Z, z1, errors = 8, 6, [0.1, 0.2, 0.3, 0.4]
	variableapec.ion_sensitivity(Z, z1, errors, Te_range=[1e4,1e9])
	
Future Plans/Potential:
=================
- Consider auto-ionization and recalculate ionization balances
- Consider resonance scattering
