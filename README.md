# variableapec
For looking at emission and line ratios as function of atomic data uncertainties.

Varies any fundamental atomic data by a set amount and re-runs the emissivity calculations to identify which lines are sensitive to the parameter(s) changed. By making use of uncertainty estimates within the AtomDB project, the routine allows identification of which processes will be of most impact in analyzing existing and upcoming high resolution data. Routines included in package: calculating new emissivities for all lines of an ion after adding uncertainty to 1 line, calculating sensitivity of one line to multiple magnitudes of uncertainties in another line, varying all lines of a line ratios (classical line ratio, g and r ratios, blended line ratios) and recalculating over range of temperatures and densities, line diagnostics over a range of temperatures and densities, partial derivatives of all lines (dE/dR and dE/E) due to changes in one or more rates, finding lines affected by more than 2% due to rate variations. 

Installation:
============
Requires PyAtomDB package and Python 3.


Usage Examples:
==============
Check new line emissivities from a 10% change in 27->1 A value for Fe XVII at 3e6 K and dens = 1:

	import variableapec
	Z, z1, up, lo, Te, dens, vary, delta_r = 26, 17, 27, 1, 3e6, 1, 'A', 0.1
	new_emiss = variableapec.get_all_new_emiss(Z, z1, up, lo, Te, dens, vary, delta_r)

Check new line emissivities from a 10% change in 23->1 direct exc rate Fe XVII at 1e6 K and dens = 1:

	import variableapec
	Z, z1, up, lo, Te, dens, vary, delta_r = 26, 17, 23, 1, 3e6, 1, 'exc', 0.1
	new_emiss = variableapec.get_all_new_emiss(Z, z1, up, lo, Te, dens, vary, delta_r)

Check sensitivity of Fe XVII 17->6 line to 3 uncertainties in 17->1 direct exc rate over 4 temperatures at dens= 1e5
	
	import variableapec
	Z, z1, up, lo, vary, errors, trans_list, temps, dens=5 = 26, 17, 17, 6, 'exc', [0.1, 0.2, 0.3], [(17,1)], [136,2e6,3e6,4e6,5e6]
	variableapec.line_sensitivity(Z, z1, up, lo, vary, errors, trans_list, temps=temps, dens=dens)
	
Run line diagnostics for O VII 2-> 1 transition (no plotting) for 8 densities over default range with 10% uncertainty on A value. 
	
	import variableapec
	Z, z1, up, lo, Te, dens, vary, delta_r, dens_range = 8, 7, 2, 1, 3e6, 1, 'A', 0.1, -1
	line_ratios = variableapec.line_diagnostics(Z, z1, up, lo, Te, dens, vary, delta_r, dens_range=dens_range, num=8, plot=False)

Run and plot line ratio diganostics for Fe XVII 3C/3D (27->1/23->1) ratio with 10% uncertainty on direct excitation rates for 10 temperatures over default temperature range. ***Specify {} or -1 for default dens and Te ranges, only provide input for either Te_range or dens_range if you want temperature or density diagnostics (default is both).***

	import variableapec
	Z, z1, up, lo, up2, lo2, Te, dens, vary, delta_r, Te_range = 26, 17, 27, 1, 23, 1, 3e6, 1, 'exc', 0.10, -1
	line_ratios = variableapec.line_ratio_diagnostics(Z, z1, 27, 1, 23, 1, Te, dens, vary, delta_r, Te_range=Te_range, num=10, plot=True)
	
Calculate and plot O VII G ratio with 10% uncertainty on direct excitation rate over 10 temperatures between 1e6-3e6K.
	
	import variableapec
	Z, z1, Te, dens, vary, delta_r, Te_range, num = 8, 7, 1e6, 1, 'exc', 0.1, (1e6,3e6), 10
	g_ratio = variableapec.g_ratio(Z, z1, Te, dens, vary, delta_r, Te_range=Te_range, num=num, need_data=True, plot=True)
	
Calculate Fe XXV R ratio (no plot) with 10% uncertainty on A value over 5 densities between (1, 1e16).
	
	import variableapec
	Z, z1, Te, dens, vary, delta_r, dens_range, num = 26, 25, 1e6, 1, 'A', 0.1, (1,1e16), 5
	r_ratio = variableapec.r_ratio(Z, z1, Te, dens, vary, delta_r, dens_range=dens_range, num=num, need_data=True, plot=False)

Generate error analysis PDF for O VI 2->1 transition at 1e6K and dens=1 with max uncertainty of 15% and a linewidth of 5eV for blended lines, output file named 'O7.pdf'.

	import variableapec
	Z, z1, up, lo, Te, dens, delta_r = 8, 7, 2, 1, 1e6, 1, 0.15
	variableapec.error_analysis(Z, z1, up, lo, Te, dens, delta_r, linewidth=5, filename='O7')
	
Future Plans/Potential:
=================
- Consider auto-ionization and recalculate ionization balances
- Consider resonance scattering
