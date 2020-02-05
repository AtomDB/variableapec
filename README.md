# variableapec
For looking at emission and line ratios as function of varying atomic data.

Varies any fundamental atomic data by a set amount and re-runs the emissivity calculations to identify which lines are sensitive to the parameter(s) changed. By making use of uncertainty estimates within the AtomDB project, the routine allows identification of which processes will be of most impact in analyzing existing and upcoming high resolution data. Routine outputs plots of sensitive emission lines and line diagnostics (emission as a function of temperature and density).

Installation:
============
Requires PyAtomDB package and Python 3.

Usage Examples:
==============
Check the emissivity sensitivity to a 10% change in the A value for 27->1 transition for Fe XVII at 3e6 K and dens = 1:

	import variableapec
	Z, z1, Te, dens, process, delta_r, transition = 26, 17, 3e6, 1, 'A', 0.1, (27,1)
	variableapec.check_sensitivity(Z, z1, Te, dens, process, delta_r, transition)

Check the emissivity sensitivity to a 10% change in the direct excitation rate for 1->23 transition for Fe XVII at 1e6 K and dens = 1:

	import variableapec
	Z, z1, Te, dens, process, delta_r, transition = 26, 17, 3e6, 1, 'A', 0.1, (1,23)
	variableapec.check_sensitivity(Z, z1, Te, dens, process, delta_r, transition)

Check the emissivity sensitivity of Fe XVII line ratio 1->27/1->23 for a 10% change in the direct excitation rate for these transitions at 3e6 K and dens =1:

	import variableapec
	Z, z1, Te, dens, process, delta_r, transition, transition_2 = 26, 17, 3e6, 1, 'exc', 0.1, (1, 27), (1, 23)
	variableapec.check_sensitivity(Z, z1, Te, dens, process, delta_r, transition, transition_2=transition_2)

Specify the range of wavelengths to plot affected emission lines to be 10-50 Angstroms. The default is 10-20 A.
	
	import variableapec
	Z, z1, Te, dens, process, delta_r, transition = 26, 17, 3e6, 1, 'A', 0.1, (1,23)
	variableapec.check_sensitivity(Z, z1, Te, dens, process, delta_r, transition, wavelen=(10, 50))
	
Specify the range of temperatures and densities to run line diagnostics on. Default temperature range is (Te/10, Te * 10) and default density range is (10, 1e17). 
	
	import variableapec
	Z, z1, Te, dens, process, delta_r, transition = 26, 17, 3e6, 1, 'A', 0.1, (1,23)
	variableapec.check_sensitivity(Z, z1, Te, dens, process, delta_r, transition, wavelen=(10, 50), Te_range=(4e5, 10e8), 		dens_range=(10e12, 10e25))
	
Check sensitivity for multiple lines individually with 20% change in direct excitation rate.
	
	import variableapec
	Z, z1, Te, dens, process, delta_r = 26, 17, 3e6, 1, 'exc', 0.2
	list=[(1,23) (1,27), (1,25)]
	variableapec.wrapper_check_sensitivity(Z, z1, Te, dens, process, delta_r, list)
	
Calculate G or R ratio with 10% uncertainty on direct excitation rate.
	
	import variableapec
	Z, z1, Te, dens, process, delta_r = 8, 7, 1e6, 1, 'exc', 0.1
	variableapec.four_line_diagnostic(Z, z1, Te, dens, process, delta_r)
	variableapec.g_ratio(Z, z1, process)
	variableapec.r_ratio(Z, z1, process)
	
Plot the sensitive emission lines of multiple transitions of an ion.
	
	import variableapec
	Z, z1, Te, dens, delta_r = 26, 17, 3e6, 1, 0.1
	A_lines = {'3C':(27,1), '3D':(23,1), '3E': (17,1), 'M2': (2,1), '3G':(3,1), '3F':(5,1)}
	exc_lines = {'3C':(1,27), '3D':(1,23), '3E': (1,17), 'M2': (1,2), '3G':(1,3), '3F':(1,5)}
	variableapec.plot_multiple_sensitivity(Z, z1, Te, dens, delta_r, A_lines, exc_lines)

	
	
Output:
=========
If successfully run, check_sensitivity() will print table of lines affected sorted by the significance of the change in emissivity. This table has columns of the original emissivity values and several partial derivatives in emissivity. The routine will also produce a scatter plot of any lines sensitive to the parameter(s) changed (see documentation for default values), a plot of emissivity as a function of density, and a plot of emissivity as a function of temperature. If two transitions are specified, the routine will produce these plots for each transition individually but also plot the line ratio as a function of density and temperature. A plot of emissivity from ionization, excitation, and recombination versus temperature will be produced for the specified transition to show the temperatures that excitation dominates. 

Future Plans/Potential:
=================
- Consider auto-ionization and recalculate ionization balances
- Consider resonance scattering
