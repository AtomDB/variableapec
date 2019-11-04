# variableapec
For looking at emission and line ratios as function of varying atomic data.

Varies any fundamental atomic data by a set amount and re-runs the emissivity calculations to identify which lines are sensitive to the parameter(s) changed. By making use of uncertainty estimates within the AtomDB project, the routine allows identification of which processes will be of most impact in analyzing existing and upcoming high resolution data. 

============
Installation:
============
Requires PyAtomDB package and Python 3.

==============
Usage Examples:
==============
Check the emissivity sensitivity to a 10$\%$ change in the A value for 27->1 transition for Fe XVII at 3e6 K and dens = 1:
	import variableapec
	variableapec.check_sensitivity(26, 17, 3e6, 1, 'A', transition=(27,1), delta_r=0.1)

Check the emissivity sensitivity to a 10$\%$ change in the direct excitation rate for 1->23 transition for Fe XVII at 3e6 K and dens = 1:
	import variableapec
	variableapec.check_sensitivity(26, 17, 3e6, 1, 'exc', transition=(1,23), delta_r=0.1)

Check the emissivity sensitivity of Fe XVII line ratio 1->27/1->23 for a 10$\%$ change in the direct excitation rate for these transitions at 3e6 K and dens =1:
	import variableapec
	variableapec.check_sensitivity(26, 17, 3e6, 1, 'exc', transition=(1, 27), transition_2=(1, 23), delta_r=0.1)

Specify the range of wavelengths to plot affected emission lines to be 10 to 50 Angstroms:
	import variableapec
	variableapec.check_sensitivity(26, 17, 3e6, 1, 'exc', transition=(1,23), wavelen=[10, 20]
	
=========
Output:
=========
If successfully run, check_sensitivity() will print table of lines affected sorted by the significance of the change in emissivity. This table has columns of the original emissivity values and several partial derivatives in emissivity. The routine will also produce a scatter plot of any lines sensitive to the parameter(s) changed (for a default range of wavelengths from 10-20 Angstroms), a plot of emissivity as a function of density, and a plot of emissivity as a function of temperature. If two transitions are specified, the routine will produce these plots for each transition individually but also plot the line ratio as a function of density and temperature. A plot of emissivity from ionization, excitation, and recombination versus temperature will be produced for the specified transition to show the temperatures that excitation dominates. 

=================
Future Plans/Potential:
=================
- Consider auto-ionization and recalculate ionization balances
- Consider resonance scattering
- User input range of temperatures/densities to run line diagnostics on
