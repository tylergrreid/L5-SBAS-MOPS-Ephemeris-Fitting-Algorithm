# L5 SBAS MOPS Ephmeris Parameter Fitting Scheme #

This matlab code fits the L5 SBAS MOPS ephemeris message parameters to precision orbit data. It also performs fit error analysis and evaluates the message performance. Specificially, this looks at the corner cases that can cause problems with the fitting algorithm convergence. This implements the algorithms outlined in Appendix B of my PhD thesis undertaken in the GPS Research Lab in the Department of Aeronautics and Astronautics at Stanford University: 

[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

    This thesis is available at the following link: https://purl.stanford.edu/dc409wn9227

This also contains functions which create the scripts necessary to run the NASA GMAT software used in creating the high precision orbit data used to fit the ephemeris parameters and perform error analysis. 


## How to Use ##

'MAIN_fit_ephemeris.m' fits the L5 SBAS MOPS ephemeris parameters to the precision orbit ephemeris data produced by the NASA GMAT ( https://gmat.gsfc.nasa.gov/ ) and does a message performance analysis. 

'MAIN_make_GMAT_script.m' creates the scripts required to run the NASA GMAT and obtain the orbit data needed for analysis (though some test data is included here already).



