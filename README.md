# L5 SBAS MOPS Ephmeris Parameter Fitting Scheme #

This matlab code fits the L5 SBAS MOPS ephemeris message parameters to precision orbit data. It also performs fit error analysis and evaluates the message performance. Specificially, this looks at the corner cases that can cause problems with the fitting algorithm convergence. This implements the algorithms outlined in Appendix B of my PhD thesis undertaken in the GPS Research Lab in the Department of Aeronautics and Astronautics at Stanford University: 

[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

This thesis is available at the following link: https://purl.stanford.edu/dc409wn9227

This also contains functions which create the scripts necessary to run the NASA GMAT software used in creating the high precision orbit data used to fit the ephemeris parameters and perform error analysis. 


## How to Use ##

'MAIN_fit_ephemeris.m' fits the L5 SBAS MOPS ephemeris parameters to the precision orbit ephemeris data produced by the NASA GMAT ( https://gmat.gsfc.nasa.gov/ ) and does message performance analysis. 

'MAIN_make_GMAT_script.m' creates the scripts required to run the NASA GMAT and obtain the orbit data needed for analysis (though some test data is included here already).


## More Info ## 

For more info, please refer to the following publications, all available on the Stanford GPS Lab website ( https://gps.stanford.edu/publications/all-publications ): 

[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

[2]	T. G. R. Reid, T. Walter, P. K. Enge, and T. Sakai, "Orbital representations for the next generation of satellite-based augmentation systems," GPS Solutions, vol. 20, no. 4, pp. 737-750, 2016.

[3]	T. Reid, T. Walter, and P. Enge, "Qualifying an L5 SBAS MOPS Ephemeris Message to Support Multiple Orbit Classes," in Proceedings of the 26th International Technical Meeting of the Satellite Division of the Institute of Navigation (ION GNSS+ 2013), Nashville, TN, 2013.

[4]	T. Reid, T. Walter, and P. Enge, "L1/L5 SBAS MOPS Ephemeris Message to Support Multiple Orbit Classes," in Proceedings of the 2013 International Technical Meeting of The Institute of Navigation, San Diego, CA, 2013.
