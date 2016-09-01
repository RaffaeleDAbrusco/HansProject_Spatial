## Collection of functions that produce random homogeneous distribution of 
## simulated coordinates in a fixed area of the sky.
##
## Written by R. D'Abrusco
##
## Last modified: 31/3/16.

generate_points_homogeneous <- function(denspoints, ra_range, dec_range) {
	# This function generates the sky coordinates of points homogeneously randomly distributed
	# with a given fixed density, in a given rectangular area of the sky, specified by 
	# the range along the right ascension and declination axes.
	# Arguments:
	#	denspoints 	->	fixed density of points [number of sources per arcmin squared].
	# 	ra_range	->	range over which right ascension values have to be generated [deg].
	#	dec_range	->	range over which declination values have to be generated [deg].
	
	# Generating the random R.A. and Dec. coordinates within the 
	# intervals specified by ra_range and dec_range, with the fixed
	# given density.
	area_points = (max(ra_range) - min(ra_range))*(max(dec_range) - min(dec_range))*3600
	numpoints_all = ceiling(area_points*denspoints)
	sim_xy_all = array(0, dim = c(numpoints_all, 2))
	sim_xy_all[, 1] = runif(numpoints_all, min = min(ra_range), max = max(ra_range))
	sim_xy_all[, 2] = runif(numpoints_all, min = min(dec_range), max = max(dec_range))
			  			   
	# Returning the structure containing the output.	
	return(sim_xy_all)
	}