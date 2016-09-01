## This function returns a description of structures in the residual map as function of the
## radial distance (expressed in two different ways) and the angular distance from the 
## major axis.
##
## Arguments:
##		name_galaxy 		<- 	Name of the galaxy.
##		class_gc			<- 	Class of globular cluster.
##		K_value 			<- 	K value.
##		num_simulations 	<-	Number of simulations for the experiment.
##	    nsteps_smooth		<-	Number of steps along the two axis.
##		sim_type 			<-	Simulation type (always "observed")
## 
## Written by R. D'Abrusco
##
## Last modified: 29/1/15.

radec_to_rtheta <- function(x_regions, 
							y_regions, 
							galaxy_center, 
							galaxy_angle) {
	# This function transforms the coordinates of the vertices of a closed
	# polygon to galactocenctric-angle coordinates.
	x_polygons_external = x_regions
	y_polygons_external = y_regions
	extrapolate_x_polygons_external_list = list()
	extrapolate_y_polygons_external_list = list()
	for (j in 1:(length(x_polygons_external) - 1)) {
		function_segment = approxfun(c(x_polygons_external[j], x_polygons_external[j + 1]), 
							 	     c(y_polygons_external[j], y_polygons_external[j + 1]))
		extrapolate_x_polygons_external_list[[j]] = seq(x_polygons_external[j], x_polygons_external[j + 1], length.out = 50)
		extrapolate_y_polygons_external_list[[j]] = function_segment(extrapolate_x_polygons_external_list[[j]])	
		}	
	extrapolate_x_polygons_external = unlist(extrapolate_x_polygons_external_list)
	extrapolate_y_polygons_external = unlist(extrapolate_y_polygons_external_list)	
	x_extrapolate_polygons_external_aux = extrapolate_x_polygons_external - galaxy_center[1]
	y_extrapolate_polygons_external_aux = extrapolate_y_polygons_external - galaxy_center[2]
	x_extrapolate_polygons_external = x_extrapolate_polygons_external_aux*cos(pi*(galaxy_angle - 90)/180) - 
									  y_extrapolate_polygons_external_aux*sin(pi*(galaxy_angle - 90)/180)
	y_extrapolate_polygons_external = x_extrapolate_polygons_external_aux*sin(pi*(galaxy_angle - 90)/180) + 
									  y_extrapolate_polygons_external_aux*cos(pi*(galaxy_angle - 90)/180)						  
	radius_extrapolate_polygons_external_aux = sqrt(x_extrapolate_polygons_external^2 + y_extrapolate_polygons_external^2)
	theta_extrapolate_polygons_external_aux = atan(abs(y_extrapolate_polygons_external/x_extrapolate_polygons_external))
	theta_extrapolate_polygons_external = vector("numeric", length = length(theta_extrapolate_polygons_external_aux))
	radius_extrapolate_polygons_external = vector("numeric", length = length(theta_extrapolate_polygons_external_aux))
	for (i in 1:length(theta_extrapolate_polygons_external_aux)) {
		if (x_extrapolate_polygons_external[i] > 0 & y_extrapolate_polygons_external[i] > 0) {
			theta_extrapolate_polygons_external[i] = pi + theta_extrapolate_polygons_external_aux[i]
			radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
			} else if (x_extrapolate_polygons_external[i] < 0 & y_extrapolate_polygons_external[i] > 0) {
				theta_extrapolate_polygons_external[i] = 2*pi - theta_extrapolate_polygons_external_aux[i]
				radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
				} else if (x_extrapolate_polygons_external[i] < 0 & y_extrapolate_polygons_external[i] < 0) { 	 
					theta_extrapolate_polygons_external[i] = theta_extrapolate_polygons_external_aux[i]
					radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
					} else {
						theta_extrapolate_polygons_external[i] = pi - theta_extrapolate_polygons_external_aux[i]
						radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
						}	  		
		}	
	theta_extrapolate_polygons_external = 180*theta_extrapolate_polygons_external/pi
	rtheta = cbind(radius_extrapolate_polygons_external, theta_extrapolate_polygons_external)
	extrapolate_xy_polygons_external = cbind(extrapolate_x_polygons_external, extrapolate_y_polygons_external)
	rtheta_extrapolate_xy_polygons_external = list()
	rtheta_extrapolate_xy_polygons_external[[1]] = rtheta
	rtheta_extrapolate_xy_polygons_external[[2]] = extrapolate_xy_polygons_external
	return(rtheta_extrapolate_xy_polygons_external)
	}
			
radec_to_rtheta_boundary <- function(x_regions, 
									 y_regions, 
									 galaxy_center, 
									 galaxy_angle) {
	# This function transforms the coordinates of the vertices of a closed
	# polygon to galactocenctric-angle coordinates.
	x_polygons_external = x_regions
	y_polygons_external = y_regions
	extrapolate_x_polygons_external_list = list()
	extrapolate_y_polygons_external_list = list()
	for (j in 1:(length(x_polygons_external) - 1)) {
		function_segment = approxfun(c(x_polygons_external[j], x_polygons_external[j + 1]), 
							 	     c(y_polygons_external[j], y_polygons_external[j + 1]))
		extrapolate_x_polygons_external_list[[j]] = seq(x_polygons_external[j], x_polygons_external[j + 1], length.out = 50)
		extrapolate_y_polygons_external_list[[j]] = function_segment(extrapolate_x_polygons_external_list[[j]])	
		}	
	extrapolate_x_polygons_external = unlist(extrapolate_x_polygons_external_list)
	extrapolate_y_polygons_external = unlist(extrapolate_y_polygons_external_list)	
	x_extrapolate_polygons_external_aux = extrapolate_x_polygons_external - galaxy_center[1]
	y_extrapolate_polygons_external_aux = extrapolate_y_polygons_external - galaxy_center[2]
	x_extrapolate_polygons_external = x_extrapolate_polygons_external_aux*cos(pi*(galaxy_angle - 90)/180) - 
									  y_extrapolate_polygons_external_aux*sin(pi*(galaxy_angle - 90)/180)
	y_extrapolate_polygons_external = x_extrapolate_polygons_external_aux*sin(pi*(galaxy_angle - 90)/180) + 
									  y_extrapolate_polygons_external_aux*cos(pi*(galaxy_angle - 90)/180)						  
	radius_extrapolate_polygons_external_aux = sqrt(x_extrapolate_polygons_external^2 + y_extrapolate_polygons_external^2)
	theta_extrapolate_polygons_external_aux = atan(abs(y_extrapolate_polygons_external/x_extrapolate_polygons_external))
	theta_extrapolate_polygons_external = vector("numeric", length = length(theta_extrapolate_polygons_external_aux))
	radius_extrapolate_polygons_external = vector("numeric", length = length(theta_extrapolate_polygons_external_aux))
	for (i in 1:length(theta_extrapolate_polygons_external_aux)) {
		if (x_extrapolate_polygons_external[i] > 0 & y_extrapolate_polygons_external[i] > 0) {
			theta_extrapolate_polygons_external[i] = pi + theta_extrapolate_polygons_external_aux[i]
			radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
			} else if (x_extrapolate_polygons_external[i] < 0 & y_extrapolate_polygons_external[i] > 0) {
				theta_extrapolate_polygons_external[i] = 2*pi - theta_extrapolate_polygons_external_aux[i]
				radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
				} else if (x_extrapolate_polygons_external[i] < 0 & y_extrapolate_polygons_external[i] < 0) { 	 
					theta_extrapolate_polygons_external[i] = theta_extrapolate_polygons_external_aux[i]
					radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
					} else {
						theta_extrapolate_polygons_external[i] = pi - theta_extrapolate_polygons_external_aux[i]
						radius_extrapolate_polygons_external[i] = radius_extrapolate_polygons_external_aux[i]
						}	  		
		}	
	theta_extrapolate_polygons_external = 180*theta_extrapolate_polygons_external/pi
	idx_lower_theta = which(theta_extrapolate_polygons_external == min(theta_extrapolate_polygons_external))
	theta_to_min = theta_extrapolate_polygons_external[1:idx_lower_theta]
	theta_from_max = theta_extrapolate_polygons_external[(idx_lower_theta + 1):length(theta_extrapolate_polygons_external)]
	theta_link = c(0, 360)
	theta_new = c(theta_to_min, theta_link, theta_from_max)
	radius_to_min = radius_extrapolate_polygons_external[1:idx_lower_theta]
	radius_from_max = radius_extrapolate_polygons_external[(idx_lower_theta + 1):length(theta_extrapolate_polygons_external)]
	radius_link = c(0, 0)
	radius_new = c(radius_to_min, radius_link, radius_from_max)
	#rtheta = cbind(radius_extrapolate_polygons_external, theta_extrapolate_polygons_external)
	extrapolate_xy_polygons_external = cbind(extrapolate_x_polygons_external, extrapolate_y_polygons_external)
	rtheta = cbind(radius_new, theta_new)
	rtheta_extrapolate_xy_polygons_external = list()
	rtheta_extrapolate_xy_polygons_external[[1]] = rtheta
	rtheta_extrapolate_xy_polygons_external[[2]] = extrapolate_xy_polygons_external
	return(rtheta_extrapolate_xy_polygons_external)
	}				
							
structures_radius <- function(name_galaxy, 
							  class_gc, 
							  K_value = 9, 
							  num_simulations = 500, 
							  nsteps_smooth = 36, 
							  sim_type = "observed", 
							  mode = "all") {
	
	# Loading required libraries.
	source("./raycasting.R")			# Used to determine grid points within 	
										# the overall field region.
	source("./auxiliary.R")				# Loading the auxiliary functions.	
	suppressMessages(library(plotrix, quietly = TRUE))	# Used to plot other shapes and lines.	
	suppressMessages(library(MASS, quietly = TRUE))		# Used to compute 2D densities.				
	cat("		So far so good! Libraries loaded!\n")

	# Setting the paths to the location of the files (data and plots)
	# that will be useful later.
	path_data = "~/Desktop/Spatial/Data/"			# Path of the folder containing the 
													# catalog listing the GCs.
	path_plot = "~/Desktop/Spatial/Plots/"			# Path of the folder containing the
													# plots produced in this program.
	factor_galaxy = list("ngc4472" = 3.5, 
						 "ngc4486" = 2, 
						 "ngc4649" = 3, 
						 "ngc4406" = 4, 
						 "ngc4382" = 3, 
						 "ngc4374" = 3, 
						 "ngc4526" = 3.5,
						 "ngc4365" = 3.5, 
						 "ngc4621" = 3, 
						 "ngc4552" = 3)
	cat("		So far so good! Paths of data and plots set!\n")												
													
	# Loading the file containing the data of the session.
	name_file_session = paste(path_data, name_galaxy, "/Sessions/gc_nsims_", 
							  num_simulations, "_nsteps_", nsteps_smooth, 
							  "_simtype_", sim_type, ".Rdata", sep = "")
	load(name_file_session)
	cat("		So far so good! Image of a previous session imported!\n")
	
	# Evaluating the values of the variables necessary to plot the residual map.
	f = which(k_vector == K_value)
	residuals_matrix_all = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_red = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_blue = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))
	residuals_matrix_highl = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_lowl = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))
	errresiduals_matrix_all_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_red_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_blue_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_highl_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_lowl_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_all = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_red = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_blue = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))			
	errresiduals_matrix_highl = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_lowl = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))			
	residuals_matrix_all = apply(res_all[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_red = apply(res_red[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_blue = apply(res_blue[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_highl = apply(res_highl[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_lowl = apply(res_lowl[[f]], c(2, 3), sum)/num_simulations
	errresiduals_matrix_all_aux = apply(errres_all[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_red_aux = apply(errres_red[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_blue_aux = apply(errres_blue[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_highl_aux = apply(errres_highl[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_lowl_aux = apply(errres_lowl[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_all = sqrt(errresiduals_matrix_all_aux)/num_simulations
	errresiduals_matrix_red = sqrt(errresiduals_matrix_red_aux)/num_simulations
	errresiduals_matrix_blue = sqrt(errresiduals_matrix_blue_aux)/num_simulations
	errresiduals_matrix_highl = sqrt(errresiduals_matrix_highl_aux)/num_simulations
	errresiduals_matrix_lowl = sqrt(errresiduals_matrix_lowl_aux)/num_simulations	
	if (name_galaxy == "ngc4649") {
		num_sigmas_threshold = c(0.5, 1, 2, 3)
		idx_residuals_matrix_all_pos = list(); idx_residuals_matrix_all_neg = list()
		idx_residuals_matrix_red_pos = list(); idx_residuals_matrix_red_neg = list()
		idx_residuals_matrix_blue_pos = list(); idx_residuals_matrix_blue_neg = list()
		idx_residuals_matrix_highl_pos = list(); idx_residuals_matrix_highl_neg = list()
		idx_residuals_matrix_lowl_pos = list(); idx_residuals_matrix_lowl_neg = list()								
		idx_residuals_matrix_all_neg[[f]] = list(); idx_residuals_matrix_all_pos[[f]] = list()
		idx_residuals_matrix_red_pos[[f]] = list(); idx_residuals_matrix_red_neg[[f]] = list()
		idx_residuals_matrix_blue_pos[[f]] = list(); idx_residuals_matrix_blue_neg[[f]] = list()
		idx_residuals_matrix_highl_pos[[f]] = list(); idx_residuals_matrix_highl_neg[[f]] = list()
		idx_residuals_matrix_lowl_pos[[f]] = list(); idx_residuals_matrix_lowl_neg[[f]] = list()
		for (z in seq(length(num_sigmas_threshold))) {
			idx_residuals_matrix_all_pos[[f]][[z]] = which(num_sigmas_all[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_all_neg[[f]][[z]] = which(num_sigmas_all[[f]] <= -num_sigmas_threshold[z])
			idx_residuals_matrix_red_pos[[f]][[z]] = which(num_sigmas_red[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_red_neg[[f]][[z]] = which(num_sigmas_red[[f]] <= -num_sigmas_threshold[z])
			idx_residuals_matrix_blue_pos[[f]][[z]] = which(num_sigmas_blue[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_blue_neg[[f]][[z]] = which(num_sigmas_blue[[f]] <= -num_sigmas_threshold[z])
			idx_residuals_matrix_highl_pos[[f]][[z]] = which(num_sigmas_highl[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_highl_neg[[f]][[z]] = which(num_sigmas_highl[[f]] <= -num_sigmas_threshold[z])
			idx_residuals_matrix_lowl_pos[[f]][[z]] = which(num_sigmas_lowl[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_lowl_neg[[f]][[z]] = which(num_sigmas_lowl[[f]] <= -num_sigmas_threshold[z])
			}	
		idx_within_HSTfield = idx_gridbins_galaxies
		}
	cat("		So far so good! Variables necessary to produce the plot read!\n")		
	
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	if (class_gc == "all") {
		num_sigmas = num_sigmas_all
		idx_residuals_matrix_neg = idx_residuals_matrix_all_neg
		idx_residuals_matrix_pos = idx_residuals_matrix_all_pos
		} else if (class_gc == "red") {
			num_sigmas = num_sigmas_red 
			idx_residuals_matrix_neg = idx_residuals_matrix_red_neg
			idx_residuals_matrix_pos = idx_residuals_matrix_red_pos			
			} else if (class_gc == "blue") {
				num_sigmas = num_sigmas_blue
				idx_residuals_matrix_neg = idx_residuals_matrix_blue_neg
				idx_residuals_matrix_pos = idx_residuals_matrix_blue_pos
				} else if (class_gc == "highl") {
					num_sigmas = num_sigmas_highl
					idx_residuals_matrix_neg = idx_residuals_matrix_highl_neg
					idx_residuals_matrix_pos = idx_residuals_matrix_highl_pos
					} else if (class_gc == "lowl") {
						num_sigmas = num_sigmas_lowl
						idx_residuals_matrix_neg = idx_residuals_matrix_lowl_neg
						idx_residuals_matrix_pos = idx_residuals_matrix_lowl_pos						
						}
	# Calculating the distance of the center of each pixel from the center
	# of the galaxy, and the angular distance from the major axis of the 
	# galaxy. These variables will be used later to calculate properties of 
	# the merged structures.
	dist_grid_dens = sqrt((grid_dens[, 1] - galaxy1_center[1])^2 + 
						  (grid_dens[, 2] - galaxy1_center[2])^2)		
	coords_p1 = return_extreme_axis(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  				    galaxy1_angle)[1, ]
	coords_p2 = return_extreme_axis(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  				    galaxy1_angle)[2, ]			  				    
	x_grid_dens_aux = grid_dens[, 1] - galaxy1_center[1]
	y_grid_dens_aux = grid_dens[, 2] - galaxy1_center[2]
	x_grid_dens = x_grid_dens_aux*cos(pi*(galaxy1_angle - 90)/180) - y_grid_dens_aux*sin(pi*(galaxy1_angle - 90)/180)
	y_grid_dens = x_grid_dens_aux*sin(pi*(galaxy1_angle - 90)/180) + y_grid_dens_aux*cos(pi*(galaxy1_angle - 90)/180)
	radius_grid_dens = sqrt(x_grid_dens^2 + y_grid_dens^2)
	theta_grid_dens_aux = atan(abs(y_grid_dens/x_grid_dens))
	theta_grid_dens = vector("numeric", length = length(theta_grid_dens_aux))
	for (i in 1:length(theta_grid_dens_aux)) {
		if (x_grid_dens[i] > 0 & y_grid_dens[i] > 0) {
			theta_grid_dens[i] = pi + theta_grid_dens_aux[i]
			} else if (x_grid_dens[i] < 0 & y_grid_dens[i] > 0) {
				theta_grid_dens[i] = 2*pi - theta_grid_dens_aux[i]
				} else if (x_grid_dens[i] < 0 & y_grid_dens[i] < 0) { 	 
					theta_grid_dens[i] = theta_grid_dens_aux[i]
					} else {
						theta_grid_dens[i] = pi - theta_grid_dens_aux[i]
						}	  		
		}	
	theta_grid_dens = 180*theta_grid_dens/pi	
	idx_residuals_matrix_neg_within = list(); idx_residuals_matrix_pos_within = list()
	idx_residuals_matrix_neg_within[[f]] = list()
	idx_residuals_matrix_pos_within[[f]] = list()
	for (l in 1:4) {				  			
		idx_residuals_matrix_neg_within[[f]][[l]] = intersect(idx_residuals_matrix_neg[[f]][[l]], idx_within_HSTfield)
		idx_residuals_matrix_pos_within[[f]][[l]] = intersect(idx_residuals_matrix_pos[[f]][[l]], idx_within_HSTfield)
		}	
	cat("		So far so good! Polar coordinates of the pixel centers calculated!\n")		
	
	# Interpolating the external contours of the footprints of the observations
	# to evaluate the polar coordinates of the associated vertices. These coordinates
	# will be returned by the function.
	if (name_galaxy == "ngc4472") {	
		polygons_external$x = c(187.4793, 187.4536, 187.4027, 187.3762, 
								187.3298, 187.3684, 187.3784, 187.3657,
 								187.4204, 187.4357, 187.4262, 187.4797)
		polygons_external$y = c(8.016989, 7.966688, 7.987010, 7.956322, 
								7.992714, 8.034543, 8.026278, 8.078866,
								8.097434, 8.041769, 8.039177, 8.017249)
		}	
	if (name_galaxy == c("ngc4472")) {	
		radiustheta_xy_extrapolate_polygons_external = radec_to_rtheta_boundary(polygons_external$x, 
																	   		    polygons_external$y, 
																	   		    galaxy1_center, 
																	   			galaxy1_angle)
		} else {
			radiustheta_xy_extrapolate_polygons_external = radec_to_rtheta(polygons_external$x, 
																		   polygons_external$y, 
																		   galaxy1_center, 
																		   galaxy1_angle)		
			}																		   			
	radius_extrapolate_polygons_external = radiustheta_xy_extrapolate_polygons_external[[1]][, 1]
	theta_extrapolate_polygons_external = radiustheta_xy_extrapolate_polygons_external[[1]][, 2]
	extrapolate_x_polygons_external = radiustheta_xy_extrapolate_polygons_external[[2]][, 1]
	extrapolate_y_polygons_external = radiustheta_xy_extrapolate_polygons_external[[2]][, 2]
	cat("		So far so good! Polar coordinates of the external footprint boundary calculated!\n")
	
	# Searching for structures.
	period = 36
	iaround_1sigmas = list()
	iaround_structures_1sigmas = list()
	size_iaround_structures_1sigmas = vector("numeric", length = length(idx_residuals_matrix_pos_within[[f]][[2]]))
	# Using pixels >= 1sigma as seeds, if any.
	if (length(idx_residuals_matrix_pos_within[[f]][[2]]) > 0) { 	
		for (i in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			icenter = idx_residuals_matrix_pos_within[[f]][[2]][i]
			iup = icenter - period
			idown = icenter + period
			ileft = icenter + 1
			iright = icenter - 1
			iupleft = icenter - period + 1
			iupright = icenter - period - 1
			idownleft = icenter + period + 1
			idownright = icenter + period - 1			  
			iaround_1sigmas[[i]] = c(iup, idown, ileft, iright, iupleft, iupright, idownleft, idownright)    
			iaround_1sigmas[[i]] = iaround_1sigmas[[i]][iaround_1sigmas[[i]] >= 0]  
			iaround_1sigmas[[i]] = intersect(iaround_1sigmas[[i]], idx_within_HSTfield)
			iaround_structures_1sigmas[[i]] = c(icenter, 
												iaround_1sigmas[[i]][which(num_sigmas[[f]][iaround_1sigmas[[i]]] >= 1)])
		    iaround_structures_1sigmas[[i]] = intersect(iaround_structures_1sigmas[[i]], idx_within_HSTfield)
		    }			
		}
	# Searching structures that overlap, and merging them.
	merged1_iaround_structures_1sigmas = list()
	merged2_iaround_structures_1sigmas = list()
	merged3_iaround_structures_1sigmas = list()
	merged4_iaround_structures_1sigmas = list()	
	merged5_iaround_structures_1sigmas = list()
	merged6_iaround_structures_1sigmas = list()	
	merged7_iaround_structures_1sigmas = list()
	merged8_iaround_structures_1sigmas = list()	
	merged9_iaround_structures_1sigmas = list()	
	merged10_iaround_structures_1sigmas = list()		
	merged11_iaround_structures_1sigmas = list()
	merged12_iaround_structures_1sigmas = list()				
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged1_iaround_structures_1sigmas[[j]] = iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged1_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged1_iaround_structures_1sigmas[[j]] = union(merged1_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}	
		}
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged2_iaround_structures_1sigmas[[j]] = merged1_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged2_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged2_iaround_structures_1sigmas[[j]] = union(merged2_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged3_iaround_structures_1sigmas[[j]] = merged2_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged3_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged3_iaround_structures_1sigmas[[j]] = union(merged3_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged4_iaround_structures_1sigmas[[j]] = merged3_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged4_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged4_iaround_structures_1sigmas[[j]] = union(merged4_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged5_iaround_structures_1sigmas[[j]] = merged4_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged5_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged5_iaround_structures_1sigmas[[j]] = union(merged5_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged6_iaround_structures_1sigmas[[j]] = merged5_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged6_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged6_iaround_structures_1sigmas[[j]] = union(merged6_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged7_iaround_structures_1sigmas[[j]] = merged6_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged7_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged7_iaround_structures_1sigmas[[j]] = union(merged7_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged8_iaround_structures_1sigmas[[j]] = merged7_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged8_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged8_iaround_structures_1sigmas[[j]] = union(merged8_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}														
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged9_iaround_structures_1sigmas[[j]] = merged8_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged9_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged9_iaround_structures_1sigmas[[j]] = union(merged9_iaround_structures_1sigmas[[j]], 	
															    iaround_structures_1sigmas[[k]])										   
				}
			}
		}		
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged10_iaround_structures_1sigmas[[j]] = merged9_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged10_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged10_iaround_structures_1sigmas[[j]] = union(merged10_iaround_structures_1sigmas[[j]], 	
															     iaround_structures_1sigmas[[k]])										   
				}
			}
		}	
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged11_iaround_structures_1sigmas[[j]] = merged10_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged11_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged11_iaround_structures_1sigmas[[j]] = union(merged11_iaround_structures_1sigmas[[j]], 	
															     iaround_structures_1sigmas[[k]])										   
				}
			}
		}			
	for (j in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
		merged12_iaround_structures_1sigmas[[j]] = merged11_iaround_structures_1sigmas[[j]]
		for (k in 1:length(idx_residuals_matrix_pos_within[[f]][[2]])) {
			intersect_structures = intersect(merged12_iaround_structures_1sigmas[[j]], iaround_structures_1sigmas[[k]])
			if (length(intersect_structures) > 0) {
				merged12_iaround_structures_1sigmas[[j]] = union(merged12_iaround_structures_1sigmas[[j]], 	
															     iaround_structures_1sigmas[[k]])										   
				}
			}
		}						
	sorted_merged_iaround_structures_1sigmas = lapply(merged12_iaround_structures_1sigmas, sort)	
	unique_merged_iaround_structures_1sigmas = unique(sorted_merged_iaround_structures_1sigmas)
	colors_merged_structures = sample(rainbow(length(unique_merged_iaround_structures_1sigmas)))
	cat("		So far so good! Merged structures in the residual maps determined!\n")
	
	# Determining the properties of the merged structures in the residual maps.
	average_significance_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	average_distance_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	average_theta_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	size_pixels_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	number_gc_excess_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	number_gc_total_unique_merged_structures = vector("numeric", length = length(unique_merged_iaround_structures_1sigmas))
	namegalaxy_unique_merged_structures = rep(name_galaxy, length(unique_merged_iaround_structures_1sigmas))
	thetas_unique_merged_structures = list(); range_thetas_unique_merged_structures = list()
	for (z in 1:length(unique_merged_iaround_structures_1sigmas)) {
		average_significance_unique_merged_structures[z] = mean(num_sigmas[[f]][unique_merged_iaround_structures_1sigmas[[z]]])
		average_distance_unique_merged_structures[z] = mean(dist_grid_dens[unique_merged_iaround_structures_1sigmas[[z]]])
		thetas_unique_merged_structures[[z]] = theta_grid_dens[unique_merged_iaround_structures_1sigmas[[z]]]
		range_thetas_unique_merged_structures[[z]] = range(thetas_unique_merged_structures[[z]])
		if ((range_thetas_unique_merged_structures[[z]][2] - range_thetas_unique_merged_structures[[z]][1]) > 180)	{
			thetas_unique_merged_structures[[z]][which(thetas_unique_merged_structures[[z]] > 180)] = 
				thetas_unique_merged_structures[[z]][which(thetas_unique_merged_structures[[z]] > 180)] - 360
			}		
		average_theta_unique_merged_structures[z] = mean(thetas_unique_merged_structures[[z]])
		if (average_theta_unique_merged_structures[z] < 0) {
			average_theta_unique_merged_structures[z] = average_theta_unique_merged_structures[z] + 360
			}
		size_pixels_unique_merged_structures[z] = length(unique_merged_iaround_structures_1sigmas[[z]])
		number_gc_excess_unique_merged_structures[z] = sum(numgcs_all[unique_merged_iaround_structures_1sigmas[[z]]]) - 
													   sum(sim_numgcs_all[unique_merged_iaround_structures_1sigmas[[z]]])
		number_gc_total_unique_merged_structures[z] = sum(numgcs_all[unique_merged_iaround_structures_1sigmas[[z]]])											   
		}												   
	cat("		So far so good! Global properties of the residual merged structures determined!\n")
	
	# Selection of the large extended residual structures. 
	# "Large", in this case, is defined as "containing more than 20 pixels".
	# "Intermediate" is defined as "containing more than 5 pixels but less than 20 pixels"
	idx_large_structures = which(size_pixels_unique_merged_structures >= 20) 
	idx_intermediate_structures = which(size_pixels_unique_merged_structures >= 5 & size_pixels_unique_merged_structures < 20) 
	idx_pixels_large_structures = list(); idx_pixels_intermediate_structures = list()
	size_pixels_large_structures = list(); size_pixels_intermediate_structures = list()
	for (j in 1:length(idx_large_structures)) {
		idx_pixels_large_structures[[j]] = unique_merged_iaround_structures_1sigmas[[idx_large_structures[j]]]
		size_pixels_large_structures[[j]] = length(idx_pixels_large_structures[[j]])
		}
	if (length(idx_intermediate_structures) > 0) {	
		for (h in 1:length(idx_intermediate_structures)) {
			idx_pixels_intermediate_structures[[h]] = unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[h]]]
			size_pixels_intermediate_structures[[h]] = length(idx_pixels_intermediate_structures[[h]])
			}	
		}			
	cat("		So far so good! Indices of the pixels belonging to the large and intermediate residual structures determined!\n")
	
	## Diagnostic plots.		
	# Plot of all the footprint of the observations only.
	pathplotfoot = paste(path_plot, name_galaxy, "/GCs/map_merged_footprint.pdf", sep = "")	
	pdf(pathplotfoot, width = 7, height = 7, paper = "special")
	#quartz()
	plot(grid_dens[, 1], grid_dens[, 2], cex = 0.2, col = "gray", 
		 pch = 20, 
		 xlim = c(max(grid_dens[, 1]), min(grid_dens[, 1])), 
		 axes = F, 
		 xlab = "", ylab = "")		
	points(galaxy1_center[1], galaxy1_center[2], cex = 2, pch = 16, col = "darkgreen")	
	color_palette_theta = colorRampPalette(c("red", "darkblue"))
	# Plotting angular distance from the major axis.
	points(grid_dens[, 1], 
		   grid_dens[, 2],
		   col = color_palette_theta(720)[cut(theta_grid_dens, seq(0, 360, length.out = 720))], 
		   pch = 7, 
		   cex = 1.5)	
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1)		
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle)			  
	if (name_galaxy == "ngc4649") {
		# Shape of the second galaxy.
		plot_galaxies(galaxy2_center, galaxy2_isodiameter, galaxy2_isodiameter_ratio, 
					  galaxy2_angle)
		}		
	# Plotting the interpolated points of the external footprint.
	points(extrapolate_x_polygons_external, 
		   extrapolate_y_polygons_external, 
		   pch = 20, 
		   cex = 3*radius_extrapolate_polygons_external/max(radius_extrapolate_polygons_external), 
		   col = color_palette_theta(720)[cut(theta_extrapolate_polygons_external, seq(0, 360, length.out = 720))])
	dev.off()	
	cat("		So far so good! Plot of the footprint of the HST data produced!\n")
	
	# Plot of all the merged structures from the residual map.
	pathplot = paste(path_plot, name_galaxy, "/GCs/map_merged_structures_", 
					 class_gc, "_K_", K_value, ".pdf", sep = "")	
	pdf(pathplot, width = 7, height = 7, paper = "special")
	#quartz()
	plot(grid_dens[, 1], grid_dens[, 2], cex = 0.2, col = "gray", 
		 pch = 20, 
		 xlim = c(max(grid_dens[, 1]), min(grid_dens[, 1])), 
		 axes = F, 
		 xlab = "", ylab = "")		
	points(galaxy1_center[1], galaxy1_center[2], cex = 2, pch = 16, col = "darkgreen")	
	color_palette_theta = colorRampPalette(c("red", "darkblue"))
	# Plotting angular distance from the major axis.
	points(grid_dens[, 1], 
		   grid_dens[, 2],
		   col = color_palette_theta(720)[cut(theta_grid_dens, seq(0, 360, length.out = 720))], 
		   pch = 7, 
		   cex = 1.5)	
	# Plotting the merged structures with different colors.	
	for (l in 1:length(unique_merged_iaround_structures_1sigmas)) {
		points(grid_dens[unique_merged_iaround_structures_1sigmas[[l]], 1], 
			   grid_dens[unique_merged_iaround_structures_1sigmas[[l]], 2],
		   	   col = colors_merged_structures[l], 
		       pch = 20, 
		       cex = cex_2ddens_maps)
		}	
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1)		
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle)			  
	if (name_galaxy == "ngc4649") {
		# Shape of the second galaxy.
		plot_galaxies(galaxy2_center, galaxy2_isodiameter, galaxy2_isodiameter_ratio, 
					  galaxy2_angle)
		}		
	# Plotting the interpolated points of the external footprint.
	points(extrapolate_x_polygons_external, 
		   extrapolate_y_polygons_external, 
		   pch = 20, 
		   cex = 3*radius_extrapolate_polygons_external/max(radius_extrapolate_polygons_external), 
		   col = color_palette_theta(720)[cut(theta_extrapolate_polygons_external, seq(0, 360, length.out = 720))])
	dev.off()	
	cat("		So far so good! Plot of all the merged structures produced!\n")
	
	# Producing the diagnostic plots of the position of the large and intermediate residual
	# structures in the residual plots and some variable that will be used later.
	idx_grid_dens_large_structures = list(); idx_grid_dens_intermediate_structures = list()
	chull_large_structure = list(); chull_large_structure_xy = list()
	chull_intermediate_structure = list(); chull_intermediate_structure_xy = list()
	x_squares_large_structure = list(); y_squares_large_structure = list()
	x_squares_intermediate_structure = list(); y_squares_intermediate_structure = list()
	xy_coords_large_contour = list(); xy_coords_intermediate_contour = list()
	density_large_structure = list(); density_intermediate_structure = list()
	number_gc_excess_intermediate_structure = list()
	number_gc_total_intermediate_structure = list()	
	sd_number_gc_excess_intermediate_structure = list()
	number_gc_excess_large_structure = list()
	number_gc_total_large_structure = list()
	sd_number_gc_excess_large_structure = list()	
	halfwidth_grid_kde2d = 3
	for (z in 1:length(idx_pixels_large_structures)) {
		idx_grid_dens_large_structures[[z]] = unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]]
		number_gc_total_large_structure[[z]] = sum(numgcs_all[idx_grid_dens_large_structures[[z]]])
		number_gc_excess_large_structure[[z]] = sum(numgcs_all[idx_grid_dens_large_structures[[z]]]) - 
												sum(sim_numgcs_all[idx_grid_dens_large_structures[[z]]])
		sd_number_gc_excess_large_structure[[z]] = sqrt(number_gc_excess_large_structure[[z]])									
		# Plot of the largest structures found in the residual maps (separately).
		pathplotlarge = paste(path_plot, name_galaxy, "/GCs/map_merged_structures_", 
				 		 	  class_gc, "_K_", K_value, "_largestructure_", z, ".pdf", sep = "")	
		pdf(pathplotlarge, width = 7, height = 7, paper = "special")
		plot(grid_dens[, 1], grid_dens[, 2], cex = 0.2, col = "gray", 
			 pch = 20, 
			 xlim = c(max(grid_dens[, 1]), min(grid_dens[, 1])), 
			 axes = F, 
			 xlab = "", ylab = "")		
		points(galaxy1_center[1], galaxy1_center[2], cex = 2, pch = 16, col = "darkgreen")	
		color_palette_theta = colorRampPalette(c("red", "darkblue"))
		# Plotting angular distance from the major axis.
		points(grid_dens[, 1], 
			   grid_dens[, 2],
			   col = color_palette_theta(720)[cut(theta_grid_dens, seq(0, 360, length.out = 720))], 
			   pch = 7, 
			   cex = 1.5)	
		## Plotting the pixels within the footprint of the
		## HST observations used.
		#points(grid_dens[idx_within_HSTfield, 1], 
		#	   grid_dens[idx_within_HSTfield, 2],
		#	   col = "black",
		#	   pch = 15, 
		#	   cex = 1.2)		
		# Plotting the large residual structure.	
		points(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]], 1], 
			   grid_dens[unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]], 2],
			   col = colors_merged_structures[idx_large_structures[z]], 
			   pch = 20, 
			   cex = cex_2ddens_maps)
		# Plotting a round curve around the pixels occupied by the 
		# large structure. 
		halfstep_x = 0.002; halfstep_y = 0.002
		x_squares_large_structure_pixel = list()
		y_squares_large_structure_pixel = list()
		for (k in 1:length(unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]])) {
			x_squares_large_structure_pixel[[k]] = rep(rep(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]][k], 1], 2) + 					
												   	   c(halfstep_x, -halfstep_x), 2)	# Rdius*sin(angles_radians)
			y_squares_large_structure_pixel[[k]] = rep(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_large_structures[z]]][k], 2], 2) + 					
												   c(halfstep_y, -halfstep_y)		#radius*cos(angles_radians)									  
			y_squares_large_structure_pixel[[k]] = c(y_squares_large_structure_pixel[[k]], rev(y_squares_large_structure_pixel[[k]]))
			}
		x_squares_large_structure[[z]] = unlist(x_squares_large_structure_pixel)
		y_squares_large_structure[[z]] = unlist(y_squares_large_structure_pixel)
		density_large_structure[[z]] = kde2d(x_squares_large_structure[[z]], 
											 y_squares_large_structure[[z]], 
											 lims = c(min(x_squares_large_structure[[z]]) - halfwidth_grid_kde2d*halfstep_x, 
											 		  max(x_squares_large_structure[[z]]) + halfwidth_grid_kde2d*halfstep_x, 
											 		  min(y_squares_large_structure[[z]]) - halfwidth_grid_kde2d*halfstep_y, 
											 		  max(y_squares_large_structure[[z]]) + halfwidth_grid_kde2d*halfstep_y))
		contour(density_large_structure[[z]], 
				col = "darkgreen", 
				lwd = 2, 
				nlevels = 10, 
				labels = "", 
				add = T)
		contour(density_large_structure[[z]], 
				col = "red", 
				lwd = 3, 
				levels = max(density_large_structure[[z]]$z)/factor_galaxy[[name_galaxy]], 
				labels = "", 
				add = T)
		x_grid = seq(min(x_squares_large_structure[[z]]) - halfwidth_grid_kde2d*halfstep_x, 
					 max(x_squares_large_structure[[z]]) + halfwidth_grid_kde2d*halfstep_x, 
					 length.out = 25)
		y_grid = seq(min(y_squares_large_structure[[z]]) - halfwidth_grid_kde2d*halfstep_y, 
					 max(y_squares_large_structure[[z]]) + halfwidth_grid_kde2d*halfstep_y, 
					 length.out = 25)
		xy_coords_large_contour[[z]] = contourLines(x_grid, 
							 			 	  	    y_grid, 
							 			 	  		density_large_structure[[z]]$z,
							 			 	  		levels = max(density_large_structure[[z]]$z)/factor_galaxy[[name_galaxy]])			 					 
		lines(xy_coords_large_contour[[z]][[1]]$x, 
			  xy_coords_large_contour[[z]][[1]]$y, 
			  lwd = 2, 
			  lty = 2,
			  col = "black")					  
		# Drawing the polygons of the HST fields. 
		plot_footprint(polygons, polygons_external, polygons_external_avoid, 1)		
		# Shape of the first galaxy.
		plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
					  galaxy1_angle)			  
		if (name_galaxy == "ngc4649") {
			# Shape of the second galaxy.
			plot_galaxies(galaxy2_center, galaxy2_isodiameter, galaxy2_isodiameter_ratio, 
						  galaxy2_angle)
			}		
		# Plotting the interpolated points of the external footprint.
		points(extrapolate_x_polygons_external, 
		  	   extrapolate_y_polygons_external, 
		   	   pch = 20, 
		   	   cex = 3*radius_extrapolate_polygons_external/max(radius_extrapolate_polygons_external), 
		   	   col = color_palette_theta(720)[cut(theta_extrapolate_polygons_external, seq(0, 360, length.out = 720))])
		dev.off()	
		}
	if (length(idx_pixels_intermediate_structures) > 0) {	
		for (z in 1:length(idx_pixels_intermediate_structures)) {
			idx_grid_dens_intermediate_structures[[z]] = unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]]
			number_gc_total_intermediate_structure[[z]] = sum(numgcs_all[idx_grid_dens_intermediate_structures[[z]]])
			number_gc_excess_intermediate_structure[[z]] = sum(numgcs_all[idx_grid_dens_intermediate_structures[[z]]]) - 
												    	   sum(sim_numgcs_all[idx_grid_dens_intermediate_structures[[z]]])
			sd_number_gc_excess_intermediate_structure[[z]] = sqrt(number_gc_excess_intermediate_structure[[z]])
			# Plot of the largest structures found in the residual maps (separately).
			pathplotlarge = paste(path_plot, name_galaxy, "/GCs/map_merged_structures_", 
					 		 	  class_gc, "_K_", K_value, "_intermediatestructure_", z, ".pdf", sep = "")	
			pdf(pathplotlarge, width = 7, height = 7, paper = "special")
			plot(grid_dens[, 1], grid_dens[, 2], cex = 0.2, col = "gray", 
				 pch = 20, 
				 xlim = c(max(grid_dens[, 1]), min(grid_dens[, 1])), 
				 axes = F, 
				 xlab = "", ylab = "")		
			points(galaxy1_center[1], galaxy1_center[2], cex = 2, pch = 16, col = "darkgreen")	
			color_palette_theta = colorRampPalette(c("red", "darkblue"))
			# Plotting angular distance from the major axis.
			points(grid_dens[, 1], 
				   grid_dens[, 2],
				   col = color_palette_theta(720)[cut(theta_grid_dens, seq(0, 360, length.out = 720))], 
				   pch = 7, 
				   cex = 1.5)	
			# Plotting the intermediate residual structure.	
			points(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]], 1], 
				   grid_dens[unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]], 2],
				   col = colors_merged_structures[idx_intermediate_structures[z]], 
				   pch = 20, 
				   cex = cex_2ddens_maps)
			# Plotting a round curve around the pixels occupied by the 
			# intermediate structure. 
			halfstep_x = 0.002; halfstep_y = 0.002
			x_squares_intermediate_structure_pixel = list()
			y_squares_intermediate_structure_pixel = list()
			for (k in 1:length(unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]])) {
				x_squares_intermediate_structure_pixel[[k]] = rep(rep(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]][k], 1], 2) + 					
													   	   		  c(halfstep_x, -halfstep_x), 2)	
				y_squares_intermediate_structure_pixel[[k]] = rep(grid_dens[unique_merged_iaround_structures_1sigmas[[idx_intermediate_structures[z]]][k], 2], 2) + 					
													   			  c(halfstep_y, -halfstep_y)											  
				y_squares_intermediate_structure_pixel[[k]] = c(y_squares_intermediate_structure_pixel[[k]], rev(y_squares_intermediate_structure_pixel[[k]]))
				}
			x_squares_intermediate_structure[[z]] = unlist(x_squares_intermediate_structure_pixel)
			y_squares_intermediate_structure[[z]] = unlist(y_squares_intermediate_structure_pixel)
			density_intermediate_structure[[z]] = kde2d(x_squares_intermediate_structure[[z]], 
												 		y_squares_intermediate_structure[[z]], 
												 		lims = c(min(x_squares_intermediate_structure[[z]]) - halfwidth_grid_kde2d*halfstep_x, 
												 		  		 max(x_squares_intermediate_structure[[z]]) + halfwidth_grid_kde2d*halfstep_x, 
												 		  		 min(y_squares_intermediate_structure[[z]]) - halfwidth_grid_kde2d*halfstep_y, 
												 		  		 max(y_squares_intermediate_structure[[z]]) + halfwidth_grid_kde2d*halfstep_y))
			contour(density_intermediate_structure[[z]], 
					col = "darkgreen", 
					lwd = 2, 
					nlevels = 10, 
					labels = "", 
					add = T)
			contour(density_intermediate_structure[[z]], 
					col = "red", 
					lwd = 3, 
					levels = max(density_intermediate_structure[[z]]$z)/factor_galaxy[[name_galaxy]], 
					labels = "", 
					add = T)
			x_grid = seq(min(x_squares_intermediate_structure[[z]]) - halfwidth_grid_kde2d*halfstep_x, 
						 max(x_squares_intermediate_structure[[z]]) + halfwidth_grid_kde2d*halfstep_x, 
						 length.out = 25)
			y_grid = seq(min(y_squares_intermediate_structure[[z]]) - halfwidth_grid_kde2d*halfstep_y, 
						 max(y_squares_intermediate_structure[[z]]) + halfwidth_grid_kde2d*halfstep_y, 
						 length.out = 25)
			xy_coords_intermediate_contour[[z]] = contourLines(x_grid, 
								 			 	  	    	   y_grid, 
								 			 	  			   density_intermediate_structure[[z]]$z,
								 			 	  			   levels = max(density_intermediate_structure[[z]]$z)/factor_galaxy[[name_galaxy]])			 					 
			lines(xy_coords_intermediate_contour[[z]][[1]]$x, 
				  xy_coords_intermediate_contour[[z]][[1]]$y, 
				  lwd = 2, 
				  lty = 2,
				  col = "black")					  
			# Drawing the polygons of the HST fields. 
			plot_footprint(polygons, polygons_external, polygons_external_avoid, 1)		
			# Shape of the first galaxy.
			plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
						  galaxy1_angle)			  
			if (name_galaxy == "ngc4649") {
				# Shape of the second galaxy.
				plot_galaxies(galaxy2_center, galaxy2_isodiameter, galaxy2_isodiameter_ratio, 
							  galaxy2_angle)
				}		
			# Plotting the interpolated points of the external footprint.
			points(extrapolate_x_polygons_external, 
			  	   extrapolate_y_polygons_external, 
			   	   pch = 20, 
			   	   cex = 3*radius_extrapolate_polygons_external/max(radius_extrapolate_polygons_external), 
			   	   col = color_palette_theta(720)[cut(theta_extrapolate_polygons_external, seq(0, 360, length.out = 720))])
			dev.off()	
			}	
		}
		
	# Plot of the boundaries of the HST observations in the radial distance vs angle from 
	# major axis plane. 
	pathplot2 = paste(path_plot, name_galaxy, "/GCs/footprint_boundaries_r_vs_angle.pdf", sep = "")	
	pdf(pathplot2, width = 7, height = 7, paper = "special")
	#quartz()
	plot(radius_extrapolate_polygons_external, theta_extrapolate_polygons_external, 
		 pch = 20, 
		 col = color_palette_theta(720)[cut(theta_extrapolate_polygons_external, 
		 									seq(0, 360, length.out = 720))], 
		 cex = 3*radius_extrapolate_polygons_external/max(radius_extrapolate_polygons_external), 
		 xlab = "Galacto-centric distance [deg]", 
		 ylab = "Angle from major axis [deg]")
	dev.off()		
	cat("		So far so good! Plot of the boundary of the HST footprint produced!\n")
			
	# Returning the values of the parameters of the merged residual structures.
	if (mode == "all") {
		cat(paste("	Fraction of pixels in the residual structure for galaxy ", name_galaxy, "!\n", sep = ""))
		cat(paste("	Number of pixels in the residual structures: ", length(unlist(unique_merged_iaround_structures_1sigmas))), "!\n", sep = "")
		cat(paste("	Number of pixels within the footprint of the HST data: ", length(idx_within_HSTfield)), "!\n", sep = "")
		cat(paste("	Fraction of pixels within the residual structures: ", round(100*length(unlist(unique_merged_iaround_structures_1sigmas))/length(idx_within_HSTfield), 2)), "!\n", 
					sep = "")
		cat(paste("	Number of pixels in the large spatially coherent residual structures: ", length(unlist(idx_grid_dens_large_structures))), "!\n", sep = "")				
		cat(paste("	Fraction of pixels within the large spatially coherent residual structures: ", round(100*length(unlist(idx_grid_dens_large_structures))/length(idx_within_HSTfield), 2)), "!\n", 
				  sep = "")
		data_structures = cbind(size_pixels_unique_merged_structures, average_significance_unique_merged_structures, 
					 			average_distance_unique_merged_structures, average_theta_unique_merged_structures, 
					 			namegalaxy_unique_merged_structures, number_gc_excess_unique_merged_structures, 
					 			number_gc_total_unique_merged_structures) 
		data_boundaries = cbind(radius_extrapolate_polygons_external, theta_extrapolate_polygons_external)
		data_structures_boundaries = list(structures = data_structures, boundaries = data_boundaries, 
										  galaxy_center = galaxy1_center, galaxy_angle = galaxy1_angle, 
										  idx_large_structures = idx_grid_dens_large_structures, 
										  idx_largeintermediate_structures = c(idx_large_structures, idx_intermediate_structures))
		return(data_structures_boundaries)
		cat("		So far so good! Properties of the merged structures returned!\n")
		} else if (mode == "large") {
			return(xy_coords_large_contour)
			} else if (mode == "large_intermediate") {
				xy_coords_largeintermediate_contour = list()
				number_gc_excess_largeintermediate_structure = list()
				number_gc_total_largeintermediate_structure = list()
				idx_largeintermediate_structure = list()
				size_pixels_largeintermediate_structure = list()
					for (i in 1:length(xy_coords_large_contour)) {
					xy_coords_largeintermediate_contour[[i]] = xy_coords_large_contour[[i]]
					number_gc_excess_largeintermediate_structure[[i]] = number_gc_excess_large_structure[[i]]
					number_gc_total_largeintermediate_structure[[i]] = number_gc_total_large_structure[[i]]
					idx_largeintermediate_structure[[i]] = idx_large_structures[i]
					size_pixels_largeintermediate_structure[[i]] = size_pixels_large_structures[[i]]
						}
				if (length(xy_coords_intermediate_contour) > 0) {
					for (k in 1:length(xy_coords_intermediate_contour)) {
						xy_coords_largeintermediate_contour[[k + length(xy_coords_large_contour)]] = xy_coords_intermediate_contour[[k]]
						number_gc_excess_largeintermediate_structure[[k + length(number_gc_excess_large_structure)]] = number_gc_excess_intermediate_structure[[k]]
						number_gc_total_largeintermediate_structure[[k + length(number_gc_total_large_structure)]] = number_gc_total_intermediate_structure[[k]]
						idx_largeintermediate_structure[[k + length(number_gc_excess_large_structure)]] = idx_intermediate_structures[k]
						size_pixels_largeintermediate_structure[[k + length(xy_coords_large_contour)]] = size_pixels_intermediate_structures[[k]]
						}	
					}
				variables_largeintermediate = list(xy_coords_contour = xy_coords_largeintermediate_contour, 
												   number_gc_excess = number_gc_excess_largeintermediate_structure, 
												   idx_structure = idx_largeintermediate_structure, 
												   number_gc_total = number_gc_total_largeintermediate_structure, 
												   size_structure = size_pixels_largeintermediate_structure)	
				return(variables_largeintermediate)
				} else if (mode == "large_intermediate_mrd") {
					xy_coords_largeintermediate_contour = list()
					number_gc_excess_largeintermediate_structure = list()
					sd_number_gc_excess_largeintermediate_structure = list()
					idx_pixels_largeintermediate_structure = list()
					for (i in 1:length(xy_coords_large_contour)) {
						xy_coords_largeintermediate_contour[[i]] = xy_coords_large_contour[[i]]
						number_gc_excess_largeintermediate_structure[[i]] = number_gc_excess_large_structure[[i]]
						sd_number_gc_excess_largeintermediate_structure[[i]] = sd_number_gc_excess_large_structure[[i]]
						idx_pixels_largeintermediate_structure[[i]] = idx_pixels_large_structures[[i]]
						}
					if (length(xy_coords_intermediate_contour) > 0) {
						for (k in 1:length(xy_coords_intermediate_contour)) {
							xy_coords_largeintermediate_contour[[k + length(xy_coords_large_contour)]] = xy_coords_intermediate_contour[[k]]
							number_gc_excess_largeintermediate_structure[[k + length(number_gc_excess_large_structure)]] = number_gc_excess_intermediate_structure[[k]]
							sd_number_gc_excess_largeintermediate_structure[[k + length(number_gc_excess_large_structure)]] = sd_number_gc_excess_intermediate_structure[[k]]
							idx_pixels_largeintermediate_structure[[k + length(xy_coords_large_contour)]] = idx_pixels_intermediate_structures[[k]]
							}	
						}
					variables_largeintermediate = list(xy_coords_contour = xy_coords_largeintermediate_contour, 
													   number_gc_excess = number_gc_excess_largeintermediate_structure, 
													   sd_number_gc_excess = sd_number_gc_excess_largeintermediate_structure,
													   pixels_structures = idx_pixels_largeintermediate_structure)	
					return(variables_largeintermediate)
					}
	}