## This function returns a description of the GCs structures located in the residual maps 
## generated with the minimal radial density profiles.
##
## Arguments:
##		name_galaxy 					<- 	Name of the galaxy.
##		K_value 						<- 	K value.
##		num_simulations_original 		<-	Number of simulations for the experiment 
##											performed with the natural radial density
##											profile.
##		num_simulations_mrd		 		<-	Number of simulations for the experiment 
##											performed with the rescaled radial density
##											profile (mrd).
##	    nsteps_smooth					<-	Number of steps along the two axis.
##		sim_type 						<-	Simulation type (always "observed")
## 		properties_structures_original	<- 	Structure containing the properties of the 
##											residual structures detected in the map
##											derived from the natural radial density profile.
## 
## Written by R. D'Abrusco
##
## Last modified: 15/8/14.

structures_radius_mrd <- function(name_galaxy, 
							      K_value = 9, 
							      num_simulations_original = 500, 
							      num_simulations_mrd = 500,  
							      nsteps_smooth = 36, 
							      sim_type = "observed", 
							      properties_largeintermediate_structures) {
	
	# Loading required libraries.
	source("./raycasting.R")							# Used to determine grid points within 	
														# the overall field region.
	source("./auxiliary.R")								# Loading the auxiliary functions.	
	suppressMessages(library(plotrix, quietly = TRUE))	# Used to plot other shapes and lines.	
	suppressMessages(library(MASS, quietly = TRUE))		# Used to compute 2D densities.				
	cat("		So far so good! Libraries loaded!\n")

	# Setting the paths to the location of the files (data and plots)
	# that will be useful later.
	name_galaxy_mrd = paste(name_galaxy, "_mrd", sep = "")
	class_gc = "all"
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
	name_file_session_mrd = paste(path_data, name_galaxy_mrd, "/Sessions/gc_nsims_", 
								  num_simulations_mrd, "_nsteps_", nsteps_smooth, 
								  "_simtype_", sim_type, ".Rdata", sep = "")
	load(name_file_session_mrd)
	cat("		So far so good! Image of a previous session imported!\n")
	
	# Evaluating the values of the variables necessary to plot the residual map.
	f = which(k_vector == K_value)
	residuals_matrix_all_mrd = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_all_aux_mrd = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_all_mrd = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_all_mrd = apply(res_all_mrd[[f]], c(2, 3), sum)/num_simulations
	errresiduals_matrix_all_aux_mrd = apply(errres_all_mrd[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_all_mrd = sqrt(errresiduals_matrix_all_aux_mrd)/num_simulations
	if (name_galaxy == "ngc4649") {
		num_sigmas_threshold = c(0.5, 1, 2, 3)
		idx_residuals_matrix_all_pos_mrd = list(); idx_residuals_matrix_all_neg_mrd = list()
		idx_residuals_matrix_all_neg_mrd[[f]] = list(); idx_residuals_matrix_all_pos_mrd[[f]] = list()
		for (z in seq(length(num_sigmas_threshold))) {
			idx_residuals_matrix_all_pos_mrd[[f]][[z]] = which(num_sigmas_all_mrd[[f]] >= num_sigmas_threshold[z])
			idx_residuals_matrix_all_neg_mrd[[f]][[z]] = which(num_sigmas_all_mrd[[f]] <= -num_sigmas_threshold[z])
			}	
		idx_within_HSTfield = idx_gridbins_galaxies
		}
	cat("		So far so good! Variables necessary to produce the plot read!\n")		
	
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	num_sigmas_mrd = num_sigmas_all_mrd
	idx_residuals_matrix_neg_mrd = idx_residuals_matrix_all_neg_mrd
	idx_residuals_matrix_pos_mrd = idx_residuals_matrix_all_pos_mrd
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
	idx_residuals_matrix_neg_within_mrd = list(); idx_residuals_matrix_pos_within_mrd = list()
	idx_residuals_matrix_neg_within_mrd[[f]] = list()
	idx_residuals_matrix_pos_within_mrd[[f]] = list()
	for (l in 1:4) {				  			
		idx_residuals_matrix_neg_within_mrd[[f]][[l]] = intersect(idx_residuals_matrix_neg_mrd[[f]][[l]], idx_within_HSTfield)
		idx_residuals_matrix_pos_within_mrd[[f]][[l]] = intersect(idx_residuals_matrix_pos_mrd[[f]][[l]], idx_within_HSTfield)
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
	
	# Loading the properties of the residual structures extracted from the residual maps
	# obtained from the observed radial density profile. 
	number_largeintermediate_structures = length(properties_largeintermediate_structures$pixels_structures)
	number_gc_excess_unique_merged_structures_mrd = vector("numeric", length = length(number_largeintermediate_structures))
	sd_number_gc_excess_unique_merged_structures_mrd = vector("numeric", length = length(number_largeintermediate_structures))
	for (z in seq(number_largeintermediate_structures)) {
		number_gc_excess_unique_merged_structures_mrd[z] = sum(numgcs_all_mrd[properties_largeintermediate_structures$pixels_structures[[z]]]) - 
												 	  	   sum(sim_numgcs_all_mrd[properties_largeintermediate_structures$pixels_structures[[z]]])
		sd_number_gc_excess_unique_merged_structures_mrd[z] = sqrt(number_gc_excess_unique_merged_structures_mrd[z])
		}
	colors_merged_structures = sample(rainbow(length(properties_largeintermediate_structures$pixels_structures)))
	
	## Diagnostic plots.		
	# Plot of all the footprint of the observations only.
	pathplotfoot = paste(path_plot, name_galaxy_mrd, "/GCs/map_merged_footprint.pdf", sep = "")	
	pdf(pathplotfoot, width = 7, height = 7, paper = "special")
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
	pathplot = paste(path_plot, name_galaxy_mrd, "/GCs/map_merged_structures_", 
					 class_gc, "_K_", K_value, ".pdf", sep = "")	
	pdf(pathplot, width = 7, height = 7, paper = "special")
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
	for (l in seq(number_largeintermediate_structures)) {
		points(grid_dens[properties_largeintermediate_structures$pixels_structures[[l]], 1], 
			   grid_dens[properties_largeintermediate_structures$pixels_structures[[l]], 2],
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
	number_gc_excess_large_structure = list()
	halfwidth_grid_kde2d = 3
	for (z in seq(number_largeintermediate_structures)) {
		idx_grid_dens_large_structures[[z]] = properties_largeintermediate_structures$pixels_structures[[z]]
		number_gc_excess_large_structure[[z]] = sum(numgcs_all_mrd[idx_grid_dens_large_structures[[z]]]) - 
												sum(sim_numgcs_all_mrd[idx_grid_dens_large_structures[[z]]])
		# Plot of the largest structures found in the residual maps (separately).
		pathplotlargeintermediate = paste(path_plot, name_galaxy_mrd, "/GCs/map_merged_structures_", 
				 		 	  			  class_gc, "_K_", K_value, "_largeintermediate_", z, ".pdf", sep = "")	
		pdf(pathplotlargeintermediate, width = 7, height = 7, paper = "special")
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
		# Plotting the large or intermediate residual structure.	
		points(grid_dens[properties_largeintermediate_structures$pixels_structures[[z]], 1], 
			   grid_dens[properties_largeintermediate_structures$pixels_structures[[z]], 2],
			   col = colors_merged_structures[z], 
			   pch = 20, 
			   cex = cex_2ddens_maps)
		# Plotting a round curve around the pixels occupied by the 
		# large structure. 
		halfstep_x = 0.002; halfstep_y = 0.002
		x_squares_large_structure_pixel = list()
		y_squares_large_structure_pixel = list()
		for (k in seq(properties_largeintermediate_structures$pixels_structures[[z]])) {
			x_squares_large_structure_pixel[[k]] = rep(rep(grid_dens[properties_largeintermediate_structures$pixels_structures[[z]][k], 1], 2) + 					
												   	   c(halfstep_x, -halfstep_x), 2)
			y_squares_large_structure_pixel[[k]] = rep(grid_dens[properties_largeintermediate_structures$pixels_structures[[z]][k], 2], 2) + 					
												   	   c(halfstep_y, -halfstep_y)										  
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

	# Returning the values of the parameters of the merged residual structures.
	variables_largeintermediate_mrd = list(number_gc_excess = number_gc_excess_unique_merged_structures_mrd, 
										   sd_number_gc_excess = sd_number_gc_excess_unique_merged_structures_mrd)	
	return(variables_largeintermediate_mrd)
	}