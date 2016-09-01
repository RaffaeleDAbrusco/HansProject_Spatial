## Collection of auxiliary functions performing operations required in the
## main program for the characterization of the spatial distribution of GCs
## and LMXBs.
##
## Written by R. D'Abrusco
##
## Last modified: 10/9/15.
#
#	10/9/15		->		Modified "plot_scale" to accept >= 3600 arcsec and 
#						use degrees as units. Modified "plot_labels" for 
#						an improved look of the label box.

addalpha <- function(colorname, alpha = 100) {
	newcolor <- col2rgb(colorname)
	apply(newcolor, 2, function(curcoldata){rgb(red = curcoldata[1], green = curcoldata[2],
    										    blue = curcoldata[3], alpha = alpha, maxColorValue = 255)})
	}

MHmakeRandomString <- function(n = 1, 
							   lenght = 12)  {
    randomString <- c(1:n)     # initialize vector
    for (i in 1:n)	{
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                 lenght, 
                                 replace = TRUE),
                                 collapse = "")
   		}
    return(randomString)
}

read_galaxy_metadata <- function(galaxy_name,
								 path_data, 
								 primary_galaxy) {
	# This function reads the essential information
	# about the galaxy whose GCs distribution will
	# be considered.
	path_file_galaxy = paste(path_data, primary_galaxy, "/", galaxy_name, "_general_data.dat", sep = "")
	primary_galaxy_data = read.table(path_file_galaxy)
	# Coordinates (in degrees) of the center of the primary galaxy.
	galaxy_center = c(primary_galaxy_data$V1, primary_galaxy_data$V2)
	# Mean log of the apparent major isophotal diameter (m_B = 25 B/mag^2),
	# measured in arcmin.
	galaxy_isodiameter= primary_galaxy_data$V3
	# Mean logarithm of the ratio of the D25 isophotal diameters.
	galaxy_isodiameter_ratio = primary_galaxy_data$V4
	# Position angle of the main axis of the galaxy [degrees].
	galaxy_angle = primary_galaxy_data$V5
	# Return everything.
	output = c(galaxy_center, galaxy_isodiameter, galaxy_isodiameter_ratio, 
			   galaxy_angle)
	}
	
read_field_metadata <- function(field_name,
								path_data) {
	# This function reads the essential information
	# about the field whose GCs distribution will
	# be investigated.
	path_file_field = paste(path_data, field_name, "/", field_name, "_general_data.dat", sep = "")
	primary_field_data = read.table(path_file_field)
	# Names of the galaxies in the field.
	galaxies_name = primary_field_data$V1
	# Coordinates (in degrees) of the center of the galaxies in the field.
	galaxies_center = cbind(primary_field_data$V2, primary_field_data$V3)
	# Mean log of the apparent major isophotal diameter (m_B = 25 B/mag^2),
	# measured in arcmin.
	galaxies_isodiameter = cbind(primary_field_data$V4)
	# Mean logarithm of the ratio of the D25 isophotal diameters.
	galaxies_isodiameter_ratio = cbind(primary_field_data$V5)
	# Position angle of the main axis of the galaxy [degrees].
	galaxies_angle = cbind(primary_field_data$V6)
	# RA and Dec. interval extremes.
	ra_range_min = primary_field_data$V7
	ra_range_max = primary_field_data$V8
	dec_range_min = primary_field_data$V9
	dec_range_max = primary_field_data$V10
	# Reading parameters of the binning. 
	max_radial_distance = primary_field_data$V11
	binwidth_radial_distance = primary_field_data$V12
	# Return all the parameters read.
	output = data.frame(as.character(galaxies_name), galaxies_center, galaxies_isodiameter, galaxies_isodiameter_ratio, 
			   	   	    galaxies_angle, ra_range_min, ra_range_max, dec_range_min, dec_range_max, 
			   	   	    max_radial_distance, binwidth_radial_distance)
	}	
	
read_galaxy_metadata_footprints <- function(metadata_file_name) {
	# This function reads the essential information
	# about the galaxy whose HST footprints are to be
	# characterized
	primary_galaxy_data = read.table(metadata_file_name)
	primary_galaxy_data = primary_galaxy_data[1, ]
	# Coordinates (in degrees) of the center of the primary galaxy.
	galaxy_center = c(primary_galaxy_data$V1, primary_galaxy_data$V2)
	# Mean log of the apparent major isophotal diameter (m_B = 25 B/mag^2),
	# measured in arcmin.
	galaxy_isodiameter= primary_galaxy_data$V3
	# Mean logarithm of the ratio of the D25 isophotal diameters.
	galaxy_isodiameter_ratio = primary_galaxy_data$V4
	# Position angle of the main axis of the galaxy [degrees].
	galaxy_angle = primary_galaxy_data$V5
	# Return everything.
	output = c(galaxy_center, galaxy_isodiameter, galaxy_isodiameter_ratio, 
			   galaxy_angle)
	}	
	
plot_galaxies <- function(galaxy_center, 
						  galaxy_isodiameter, 
						  galaxy_isodiameter_ratio, 
						  galaxy_angle, 
						  lwd_ellipse = 1.5, 
						  cex_center = 5, 
						  scale_galaxy = 1, 
						  color_galaxy_border = "brown", 
						  color_galaxy = NULL) {
	# This function plots the D25 shape of a galaxy, provided with the coordinates 
	# of the center of the galaxy, the D25 diameter of the galaxy, the ratio of the
	# major and minor axes and the angle of the major axis of the galaxy with RA axis.
	# Plot the center of the galaxy.
	points(galaxy_center[1], galaxy_center[2], 
	   	   pch = 3, 
	   	   cex = cex_center, 
	   	   col = "darkgreen")
	# Plot the D25 ellipse.  		
	draw.ellipse(galaxy_center[1], galaxy_center[2], 
				 a = scale_galaxy*(10^(galaxy_isodiameter - 1))/120, 
			 	 b = scale_galaxy*((10^(galaxy_isodiameter - 1))/120)/10^(galaxy_isodiameter_ratio), 
			 	 angle = 270 - galaxy_angle,
			 	 nv = 1000, 
			 	 border = color_galaxy_border,
			 	 col = color_galaxy,
			 	 lwd = lwd_ellipse, 
			 	 deg = TRUE)
	# Determining the coordinates of the points required
	# to draw the major axis of the galaxy.
	if (galaxy_angle == 0) {
		galaxy_angle = galaxy_angle - 1
		}
	slope = tan(pi/2 - galaxy_angle*pi/180)
	intercept = galaxy_center[2] - galaxy_center[1]*tan(pi/2 - galaxy_angle*pi/180)
	majoraxis = scale_galaxy*(10^(galaxy_isodiameter - 1))/120
	minoraxis = scale_galaxy*((10^(galaxy_isodiameter - 1))/120)/10^(galaxy_isodiameter_ratio)
	x_intersections = polyroot(c((galaxy_center[1]^2 - majoraxis^2 + (intercept - galaxy_center[2])^2),  
								 (2*slope*intercept - 2*galaxy_center[1] - 2*galaxy_center[2]*slope), 
								 (slope^2 + 1)))
	x_intersections = Mod(x_intersections)
	y_intersections = x_intersections*slope + intercept
	# Plot the major axis of the galaxy.
	segments(x_intersections[1], y_intersections[1], 
			 x_intersections[2], y_intersections[2], 
			 col = color_galaxy_border,
			 lwd = lwd_ellipse)	    
	}

return_extreme_axis <- function(galaxy_center, galaxy_isodiameter, galaxy_isodiameter_ratio, 
						  	 	galaxy_angle) {
	# This function returns the cartesian coordinates of the intersections of the major
	# axis of a galaxy with its D25 isophote. It is used to determine the angular distance
	# of a pixel center from the major axis using trigonometry.
	if (galaxy_angle == 0) {
		galaxy_angle = galaxy_angle - 1
		}
	slope = tan(pi/2 - galaxy_angle*pi/180)
	intercept = galaxy_center[2] - galaxy_center[1]*tan(pi/2 - galaxy_angle*pi/180)
	majoraxis = (10^(galaxy_isodiameter - 1))/120
	minoraxis = ((10^(galaxy_isodiameter - 1))/120)/10^(galaxy_isodiameter_ratio)
	x_intersections = polyroot(c((galaxy_center[1]^2 - majoraxis^2 + (intercept - galaxy_center[2])^2),  
								 (2*slope*intercept - 2*galaxy_center[1] - 2*galaxy_center[2]*slope), 
								 (slope^2 + 1)))
	x_intersections = Mod(x_intersections)
	y_intersections = x_intersections*slope + intercept
	return(cbind(x_intersections, y_intersections))
	}	
	
plot_footprint <- function(polygons, polygon_external, polygon_external_avoid,
						   only_avoid, 
						   color_ext = FALSE, 
						   alpha_ext = 255, 
						   color_border = "black",
						   lty_border = 1,
						   diff) {
	if (length(polygons) != 0) {					   
		## This function plots the polygons corresponding to each
		## field observed by HST, the overall shape of the field
		## and a few other graphical elements.
		num_polygons = length(polygons)
		# Splitting the variables to draw the single polygons.
		x_polya = list(); y_polya = list()
		x_polyb = list(); y_polyb = list()
		for (i in 1:num_polygons) {
			x_polya[[i]] = polygons[[i]]$a$x; y_polya[[i]] = polygons[[i]]$a$y 
			x_polyb[[i]] = polygons[[i]]$b$x; y_polyb[[i]] = polygons[[i]]$b$y 		
			if (only_avoid == 0) {
				polygon(x_polya[[i]], y_polya[[i]], border = "black", 
						lwd = 0.5, lty = 2, col = NA)
				polygon(x_polyb[[i]], y_polyb[[i]], border = "black", 
						lwd = 0.5, lty = 2, col = NA)				
				}
			}
		}		
	glob_poly_x = polygon_external$x; glob_poly_y = polygon_external$y 
	glob_poly_outside_x = polygon_external_avoid$x; glob_poly_outside_y = polygon_external_avoid$y 
	if (missing(color_ext) == FALSE) {
		polygon(glob_poly_outside_x, 
				glob_poly_outside_y, 
				border = NA, 
				lwd = 2.5, 
				col = addalpha(color_ext, alpha_ext))	
		} else {
			polygon(glob_poly_outside_x, 
					glob_poly_outside_y, 
					border = NA, 
					lwd = 2.5, 
					col = "white")	
			}
	# Plotting the overall polygon.
	if (missing(diff) == FALSE) {
		polygon(glob_poly_x, glob_poly_y, 
				border = color_border, 
				lwd = 1.5, 
				lty = lty_border, 
				col = NA)
		} else {
			polygon(glob_poly_x, glob_poly_y, 
					border = color_border, 
					lwd = 1.5, 
					lty = lty_border, 
					col = NA)	
			}
	}	
	
plot_footprint_single_polygons <- function(polygons, 
										   color_background, 
										   alpha_ext) {
	if (length(polygons) != 0) {					   
		## This function plots the polygons corresponding to each
		## field observed by HST, the overall shape of the field
		## and a few other graphical elements.
		num_polygons = length(polygons)
		# Splitting the variables to draw the single polygons.
		x_polya = list(); y_polya = list()
		x_polyb = list(); y_polyb = list()
		for (i in 1:num_polygons) {
			x_polya[[i]] = polygons[[i]]$a$x; y_polya[[i]] = polygons[[i]]$a$y 
			x_polyb[[i]] = polygons[[i]]$b$x; y_polyb[[i]] = polygons[[i]]$b$y 	
			polygon(x_polya[[i]], y_polya[[i]], border = "black", 
					lwd = 0.1, lty = 2, col = addalpha(color_background, alpha_ext))
			#text(mean(x_polya[[i]]), mean(y_polya[[i]]), paste(i, "-a", sep = ""))		
			polygon(x_polyb[[i]], y_polyb[[i]], border = "black", 
					lwd = 0.1, lty = 2, col = addalpha(color_background, alpha_ext))				
			#text(mean(x_polyb[[i]]), mean(y_polyb[[i]]), paste(i, "-b", sep = ""))		
			}
		}		
	}	
		
plot_bars <- function(radeg, decdeg, factor_dec_bars, flag_points = FALSE) {
	# This function plots symmetrical bars around points
	# representing a third quantity (for example, the 
	# velocity dispersion) on the scatterplot.
	color_bars = rep("black", length(radeg))
	color_bars[which(factor_dec_bars > 0)] = "darkblue"
	color_bars[which(factor_dec_bars <= 0)] = "red"
	factor_dec_bars_plot = abs(factor_dec_bars)
	arrows(radeg, decdeg - factor_dec_bars_plot, 
		   radeg, decdeg + factor_dec_bars_plot,      
		   pch = 10, 
		   lwd = 2, 
		   code = 3,
		   length = 0.1,
		   col = color_bars)      
	if (flag_points == TRUE) {
		points(radeg, decdeg, 
			   pch = 20, 
			   cex = 2, 
			   col = color_bars)	
		}	   	 	
	}		
	
random_powlaw <- function(n, 
						  x0 = 0, 
						  x1 = 5, 
						  a = 0.5, 
						  b = -1) {
	y = runif(n)
	norm_random = (b + 1)/(x1^(b + 1) - x0^(b + 1))
	random_x = ((y*(b + 1))/norm_random + x0^(b + 1))^(1/(b + 1))
	}
		
random_broken_powlaw <- function(pl_left, 
								 pl_right, 
								 npoints_left, 
								 npoints_right, 
								 x0 = 0, 
								 xmid = 1, 
								 x1 = 5) {
	# This function extracts points from a broken-powerlaw function
	# given two powerlaws defined by lm(). npoints_left will be 
	# extracted from the left powerlaw and npoints_right will be
	# extracted from the right powerlaw.
	# Left powerlaw.								
	yleft = runif(npoints_left)
	bleft = pl_left$coefficients[2]
	norm_random_left = (bleft + 1)/(xmid^(bleft + 1) - x0^(bleft + 1))
	random_x_left = ((yleft*(bleft + 1))/norm_random_left + x0^(bleft + 1))^(1/(bleft + 1))
	# Right powerlaw.
	yright = runif(npoints_right)
	bright = pl_right$coefficients[2]
	norm_random_right = (bright + 1)/(x1^(bright + 1) - xmid^(bright + 1))
	random_x_right = ((yright*(bright + 1))/norm_random_right + xmid^(bright + 1))^(1/(bright + 1))
	# Wrapping things together.	
	random_x = list(random_x_left, random_x_right)										
	}	
	
plot.normal.components <- function(mixture, 
								   component.number,
								   ...) {
	# This function plots different components of a gaussian 
	# mixture model.
	curve(mixture$lambda[component.number]*
          dnorm(x, mean = mixture$mu[component.number],
          sd = mixture$sigma[component.number]), add = TRUE, 
          ...)
	}		
	
plot_surfbrightnesses <- function(col_surf = c("brown", "magenta", "orange"), lwd_surf = 1.5) {
	# This function plots the surface brightness profiles.
	# from H and K images of the galaxy over the plots.
	plot(stepfun_surfbright_h,
	 	 do.points = FALSE,
	 	 add = TRUE, 
	 	 lty = 1, lwd = lwd_surf,
	 	 col = col_surf[1])
	plot(stepfun_surfbright_k,
		 do.points = FALSE,
		 add = TRUE, 
		 lty = 1, lwd = lwd_surf, 
		 col = col_surf[2])
	}
		
binnings_theta <- function(theta_binwidth_numlow = 4, 
						   theta_binwidth_numup = 6,
						   single_binning = 0) {
	# This function evaluates angular bins extremes. It has three arguments:
	# 	theta_binwidth_numlow 	->	min number of angular bins;
	# 	theta_binwidth_numup 	->	max number of angular bins;	
	# 	single_binning			-> 	flag. For 0 (default) a list of different
	#								bins is produced. For 1, only one
	#								angular binning is produced (with default
	#								number of bins = theta_binwidth_numlow). 
	if (single_binning == 0) {
		bins_theta_polar = list()
		theta_binwidths_nums = seq(theta_binwidth_numlow, theta_binwidth_numup)
		for (i in 1:length(theta_binwidths_nums)) {
			bins_theta_polar[[i]] = seq(from = -pi, 
										to = pi, 
										length.out = theta_binwidths_nums[i])
			}
		} else {
			bins_theta_polar = seq(from = -pi, 
								   to = pi, 
								   length.out = theta_binwidth_numlow)
			}
	return(bins_theta_polar)		
	}	
	
plot_arrows <- function(ra_start, 
						ra_end, 
						dec_start, 
						dec_end) {
	# This function plots the arrows indicating the directions of N and W
	# in the spatial scatterplots.
	# Arguments:
	#	ra_start		->	starting point of the horizontal arrow.
	#	ra_end			->	ending point of the horizontal arrow.
	#	dec_start		->	starting point of the vertical arrow.
	#	dec_end			->	ending point of the vertical arrow.			
	arrows(ra_start, dec_start, 
	   	   ra_end, dec_start, 
	   	   col = "black", 
	       length = 0.1, 
	       lwd = 2)
	arrows(ra_start, dec_start, 
	       ra_start, dec_end, 
	       col = "black", 
	       length = 0.1, 
	       lwd = 2)
	text(x = ra_end, y = dec_start, 
	     "W", 
	     pos = 4, 
	     col = "black")
	text(x = ra_start, y = dec_end, 
	     "N", 
	     pos = 3, 
	     col = "black")	   
	}

plot_arrows_new <- function(position_name, 
							ra_length, 
							dec_length) {
	# This function plots the arrows indicating the directions of N and W
	# in the spatial scatterplots with a different set of arguments than
	# the older function 'plot_arrows'.
	# Arguments:
	#	position_name	->	the name of the position where the arrow
	#						will be placed (usually "bottomleft").
	#	ra_length		->	length of the horizontal arrow.
	#	dec_length		->	length of the vertical arrow.
	position <- legend(position_name,
					   "Arrows", 
					   bty = "n",
					   inset = 0.05, 
					   plot = FALSE)
	arrows(position$rect$left, position$rect$top - position$rect$h, 
	   	   position$rect$left - ra_length, position$rect$top - position$rect$h, 
	   	   col = "black", 
	       length = 0.1, 
	       lwd = 2)
	arrows(position$rect$left, position$rect$top - position$rect$h, 
	       position$rect$left, position$rect$top - position$rect$h + dec_length, 
	       col = "black", 
	       length = 0.1, 
	       lwd = 2)
	text(x = position$rect$left - ra_length, y = position$rect$top - position$rect$h, 
	     "W", 
	     pos = 4, 
	     col = "black")
	text(x = position$rect$left, y = position$rect$top - position$rect$h + dec_length, 
	     "N", 
	     pos = 3, 
	     col = "black")	   
	}		
	
plot_scale <- function(position_name, 
					   num_arcseconds, 
					   inset_value) {
	# This function plots a segment showing the length corresponding 
	# to the num_arcseconds arcseconds.
	#
	# Arguments:
	#	position_name	->	position of the segment (like in command "legend").
	#	num_arcseconds	->	length of the segments expresses in number of 
	#						arcseconds.
	#	inset_value		->	inset value
	#
	# Establishing the plot coordinates corresponding to the 
	# position in "position_name".
	position <- legend(position_name,
					   "A", 
					   bty = "n",
					   inset = inset_value, 
					   plot = FALSE)	
	# Plotting the segment corresponding to
	# the number of arcseconds required.
	if (grepl("right", position_name) == TRUE) {
		start_segment_y = position$rect$top - position$rect$h/2
		start_segment_x = position$rect$left + position$rect$w/2
		end_segment_x = start_segment_x + (num_arcseconds/3600) # expressed in degrees.
		start_text_x = start_segment_x + (num_arcseconds/3600)/2
		} else {
			start_segment_y = position$rect$top - position$rect$h/2
			start_segment_x = position$rect$left + position$rect$w/2
			end_segment_x = start_segment_x - (num_arcseconds/3600) # expressed in degrees.
			start_text_x = start_segment_x - (num_arcseconds/3600)/2
			}
	segments(start_segment_x, 
			 start_segment_y, 
			 end_segment_x, 
			 start_segment_y, 
			 col = "black",
			 lwd = 2,
			 lty = 1)	
	start_text_y = position$rect$top - 0.8*position$rect$h	
	string_angular_unit = " arcsecs"
	num_angular_unit = num_arcseconds
	if (num_arcseconds >= 60 & num_arcseconds < 3600) {
		num_angular_unit = num_arcseconds/60
		string_angular_unit = " arcmin"
		} else if (num_arcseconds >= 3600) {
			num_angular_unit = num_arcseconds/3600
			string_angular_unit = " deg"
			}
	text(start_text_x, start_text_y, 
		 paste(num_angular_unit, string_angular_unit, sep = "")) 			   
	}		
	
plot_labels <- function(position_name,
						sample_name,
						inset_value,
						cex_value,
						k_value,
						sigma_flag, 
						sigma_value, 
						col_background = "white") {	
	# Function that plots the labels of the positional
	# scatterplots, density and residual maps, given a
	# position.
	#
	# Arguments:
	#	position_name		->		"bottomleft","bottomright",etc.
	#	sample_name			-> 		string to appear in the first
	#								row of the label. For example "Blue GCs"
	#	inset_value			-> 		value of the inset parameter of legend command.
	#	cex_value			->		font size of the text.
	#	k_value				->		value of the K parameter of the KNN method.
	#	sigma_flag			->		flag (1 or 0): if 1, a third row will be added
	#								to the label, containing the sigma value
	#	sigma_value			->		value of the sigma for the sigma maps.	
	#	col_background		-> 		Background color of the box containing 
	#								the text.
	position <- legend(position_name,
					   sample_name, 
					   bty = "n",
					   cex = cex_value,
					   inset = inset_value, 
					   plot = FALSE)
	if (sigma_flag == 1) {
		left_bottom = position$rect$top - 2.4*position$rect$h
		} else {
			left_bottom = position$rect$top - 1.7*position$rect$h
			}	
	rect(position$rect$left, 
		 left_bottom, 
		 position$rect$left + 0.8*position$rect$w, 
		 position$rect$top, 
		 col = col_background)
	text(x = position$rect$left, 
		 y = position$rect$top - 0.5*position$rect$h,  	   
		 sample_name, 
		 pos = 4, 
		 cex = cex_value)
	text(x = position$rect$left, 
		 y = position$rect$top - 1.2*position$rect$h,  	   
		 paste("K = ", k_value, sep = ""),
		 pos = 4, 
		 cex = cex_value)	
	if (sigma_flag == 1) {	 	
		text(x = position$rect$left, 
			 y = position$rect$top - 1.9*position$rect$h,  	   
			 paste(">", sigma_value, sep = ""),
			 pos = 4, 
			 cex = cex_value)		
		}	 	  						
	}

plot_labels_noK <- function(position_name,
							sample_name,
							inset_value,
							cex_value,
							col_background = "white") {	
	# Function that plots the labels of the positional
	# scatterplots, density and residual maps, given a
	# position.
	#
	# Arguments:
	#	position_name		->		"bottomleft","bottomright",etc.
	#	sample_name			-> 		string to appear in the first
	#								row of the label. For example "Blue GCs"
	#	inset_value			-> 		value of the inset parameter of legend command.
	#	cex_value			->		font size of the text.
	#	col_background		-> 		Background color of the box containing 
	#								the text.
	position <- legend(position_name,
					   sample_name, 
					   bty = "n",
					   cex = cex_value,
					   inset = inset_value, 
					   plot = FALSE)
	left_bottom = position$rect$top - 1.7*position$rect$h
	rect(position$rect$left, 
		 left_bottom, 
		 position$rect$left + 0.8*position$rect$w, 
		 position$rect$top, 
		 col = col_background)
	text(x = position$rect$left, 
		 y = position$rect$top - 0.5*position$rect$h,  	   
		 sample_name, 
		 pos = 4, 
		 cex = cex_value)
	}

poligonify <- function(x_open_polygon, 
					   y_open_polygon) {
	# This function produces a closed polygon from an open broken curve.
	# Arguments:
	#	x_open_polygon	->	variable containing the x coordinates of the curve.
	#	y_open_polygon	->	variable containing the y coordinates of the curve.
	x_polygon = c(0, x_open_polygon[order(y_open_polygon)], 0)
	y_polygon = c(sort(y_open_polygon)[1], sort(y_open_polygon), 
				  sort(y_open_polygon)[length(y_open_polygon)])
	return(cbind(x_polygon, y_polygon))
	}
	
externalify <- function(x_coords, 
						y_coords) {
	# This function produces the boundaries of a very large box containing
	# the whole Virgo cluster for drawing purposes.	
	x_coords = x_coords[-length(x_coords)]
	y_coords = y_coords[-length(y_coords)]
	idx_max_y_coords = which(y_coords == max(y_coords))
	x_coords_shuffled = c(x_coords[idx_max_y_coords[1]], x_coords[(idx_max_y_coords[1] + 1):length(y_coords)], 
						  x_coords[1:(idx_max_y_coords[1] - 1)], x_coords[idx_max_y_coords[1]])				
	y_coords_shuffled = c(y_coords[idx_max_y_coords[1]], y_coords[(idx_max_y_coords[1] + 1):length(y_coords)], 
						  y_coords[1:(idx_max_y_coords[1] - 1)], y_coords[idx_max_y_coords[1]])				
	x_coords_outer = c(x_coords_shuffled, x_coords_shuffled[1], 100, 100, 200, 200, x_coords_shuffled[1], x_coords_shuffled[1])
	y_coords_outer = c(y_coords_shuffled, 30, 30, 0, 0, 30, 30, y_coords_shuffled[1])		
	return(cbind(x_coords_outer, y_coords_outer))			
	}		

externalify2 <- function(x_coords1, 
						 y_coords1, 
						 x_coords2, 
						 y_coords2) {
	# This function produces the boundaries of a very large box containing
	# the whole Virgo cluster for drawing purposes, with two large structures.	
	# First structure.
	#x_coords1 = chull_coords[[1]][, 1]
	#y_coords1 = chull_coords[[1]][, 2]
	x_coords1 = x_coords1[-length(x_coords1)]
	y_coords1 = y_coords1[-length(y_coords1)]
	idx_max_y_coords1 = which(y_coords1 == max(y_coords1))
	x_coords_shuffled1 = c(x_coords1[idx_max_y_coords1[1]], x_coords1[(idx_max_y_coords1[1] + 1):length(y_coords1)], 
						   x_coords1[1:(idx_max_y_coords1[1] - 1)], x_coords1[idx_max_y_coords1[1]])				
	y_coords_shuffled1 = c(y_coords1[idx_max_y_coords1[1]], y_coords1[(idx_max_y_coords1[1] + 1):length(y_coords1)], 
						   y_coords1[1:(idx_max_y_coords1[1] - 1)], y_coords1[idx_max_y_coords1[1]])				
	#x_coords_outer1 = c(x_coords_shuffled1, x_coords_shuffled1[1], 100, 100, 200, 200, x_coords_shuffled1[1], x_coords_shuffled1[1])
	#y_coords_outer1 = c(y_coords_shuffled1, 30, 30, 0, 0, 30, 30, y_coords_shuffled1[1])		
	# Second structure.	
	#x_coords2 = chull_coords[[2]][, 1]
	#y_coords2 = chull_coords[[2]][, 2]
	x_coords2 = x_coords2[-length(x_coords2)]
	y_coords2 = y_coords2[-length(y_coords2)]
	idx_min_y_coords2 = which(y_coords2 == max(y_coords2))		
	x_coords_shuffled2 = c(x_coords2[idx_min_y_coords2[1]:length(y_coords2)], x_coords2[1:(idx_min_y_coords2[1] - 1)], 
						   x_coords2[idx_min_y_coords2[1]])				
	y_coords_shuffled2 = c(y_coords2[idx_min_y_coords2[1]:length(y_coords2)], y_coords2[1:(idx_min_y_coords2[1] - 1)], 
						   y_coords2[idx_min_y_coords2[1]])		
	#x_coords_outer2 = c(x_coords_shuffled1, x_coords_shuffled1[1], 100, 100, 200, 200, x_coords_shuffled1[1], x_coords_shuffled1[1])
	#y_coords_outer2 = c(y_coords_shuffled1, 30, 30, 0, 0, 30, 30, y_coords_shuffled1[1])					   	
	# Combining the two structures.
	x_coords_outer = c(x_coords_shuffled1, x_coords_shuffled1[1], 100, 100, x_coords_shuffled2[1], 
					   x_coords_shuffled2, x_coords_shuffled2[1], 200, 200, x_coords_shuffled1[1], 
					   x_coords_shuffled1[1])
	y_coords_outer = c(y_coords_shuffled1, 30, 30, 0, 0, 
					   y_coords_shuffled2, 0, 0, 30, 30, y_coords_shuffled1[1])
	return(cbind(x_coords_outer, y_coords_outer))				   
	}							
	
coplot_annuli <- function(galaxy_center, 
		   				  range_ra, 
		   				  range_dec, 
		   				  max_radius_fit_arcmin, 
		   				  binwidth_r_gc, 
		   				  galaxy_isodiameter_ratio, 
		   				  galaxy_angle, 
		   				  num_knots_grid) {
    # Plotting the grid used to calculate the area of the intersection
    # of the ellipses with the HST field.
  	#abline(v = seq(min(range_ra), max(range_ra), length.out = num_knots_grid), 
    #	   col = "grey70")
    #abline(h = seq(min(range_dec), max(range_dec), length.out = num_knots_grid), 
    #	   col = "grey70")	    	 
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
					"white", 0)	
	# Shape of the first galaxy.
	# Circles (or ellipses) corrisponding to the radial bins.	 
	#x_seeds = seq(min(range_ra), max(range_ra), length.out = num_knots_grid)	
	#y_seeds = seq(min(range_dec), max(range_dec), length.out = num_knots_grid)	
	# Defining the grid.
	#xy_grid = expand.grid(x_seeds, y_seeds)	
	bins_r_polar = seq(from = 0, to = max_radius_fit_arcmin/60, by = binwidth_r_gc)							
	colors_ellipses = sample(colours(), length(bins_r_polar) - 1)
	for (j in 2:length(bins_r_polar)) {
		#points(xy_grid[intersections[[j - 1]]$indices, 1], 
		#	   xy_grid[intersections[[j - 1]]$indices, 2], 
		#	   col = colors_ellipses[j - 1], 
		#	   cex = 2, 
		#	   pch = 20)	
		draw.ellipse(galaxy_center[1], galaxy_center[2], 
					 a = bins_r_polar[j], 
					 b = bins_r_polar[j]/(10^(galaxy_isodiameter_ratio)), 
					 angle = 270 - galaxy_angle,
					 nv = 1000, 
					 border = "grey60",
					 lwd = 0.2, 
					 deg = TRUE)			 
		}	
	}	

arrow_starting_point <- function(x_range, 
							     y_range, 
							     side) {
	# Generates the starting point of an arrow given 
	# the range covered by the x and y sides of a box
	# and the side where the starting point should be
	# located (top, bottom, right, left).
	xy_starting_point = vector(mode = "numeric", length = 2) 
	if (side == "top") {
		xy_starting_point[1] = (max(x_range) - min(x_range))/2 + min(x_range)
		xy_starting_point[2] = max(y_range)
		}
	if (side == "bottom") {
		xy_starting_point[1] = (max(x_range) - min(x_range))/2 + min(x_range)
		xy_starting_point[2] = min(y_range)
		}
	if (side == "right") {
		xy_starting_point[1] = min(x_range)
		xy_starting_point[2] = (max(y_range) - min(y_range))/2 + min(y_range)
		}
	if (side == "left") {
		xy_starting_point[1] = max(x_range)
		xy_starting_point[2] = (max(y_range) - min(y_range))/2 + min(y_range)
		}
	return(xy_starting_point)						     
	}		

plot_variable_isophotes <- function(radeg_gc, decdeg_gc, 
									range_gc, range_dec,
									shape_class_lum, 
       								intersections, 
       								num_knots_grid, 
       								radii_values, 
       								galaxy_center, 
       								radii_values_sim) {
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.75, 0.75, 0.20, 0.20), mgp = c(1.8, 0.2, 0), las = 1)	
	plot(radeg_gc, decdeg_gc, 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
	     xlim = c(54.76, 54.48), 
    	 ylim = c(-35.58, -35.36), 
    	 type = "n")
    # Plotting the GCs positions.
	points(radeg_gc, decdeg_gc,      
		   pch = shape_class_lum, 
		   cex = cex_gc/3, 
		   col = "black") 	 
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
					"white", 0)		
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length, dec_length)	 	   
	# Shape of the D25 isophote of the galaxy.
	plot_galaxies(galaxy_center, galaxy_isodiameter, galaxy_isodiameter_ratio, galaxy_angle)
	# Circles (or ellipses) corrisponding to the radial bins.	 
	x_seeds = seq(min(range_ra), max(range_ra), length.out = num_knots_grid)	
	y_seeds = seq(min(range_dec), max(range_dec), length.out = num_knots_grid)	
	# Defining the grid.
	xy_grid = expand.grid(x_seeds, y_seeds)								
	colors_ellipses = sample(colours(), length(radii_values))
	for (j in 2:length(radii_values)) {
		draw.ellipse(galaxy_center[1], galaxy_center[2], 
					 a = radii_values[j],
					 b = interpolate_function_b(radii_values[j]),
					 angle = interpolate_function_angle(radii_values[j]),
					 nv = 1000, 
					 border = "grey60",
					 lwd = 2, 
					 deg = TRUE)	
		if (missing(radii_values_sim) == FALSE) {
			draw.ellipse(galaxy_center[1], galaxy_center[2], 
						 a = radii_values_sim[j],
						 b = interpolate_function_b(radii_values_sim[j]),
						 angle = interpolate_function_angle(radii_values_sim[j]),
						 nv = 1000, 
						 border = "grey60",
						 lwd = 2, 
						 lty = 2,
						 deg = TRUE)
			} else {
				points(xy_grid[intersections[[j - 1]]$indices, 1], 
					   xy_grid[intersections[[j - 1]]$indices, 2], 
					   col = colors_ellipses[j - 1], 
					   cex = 0.5, 
					   pch = 20)	
				}	 		 
		}	
	nty = par()$yaxp
	nty[3] = nty[3]*5
	ntx = par()$xaxp
	ntx[3] = ntx[3]*5
	par(xaxp = ntx, yaxp = nty, tcl = 0.2)
	axis(1, label = FALSE)
	axis(2, label = FALSE)
	axis(3, label = FALSE)
	axis(4, label = FALSE)
	}
	
radial_range_concentric_annuli <- function(a_small, b_small, pa_small,  # Major and minor axes and P.A. of the larger isophote
										   a_large, b_large, pa_large,	# Major and minor axes and P.A. of the smaller isophote
										   angle) {						# Angular position of the source [radians]
	diff_angle_small = angle - pi*pa_small/180		# Angular position relative to the P.A. of the smaller osophote
	diff_angle_large = angle - pi*pa_large/180
	range_radius = a_large*b_large/sqrt((a_large*sin(diff_angle_large))^2  + (b_large*cos(diff_angle_large))^2) - 
				   a_small*b_small/sqrt((a_small*sin(diff_angle_small))^2  + (b_small*cos(diff_angle_small))^2)										   	
	}		
	
plot_line_grid <- function(range_ra, 
						   range_dec,
						   binwidth_ra,
						   binwidth_dec = binwidth_ra) {
	abline(v = seq(min(range_ra), max(range_ra), by = binwidth_ra), 
	   	   col = "darkgray", 	
	       lwd = 0.3, 
	       lty = 2)	     
	abline(h = seq(min(range_dec), max(range_dec), by = binwidth_dec), 
	       col = "darkgray", 
	       lwd = 0.3, 
	       lty = 2)					   
	}						   									   	
						     	