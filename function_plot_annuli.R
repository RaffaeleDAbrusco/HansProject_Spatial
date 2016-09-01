## Function that creates the plot of the positions of the GCs in R.A. vs Dec.
## with the elliptical shapes used to evaluate the observed radial density of GCs. 
##
## Written by R. D'Abrusco
##
## Last modified: 31/8/16.

plot_annuli <- function(galaxy_name, 
						path_plot, 
						radeg_gc, 
						decdeg_gc, 
						range_ra,
						range_dec,
						color_class_color, 
       					intersections, 
       					num_knots_grid, 
       					polygons_annuli, 
       					polygons_external_annuli, 
       					polygons_external_avoid_annuli, 
       					ra_length_arrow, 
       					dec_length_arrow, 
       					galaxy1_center, 
       					galaxy1_isodiameter, 
       					galaxy1_isodiameter_ratio, 
       					galaxy1_angle, 
       					glob_poly_back_ra_annuli, glob_poly_back_dec_annuli, 
       					back_color_annuli, back_transparency_annuli) {
	pathpdfannuli = paste(path_plot, "positions_ellipse.pdf", sep = "")
	pdf(pathpdfannuli, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(rev(radeg_gc), decdeg_gc, 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
		 xlim = rev(range_ra), 
		 ylim = range_dec,
    	 type = "n", 
    	 axes = F)
    # Plotting the grid used to calculate the area of the intersection
    # of the ellipses with the HST field.
    abline(v = seq(min(range_ra), max(range_ra), length.out = num_knots_grid), 
    	   col = "grey70")
    abline(h = seq(min(range_dec), max(range_dec), length.out = num_knots_grid), 
    	   col = "grey70")	    	 
	# Drawing the footprint of the observations. 
	plot_footprint(polygons_annuli, polygons_external_annuli, polygons_external_avoid_annuli, 1, 
			   	   "white", 255)	
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)	
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, galaxy1_angle)
	# Circles (or ellipses) corrisponding to the radial bins.	 
	x_seeds = seq(max(range_ra), min(range_ra), length.out = num_knots_grid)	
	y_seeds = seq(min(range_dec), max(range_dec), length.out = num_knots_grid)	
	# Defining the grid.
	xy_grid = expand.grid(x_seeds, y_seeds)								
	colors_ellipses = sample(colours(), length(bins_r_polar) - 1)
	for (j in 2:length(bins_r_polar)) {
		points(xy_grid[intersections[[j - 1]]$indices, 1], 
			   xy_grid[intersections[[j - 1]]$indices, 2], 
			   col = colors_ellipses[j - 1], 
			   cex = 1, 
			   pch = 20)	
		draw.ellipse(galaxy1_center[1], galaxy1_center[2], 
					 a = bins_r_polar[j], 
					 b = bins_r_polar[j]/(10^(galaxy1_isodiameter_ratio)), 
					 angle = 270 - galaxy1_angle,
					 nv = 1000, 
					 border = "grey60",
					 lwd = 2, 
					 deg = TRUE)			 
		}	
	# Plotting the points.
	points(radeg_gc, decdeg_gc,      
		   pch = 20, 
		   cex = cex_gc/3, 
		   col = color_class_color)
	# Plotting the area used to estimate the background.
	polygon(glob_poly_back_ra_annuli, glob_poly_back_dec_annuli, 
			col = addalpha(back_color_annuli, back_transparency_annuli))	  
	# Embellishing the axes.	
	axis(1, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.1), 
			label = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.1), 
			tcl = 0.5)
	axis(1, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.025), 
			label = FALSE, 
			tcl = 0.2)
	axis(2, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.1), 
			label = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.1), 
			tcl = 0.5)
	axis(2, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.025), 
			label = FALSE, 
			tcl = 0.2)
	axis(3, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.1), 
			label = FALSE, 
			tcl = 0.5)
	axis(3, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.025), 
			label = FALSE, 
			tcl = 0.2)
	axis(4, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.1), 
			label = FALSE, 
			tcl = 0.5)
	axis(4, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.025), 
			label = FALSE, 
			tcl = 0.2)
	dev.off()	   
	}
