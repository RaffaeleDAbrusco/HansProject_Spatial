## This program characterizes the spatial distribution of 
## globular clusters (GCs) from the VEGAS (ground-based) data
## around the galaxy NGC3115 with KNN and simulations.
##
## Written by R. D'Abrusco
##
## Last modified: 31/08/16.

## Importing required libraries.
suppressMessages(library(rgl, quietly = TRUE))
suppressMessages(library(gplots, quietly = TRUE))
suppressMessages(library(HiveR, quietly = TRUE))	# Used to convert between polar and cartesian 
													# coordinates.
suppressMessages(library(shape, quietly = TRUE))	# Used to draw shapes in the plots.
suppressMessages(library(Hmisc, quietly = TRUE))	# Used to add subplots to the plots.
suppressMessages(library(plotrix, quietly = TRUE))	# Used to plot other shapes and lines.
suppressMessages(library(gplots, quietly = TRUE))	# Used to plot error bars in scatterplots.
suppressMessages(library(fields, quietly = TRUE))	# Used to smooth the 2D distributions.
suppressMessages(library(FITSio, quietly = TRUE))	# Used to save smoothed 2D maps as fits files.
suppressMessages(library(emdbook, quietly = TRUE))	# Used to created logspaced sequences.
suppressMessages(library(mixtools, quietly = TRUE))	# Used to fit gaussian mixtures to 1D data.	
suppressMessages(library(mclust, quietly = TRUE))	# Used to estimate 2D density with an alternate method.	
suppressMessages(require(MASS))						# Required to create the 2D density maps with kde2d.																
#suppressMessages(require(celestial))				# Required to calculate coverage with proposed pointings.	
source("./function_raycasting.R")					# Loading functions to determine grid points within 	
													# the overall field region.
source("./function_auxiliary.R")					# Loading the auxiliary functions.		
source("./function_structures_radius.R")			# Loading the function producing the
													# properties of merged spatial structures.	
source("./function_generate_points_ellipse.R")		# Loading the function that produces
													# simulated elliptical distributions of
													# points with a given geometry.															
source("./function_generate_points_homogeneous.R")	# Loading the function that produces
													# random homogeneous simulated coordinates
													# with given density and region on the sky.	
source("./function_plot_annuli.R")					# Loading the function that produces
													# a diagnostic plot of the distribution of 
													# candidate GCs with the elliptical annuli
													# and the grid used to calculate approximately
													# the area of the elliptical annuli.
source("./function_image_residual_maps.R")			# Produces plottable images from vector of sigma values of
													# the residual maps. 
source("./function_image_residual_maps_blurred.R")	# Produces plottable images from vector of sigma values of
													# the smoothed residual maps. 												
cat("So far so good! Libraries loaded!\n")

## Setting the main technical parameters of the program.
sim_flag = "yes"					    # Simulations flag: if "yes", simulations will 
										# be performed, otherwise only plots of observed
										# data will be produced.
plotonly_flag = "yes"					# Plot-only flag: if "yes", the data will be read
										# from the session image and only the plots will
										# be created: no final file will be written.
num_simulations = 3						# Number of simulations performed.
sim_type = "observed"					# Type of model to be used for the simulation.
plot_pixels = "no"						# Pixels plot flag: if "yes", the plots of the 
										# histograms of the density value distribution for
										# each pixel will be created.
profiles_angslices_flag = "no"			# Radial profiles in angular slices flag: if "yes", 
										# the plots of the radial profile of the GCs in 
										# angular slices will be created. 
chull_flag = "no"						# Convex hull flag: if "yes", coordinates of the
										# convex hull(s) of the large structure(s) will 
										# be retrieved and plotted on the "average...._thresholds_sigmas.pdf"
										# plots.	
typechull_flag = "large_intermediate"	# Flag indicating the type of outlines to be plotted.																								
color_red_gc = "red"					# Color used to display/plot red GCs.
color_blue_gc = "darkblue"				# Color used to display/plot blue GCs.
cex_center_galaxy = 1					# Size of the symbol showing the center of the galaxy.
basic_colors_residual = c("darkgreen", "white", "orange")		# Colors used to plot the residual maps.
cat("So far so good! Parameters of the program set!\n")

## Setting the paths to the location of the files (data and plots)
## that will be useful later.
path_data = "~/Desktop/Spatial/Data/"			# Path of the folder containing the 
												# catalog listing the GCs.
root_path_plot = "~/Desktop/Spatial/Plots/" 	# Path of the folder containing the
												# plots produced in this program.

## Setting name and physical parameters for the target galaxy.
# Name.
galaxy = "ngc3115"								# Name of the primary galaxy.
secondary_galaxy = "none"						# Name of the secondary galaxy.
# Color and magnitude thresholds for the GCs distribution.
color_threshold = 0.75							# g-z color threshold(s) to select red and
											    # blue GCs.
mag_z_threshold = 23							# z mag threshold(s) separating luminous 
												# and faint GCs.
cat("So far so good! Name of the galaxy and photometric parameters of its GC system set!\n")

# Determining (and creating, if necessary) the folders that 
# will contain the plots.
path_plot = paste(root_path_plot, galaxy, "/GCs/", sep = "")
path_plot_simulations = paste(path_plot, "/Simulations/", sep = "")
dir.create(path_plot, showWarnings = FALSE)	
dir.create(path_plot_simulations, showWarnings = FALSE)	
cat("So far so good! Directories containing the plots determined and/or created!\n")

## Reading the basic metadata about the galaxy whose GCs spatial distribution
## will be investigated.											
metadata_primary_galaxy = read_galaxy_metadata(galaxy, path_data, galaxy)			
galaxy1_center = as.numeric(metadata_primary_galaxy[1:2])
galaxy1_isodiameter = as.numeric(metadata_primary_galaxy[3])
galaxy1_isodiameter_ratio = as.numeric(metadata_primary_galaxy[4])
galaxy1_angle = as.numeric(metadata_primary_galaxy[5])
# Reading the metadata of the secondary galaxy (if any).
if (secondary_galaxy != "none") {
	metadata_secondary_galaxy = read_galaxy_metadata(secondary_galaxy, path_data, galaxy)			
	galaxy2_center = as.numeric(metadata_secondary_galaxy[1:2])
	galaxy2_isodiameter = as.numeric(metadata_secondary_galaxy[3])
	galaxy2_isodiameter_ratio = as.numeric(metadata_secondary_galaxy[4])
	galaxy2_angle = as.numeric(metadata_secondary_galaxy[5])
	}												
cat("So far so good! Photometric parameters of the galaxy optical image read!\n")

# Reading the convex hull coordinates of the large structure(s)
if (chull_flag == "yes") {
	chull_coordinates = structures_radius(galaxy, 
										  "all", 
										  9, 
										  num_simulations, 
										  36, 
										  "observed", 
										  typechull_flag)
	cat("So far so good! Data of the convex hull(s) read!\n")										  
	}

## Setting graphical parameters that will be used later for the plots
## of the GCs spatial distribution.
# Shape of the symbols used for all the plots for all, matched and unmatched GCs. 
symbol_gc = 20; symbol_unmatched = 15; symbol_matched = 17
# Size of the symbols used for all the plots for all, matched and unmatched GCs. 
cex_gc = 1; cex_unmatched = 1; cex_matched = 1
# Shape of the symbols used for the red and blue GCs. 
symbol_red = 12; symbol_blue = 17
# Size of the symbols used for red and blue GCs. 
cex_red = 1; cex_blue = 1
# Type of the line used for matched and unmatched GCs. 
ltype_matched = 4; ltype_unmatched = 6; lwd_gc = 2
# Type of the line used for red and blue GCs. 
ltype_gc = 1; ltype_red = 3; ltype_blue = 5
# Line width of the lines used for red and blue GCs. 
lwd_red = 0.8; lwd_blue = 0.8
# Line width of the lines used for matched and unmatched GCs. 
lwd_matched = 0.8; lwd_unmatched = 0.8
# Size of the symbols used in the average 2D density residual
# maps.
size_position = 1						# Size of the symbols used for the position
										# of the GCs in the plots.
size_legend = 1							# Size of the symbols used for the legends.										
alphavalue = 35							# Alpha value used for the plots.				  										
binwidth_r_gc = 0.003125*6				# Bin width of the radius for the radial profile.
theta_binwidth_numlow = 4				# Lower limit of the angular coordinate binning.
theta_binwidth_numup = 8				# Upper limit of the angular coordinate binning.
radius_avoidance_deg = 0.0				# Radius of the central avoidance regions
										# caused by the incompleteness of the observations 
										# [deg]. 
radius_avoidance_arcmin = 0.0			# As above [arcmin].
max_radius_fit_arcmin = 20				# Maximum radial distance used to fit the 
										# radial density profiles of the GCs [arcmin]. 	
single_binning = 0						# Single binning flag.
binwidth_theta_gc = pi/4				# Bin width of the theta for the slicing of the 
										# GCs spatial distribution (in radians).
cex_labels_axes = 1						# Size of the axes labels of the plots.
cex_axes = 0.5							# Size of the axes ticks of the plots.
mgp_axes = c(1, 0.1, 0)					# Positioning of the axes elements of the plots.
flag_surfbright = 1						# Flag of the surface brightness profiles. If set 
										# to 1, the relative curvers are plotted.
nsteps_smooth = 36						# Number of steps for the smoothing of the 2D
										# KNN density.
k_vector = c(2, 5, 9)
#k_vector = c(2, 3, 4, 5, 6, 7, 8, 9)  	# Values of the K parameter for the KNN 2D estimation.
index_dens_threshold = 8				# Index of the density contour used to separate
										# the sources in the high-density regions from 
										# the other sources.
index_res_threshold = 6					# Index of the residual map contour used to separate
										# the sources in the high-density regions from 
										# the other sources.										
n_breaks_rsplit_sim = 20 				# Number of bins of the histogram used to derive
										# the simulated distributions of GCs.						
pixel_percentile = .9					# Percentile of the distribution of 2D density
										# pixels (normalized to 1).	
num_sigmas_threshold = c(0.5, 1, 2, 3) 	# Number of sigmas used to select the pixels
										# with the peaks of the residuals distribution.	
num_knots_areas = 100					# Number of knots of the grid used to evaluate
										# the areas of the regions observed taking into	
										# account the fields of views, the shape of the  
										# elliptical annuli and secondary galaxies, if present.	
points_average_residuals = "no"			# Flag determining whether to plot the position
										# of the sources in the average residuals plots.
range_ra = c(151 - 0.25, 151.6 + 0.25)	# Range along the R.A. for plots. 
range_dec = c(-8.1 - 0.25, -7.5 + 0.25)	# Range along the Dec. for plots.
ra_length_arrow = 0.1					# Length of the horizontal arrow [deg].
dec_length_arrow = 0.1					# Length of the vertical arrow [deg].
cat("So far so good! Graphical parameters of the program set!\n")

## Loading the file containing the data of the session if 
## sim_flag is equal to "no".
if (sim_flag == "no" & plotonly_flag == "yes") {
	name_file_session = paste(path_data, galaxy, "/Sessions/gc_nsims_", 
							  num_simulations, "_nsteps_", nsteps_smooth, 
							  "_simtype_", sim_type, ".Rdata", sep = "")
	load(name_file_session)
	cat("So far so good! Image of a previous session imported!\n")
	# Setting (again) the values of the main program parameters.
	sim_flag = "no"					# Simulations flag: if "yes", simulations will 
									# be performed, otherwise only plots of observed
									# data will be produced.
	plotonly_flag = "yes"			# Plot-only flag: if "yes", the data will be read
									# from the session image and only the plots will
									# be created. 
	# Put here graphical parameters of the plots that need be changed.	
	index_dens_threshold = 8			# Index of the density contour used to separate
										# the sources in the high-density regions from 
										# the other sources.																										  																																														
	index_res_threshold = 6				# Index of the residual map contour used to separate
										# the sources in the high-density regions from 
										# the other sources.	
	points_average_residuals = "no"		# Flag determining whether to plot the position
										# of the sources in the average residuals plots.
	max_radius_fit_arcmin = 20			# Maximum radial distance shown in the 
										# plot of the radial density profile. 
	num_sigmas_threshold = c(0.5, 1, 2, 3) 	# Number of sigmas used to select the pixels
											# with the peaks of the residuals distribution.											
	range_ra = c(151 - 0.2, 151.6 + 0.2)	# Range along the R.A. for plots. 
	range_dec = c(-8.1 - 0.2, -7.5 + 0.2)	# Range along the Dec. for plots.
	ra_length_arrow = 0.1					# Length of the horizontal arrow [deg].
	dec_length_arrow = 0.1					# Length of the vertical arrow [deg].
	cat("So far so good! Main parameters of the program set again!\n")
	}	

## Reading the main catalog of GCs in the galaxy considered. 
filename_catalog_gc = paste(path_data, galaxy, "/gc_", galaxy, ".dat", sep = "")
data_table_catalog_gc = read.table(filename_catalog_gc, skip = 1)
id_gc = data_table_catalog_gc[, 1]		    	# ID of the source.
radeg_gc = data_table_catalog_gc[, 2]			# RA in degrees.
decdeg_gc = data_table_catalog_gc[, 3]			# Dec in degrees.
gal_dist = data_table_catalog_gc[, 4]			# Galacto-centric distance.
u_gc = data_table_catalog_gc[, 5] 				# u magnitude.
erru_gc = data_table_catalog_gc[, 6] 			# Error on u magnitude.
g_gc = data_table_catalog_gc[, 7] 				# g magnitude.
errg_gc = data_table_catalog_gc[, 8] 			# Error on g magnitude.
i_gc = data_table_catalog_gc[, 9] 				# i magnitude.
erri_gc = data_table_catalog_gc[, 10] 			# Error on i magnitude.
cat("So far so good! Catalog of the first GCs list read!\n")

## Creating a few split list of variables based on factors, that will be used later.
### Matched and unmatched to LMXBs.
## Defining the color of the GCs.
gmi_gc = g_gc - i_gc
colorclass_gc = ifelse(gmi_gc > color_threshold, color_red_gc, color_blue_gc)
factor_colorclass_gc = as.factor(colorclass_gc)
radeg_gc_colorclass = split(radeg_gc, factor_colorclass_gc)
decdeg_gc_colorclass = split(decdeg_gc, factor_colorclass_gc)
id_gc_colorclass = split(id_gc, factor_colorclass_gc)
id_gc_red = id_gc_colorclass[[2]]
id_gc_blue = id_gc_colorclass[[1]]
luminclass_gc = ifelse(g_gc > 15, "black", "black")
cat("So far so good! Color and other variables defined!\n")

## Defining the polygons associated to the footprint of the data used to
## extract the candidate GCs.
glob_poly_x = c(150.8, 151.8, 151.8, 150.8, 150.8)
glob_poly_y = c(-8.3, -8.3, -7.3, -7.3, -8.3)
polygons_external = list(x = glob_poly_x, y = glob_poly_y) 
glob_poly_xy = cbind(glob_poly_x, glob_poly_y)
glob_poly_outside_x = c(150.8, 151.8, 151.8, 150.8, 150.8, 145, 145, 155, 155, 145, 145, 150.8)  				
glob_poly_outside_y = c(-8.3, -8.3, -7.3, -7.3, -8.29999, -8.29999, -5, -5, -12, -12, -8.3, -8.3)
polygons_external_avoid = list(x = glob_poly_outside_x, y = glob_poly_outside_y) 			
cat("So far so good! Polygons defining the footprints of the observed region defined!\n")

## Defining the background region and calculating the density of 
## all, red and blue candidate GCs in the background region.
range_back_ra = c(151.1, 150.8)
range_back_dec = c(-7.6, -7.3)
glob_poly_back_ra = c(rep(range_back_ra[1], 2), rep(range_back_ra[2], 2), range_back_ra[1])
glob_poly_back_dec = c(range_back_dec, range_back_dec[2], rep(range_back_dec[1], 2))
idx_gc_all_back = which(radeg_gc >= min(range_back_ra) & radeg_gc <= max(range_back_ra) & 
						decdeg_gc <= max(range_back_dec) & decdeg_gc >= min(range_back_dec))
idx_gc_red_back = which(radeg_gc >= min(range_back_ra) & radeg_gc <= max(range_back_ra) & 
						decdeg_gc <= max(range_back_dec) & decdeg_gc >= min(range_back_dec) & 
						gmi_gc > color_threshold)
idx_gc_blue_back = which(radeg_gc >= min(range_back_ra) & radeg_gc <= max(range_back_ra) & 
						 decdeg_gc <= max(range_back_dec) & decdeg_gc >= min(range_back_dec) & 
						 gmi_gc <= color_threshold)
dens_cont_all = length(idx_gc_all_back)/((max(range_back_ra) - min(range_back_ra))*(max(range_back_dec) - min(range_back_dec)))/3600
dens_cont_red = length(idx_gc_red_back)/((max(range_back_ra) - min(range_back_ra))*(max(range_back_dec) - min(range_back_dec)))/3600
dens_cont_blue = length(idx_gc_blue_back)/((max(range_back_ra) - min(range_back_ra))*(max(range_back_dec) - min(range_back_dec)))/3600
cat("So far so good! Density of contaminants from the background region calculated!\n")

## Creating the plot of the positions of the GCs.
radeg_lims = rev(c(min(radeg_gc) - min(radeg_gc)/50000, max(radeg_gc) + max(radeg_gc)/50000))
decdeg_lims = c(min(decdeg_gc) - min(decdeg_gc)/10000, max(decdeg_gc) + max(decdeg_gc)/10000)
# Type of GCs (red or blue) is encoded in the color of the symbol.
pathpdf1 = paste(path_plot, "positions.pdf", sep = "")
pdf(pathpdf1, width = 7, height = 7, paper = "special")
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
	mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
plot(radeg_gc, decdeg_gc, 
	 xlab = "Right Ascension [deg]", 
     ylab = "Declination [deg]",
     xlim = rev(range_ra), 
     ylim = range_dec, 
     pch = 20, 
     cex = cex_gc,
     col = colorclass_gc, 
     axes = F)
# Drawing the footprints. 
polygons = list()
plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
			   "white", 255)	
# Background.
polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  
# Drawing arrows indicating the North and West directions. 	 
plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)
# Shape of the primary galaxy.
plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
			  galaxy1_angle, cex_center = cex_center_galaxy)
# Grid of lines.
plot_line_grid(range_ra, range_dec, 0.05)
# Plotting the angular scale.
plot_scale("topleft", 180, 0)		
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
cat("So far so good! Plot of the positions of the candidate GCs created!\n")

## Histograms of the g-i color of all GCs for 
## different slices of galactocentric radii.
# Histogram of the g-i color of all GCs (all GCs).
gaussian.mixtures = normalmixEM(gmi_gc, 
								k = 2, 
								mu = NULL, 
								maxit = 1000, 
								epsilon = 1E-5)
pathpdfhistogmicol = paste(path_plot, "histo_gmi_color.pdf", sep = "")	
pdf(pathpdfhistogmicol, width = 7, height = 7, paper = "special")
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 1, 
	mgp = c(2, 0.3, 0), mai = c(1.02, 0.82, 0.2, 0.2), las = 1, bg = "white")
histo_all <- hist(gmi_gc, 
	 			  border = "black", 
	 			  xlab = "g-i", 
	 			  ylab = "", 
	 			  plot = TRUE, 
	 			  main = "", 
	 			  prob = TRUE, 
	 			  yaxt = "no", 
	 			  breaks = 20)
lines(seq(min(gmi_gc), max(gmi_gc), length.out = 100), 
	  exp((-(seq(min(gmi_gc), max(gmi_gc), length.out = 100) - gaussian.mixtures$mu[1])^2)/gaussian.mixtures$sigma[1]), 
	  col = color_blue_gc, 
	  lwd = 2, 
	  lty = 1)	 
lines(seq(min(gmi_gc), max(gmi_gc), length.out = 100), 
	  exp((-(seq(min(gmi_gc), max(gmi_gc), length.out = 100) - gaussian.mixtures$mu[2])^2)/gaussian.mixtures$sigma[2]), 
	  col = color_red_gc, 
	  lwd = 2, 
	  lty = 1)	 
lines(seq(min(gmi_gc), max(gmi_gc), length.out = 100), 
	  exp((-(seq(min(gmi_gc), max(gmi_gc), length.out = 100) - gaussian.mixtures$mu[2])^2)/gaussian.mixtures$sigma[2]) + 
	  exp((-(seq(min(gmi_gc), max(gmi_gc), length.out = 100) - gaussian.mixtures$mu[1])^2)/gaussian.mixtures$sigma[1]), 
	  col = "black", 
	  lwd = 4, 
	  lty = 1)		  	  			  
dev.off()
cat("So far so good! Histograms of the magnitude and color distributions of GCs created!\n")

## Converting the cartesian coordinates of the GCs to polar coordinates to determine 
## slices of the GC spatial distribution and radial profiles.
# Determining the angular and radial binnings to be
# applied to the different datasets.
galaxy_a = (10^(galaxy1_isodiameter - 1))/120
galaxy_b = ((10^(galaxy1_isodiameter - 1))/120)/10^(galaxy1_isodiameter_ratio) 
bins_r_polar = seq(from = 0, to = (max_radius_fit_arcmin/60 + binwidth_r_gc), by = binwidth_r_gc)
bins_r_polar_sim = bins_r_polar
#bins_r_polar_sim = seq(from = 0, to = max_radius_fit_arcmin/60, by = binwidth_r_gc*1.1)
# Calculating the polar coordinates of all sources in one of the 
# coordinates systems ("natural" and "galaxy" allowed so far).
bins_theta_polar = binnings_theta(theta_binwidth_numlow, 
								  theta_binwidth_numup, 
								  single_binning)
bins_theta_polar_axes = c(-pi, -3/4*pi, -pi/2, -1/4*pi, 0, 1/4*pi,pi/2, 3/4*pi, pi)		
# Determining the type of coordinate system to be used.								  
azimuthtype = "galaxy"
if (azimuthtype == "natural") {
	theta_gc = atan2((decdeg_gc - galaxy1_center[2]), (radeg_gc - galaxy1_center[1]))							  
	} else {
		rot_angle = pi/2 - galaxy1_angle*pi/180
		rot = matrix(c(cos(rot_angle), -sin(rot_angle), sin(rot_angle), cos(rot_angle)), 
					 byrow = TRUE, 2, 2)
		radegdecdeg_norm_gc = cbind(radeg_gc - galaxy1_center[1], decdeg_gc - galaxy1_center[2])
		radegdecdeg_normrotate_gc = radegdecdeg_norm_gc%*%rot
		theta_gc = atan2(radegdecdeg_normrotate_gc[, 2], radegdecdeg_normrotate_gc[, 1])
		}
r_gc = sqrt((radeg_gc - galaxy1_center[1])^2 + (decdeg_gc - galaxy1_center[2])^2)
# Splitting the GCs coordinates in red and blue sources.
radeg_gc_colorclass = split(radeg_gc, factor_colorclass_gc)
decdeg_gc_colorclass = split(decdeg_gc, factor_colorclass_gc)
radeg_gc_red = radeg_gc_colorclass[[2]]
radeg_gc_blue = radeg_gc_colorclass[[1]] 
decdeg_gc_red = decdeg_gc_colorclass[[2]]
decdeg_gc_blue = decdeg_gc_colorclass[[1]] 
theta_gc_colorclass = split(theta_gc, factor_colorclass_gc)
r_gc_colorclass = split(r_gc, factor_colorclass_gc)
theta_gc_red = theta_gc_colorclass[[2]]
theta_gc_blue = theta_gc_colorclass[[1]]
r_gc_red = r_gc_colorclass[[2]]
r_gc_blue = r_gc_colorclass[[1]]
cat("So far so good! Polar coordinates and binning(s) of the polar coordinates calculated!\n")

# Calculating the area of the elliptical annuli in arcmin2.
stuff_intersections = list()
bincenters_r_polar = (bins_r_polar[1:length(bins_r_polar)] + binwidth_r_gc/2)*60
area_gc_rsplit = vector(mode = "numeric", length(bincenters_r_polar))
area_gc_rsplit_intersections = vector(mode = "numeric", length(bincenters_r_polar))
average_radial_distance_gc_rsplit = vector(mode = "numeric", length = length(bincenters_r_polar))
a_r_polar = bins_r_polar
b_r_polar = bins_r_polar/(10^(galaxy1_isodiameter_ratio))
#load(paste(path_data, galaxy, "stuff_intersections.Rdat", sep = ""))
for (k in 1:length(bincenters_r_polar)) {
	cat(paste("		So far so good! Calculating the area of the ", k, "-th radial bin!\n", sep = ""))
	stuff_intersections[[k]] = area_intersections(rev(range_ra), 
												  range_dec,
												  num_knots_areas,
												  glob_poly_xy,
												  galaxy1_center,
												  90 - galaxy1_angle,	
												  c(a_r_polar[k + 1], b_r_polar[k + 1]), 
												  c(a_r_polar[k], b_r_polar[k]))	
	area_gc_rsplit_intersections[k] = 3600*stuff_intersections[[k]]$area 													 
	area_gc_rsplit[k] = pi*3600*((a_r_polar[k + 1]*b_r_polar[k + 1]) - (a_r_polar[k]*b_r_polar[k]))
	average_radial_distance_gc_rsplit[k] = 60*(a_r_polar[k + 1] + a_r_polar[k])/2
	}
save(stuff_intersections, 
	 file = paste(path_data, galaxy, "stuff_intersections.Rdat", sep = ""))		
idx_area_gc_rsplit_na = which(is.na(area_gc_rsplit))	 
area_gc_rsplit = na.omit(area_gc_rsplit)
average_radial_distance_gc_rsplit = na.omit(average_radial_distance_gc_rsplit)
area_gc_rsplit_intersections = area_gc_rsplit_intersections[-idx_area_gc_rsplit_na]
base_percent_area = round(area_gc_rsplit/area_gc_rsplit_intersections, digits = 1)	
base_percent_area[which(base_percent_area == Inf)] = 0
base_percent_area[which(is.na(base_percent_area))] = 0
percent_area = 2*base_percent_area*seq(1, 2.5*max(base_percent_area), length.out = length(base_percent_area))
## Plotting the diagnostic plot showing the position of the 
## candidate GCs, the elliptical annuli used to estimate the 
## radial density and the grid used to calculate the area of the
## effective elliptical annuli.
plot_annuli(galaxy, 
			path_plot,
			radeg_gc, decdeg_gc, 
			range_ra, range_dec, 
			colorclass_gc, 
			stuff_intersections, 
			num_knots_areas, 
			polygons, polygons_external, polygons_external_avoid,
			ra_length_arrow, dec_length_arrow,
			galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, galaxy1_angle, 
			glob_poly_back_ra, glob_poly_back_dec, 
			"darkgray", 155)  			
## Plot of the not corrected area of the intersections of the elliptical annuli
## with the footprint of the observed region vs the corrected areas of the 
## elliptical annuli.		
pathpdfcomparisonarea = paste(path_plot, "comparison_areas.pdf", sep = "")
pdf(pathpdfcomparisonarea, width = 7, height = 7, paper = "special")
layout(matrix(c(1,2), 2, 1, byrow = TRUE), height = c(3, 1))
# Upper panel.
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
	mai = c(0, 0.5, 0.2, 0.2), mgp = c(1, 0, 0)) 
plot(average_radial_distance_gc_rsplit, 
	 round(area_gc_rsplit_intersections/area_gc_rsplit*100, 1),
	 ylab = "Fraction of area covered by observations", 
	 xlab = "Galactocentric distance", 
	 xlim = range(average_radial_distance_gc_rsplit, na.rm = T),
	 ylim = c(50, 120),
	 pch = 20, 
	 cex = 1, 
	 type = "b")
abline(h = c(100, 90),  
	   col = "black", 
	   lwd = 2)	 
# Other graphical elements.
abline(v = (0.1*10^galaxy1_isodiameter + (0.1*10^galaxy1_isodiameter)/10^(galaxy1_isodiameter_ratio))/4, 
	   col = "black",
	   lty = 2,
	   lwd = 1)	    
# Lower panel.
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
	mai = c(0.5, 0.5, 0, 0.2), mgp = c(1, 0, 0))
cum_area_rsplit_intersections = cumsum(area_gc_rsplit_intersections)
cum_area_rsplit = cumsum(area_gc_rsplit)
fraction_area_gc_rsplit = round(100*cum_area_rsplit_intersections/cum_area_rsplit, 2)
plot(average_radial_distance_gc_rsplit, fraction_area_gc_rsplit, 
	 xlab = "Average major axis of elliptical annuli [arcmin]", 
	 ylab = "% cum. area", 
	 xlim = range(average_radial_distance_gc_rsplit, na.rm = T),
	 ylim = c(50, 120),
	 pch = 20, 
	 cex = 1, 
	 type = "b", 
	 col = color_red_gc)
# Other graphical elements.
abline(v = (0.1*10^galaxy1_isodiameter + (0.1*10^galaxy1_isodiameter)/10^(galaxy1_isodiameter_ratio))/4, 
	   col = "black",
	   lty = 2,
	   lwd = 1)	 
abline(h = seq(0, 100, by = 5), 
	   col = "gray70", 
	   lty = 2, 
	   lwd = 0.5)		   
abline(h = seq(0, 100, by = 10), 
	   col = "gray30", 
	   lty = 1, 
	   lwd = 1)	      	
dev.off()	
#write.table(cbind(average_radial_distance_gc_rsplit, area_gc_rsplit, area_gc_rsplit_intersections), 
#			paste(path_data, galaxy, "/gc_areas_ellipses.dat", sep = ""), 
#			row.names = F, col.names = c("AverageRadius", "AreaAnnuli", "AreaAnnuliFootprint"))	
## Replacing the uncorrected areas with the corrected areas.
#area_gc_rsplit = area_gc_rsplit_intersections
cat("So far so good! Areas of the elliptical annuli corrected for the intersections area calculated!\n")

## Counting the observed candidate GCs in elliptical annuli.
a_r_polar = bins_r_polar
galaxy1_ellipt = 1 - (bins_r_polar[2]/(10^(galaxy1_isodiameter_ratio))/bins_r_polar[2])
# Defining variables.
idx_ellipse_gc_all_rsplit = list()
structure_ellipse_gc_all_rsplit = list()
num_ellipse_gc_all_rsplit = vector(mode = "numeric", length = length(bins_r_polar) - 1)
num_ellipse_gc_all_rsplit_corr = vector(mode = "numeric", length = length(bins_r_polar) - 1)
idx_ellipse_gc_red_rsplit = list()
structure_ellipse_gc_red_rsplit = list()
num_ellipse_gc_red_rsplit = vector(mode = "numeric", length = length(bins_r_polar) - 1)
num_ellipse_gc_red_rsplit_corr = vector(mode = "numeric", length = length(bins_r_polar) - 1)
idx_ellipse_gc_blue_rsplit = list()
structure_ellipse_gc_blue_rsplit = list()
num_ellipse_gc_blue_rsplit = vector(mode = "numeric", length = length(bins_r_polar) - 1)
num_ellipse_gc_blue_rsplit_corr = vector(mode = "numeric", length = length(bins_r_polar) - 1)
## Number of candidate GCs not corrected for contamination.
# All GCs.
structure_ellipse_gc_all_rsplit = points_in_ellipses(a_r_polar, galaxy1_ellipt, galaxy1_angle, 
				   									 radeg_gc, decdeg_gc, 
				   									 galaxy1_center[1], galaxy1_center[2])
idx_ellipse_gc_all_rsplit =	structure_ellipse_gc_all_rsplit$x
num_ellipse_gc_all_rsplit = structure_ellipse_gc_all_rsplit$y										 	   
# Red GCs.
structure_ellipse_gc_red_rsplit = points_in_ellipses(a_r_polar, galaxy1_ellipt, galaxy1_angle, 
				   									 radeg_gc_red, decdeg_gc_red, 
				   									 galaxy1_center[1], galaxy1_center[2])
idx_ellipse_gc_red_rsplit =	structure_ellipse_gc_red_rsplit$x
num_ellipse_gc_red_rsplit = structure_ellipse_gc_red_rsplit$y										 	   
# Blue GCs.
structure_ellipse_gc_blue_rsplit = points_in_ellipses(a_r_polar, galaxy1_ellipt, galaxy1_angle, 
				   									 radeg_gc_blue, decdeg_gc_blue, 
				   									 galaxy1_center[1], galaxy1_center[2])
idx_ellipse_gc_blue_rsplit = structure_ellipse_gc_blue_rsplit$x
num_ellipse_gc_blue_rsplit = structure_ellipse_gc_blue_rsplit$y		
## Number of candidate GCs corrected for contamination from the background region.
num_ellipse_cont_all_rsplit = round(area_gc_rsplit_intersections*dens_cont_all, 0)
num_ellipse_cont_red_rsplit = round(area_gc_rsplit_intersections*dens_cont_red, 0)
num_ellipse_cont_blue_rsplit = round(area_gc_rsplit_intersections*dens_cont_blue, 0)
num_ellipse_gc_all_rsplit_corr = num_ellipse_gc_all_rsplit - num_ellipse_cont_all_rsplit
num_ellipse_gc_red_rsplit_corr = num_ellipse_gc_red_rsplit - num_ellipse_cont_red_rsplit
num_ellipse_gc_blue_rsplit_corr = num_ellipse_gc_blue_rsplit - num_ellipse_cont_blue_rsplit	
num_ellipse_gc_all_rsplit_corr[which(num_ellipse_gc_all_rsplit_corr < 0)] = 0	
num_ellipse_gc_red_rsplit_corr[which(num_ellipse_gc_red_rsplit_corr < 0)] = 0	
num_ellipse_gc_blue_rsplit_corr[which(num_ellipse_gc_blue_rsplit_corr < 0)] = 0	
# Calculating the observed densities of GCs in elliptical annulis (i.e. 
# the number in elliptical annuli divided by the area of each annulus in arcmin2).	
if (sim_flag == "yes") {
	num_gc_all_rsplit = num_ellipse_gc_all_rsplit_corr
	num_gc_red_rsplit = num_ellipse_gc_red_rsplit_corr
	num_gc_blue_rsplit = num_ellipse_gc_blue_rsplit_corr
	num_gc_all_rsplit_sim = num_ellipse_gc_all_rsplit_corr
	num_gc_red_rsplit_sim = num_ellipse_gc_red_rsplit_corr
	num_gc_blue_rsplit_sim = num_ellipse_gc_blue_rsplit_corr
	dens_gc_all_rsplit = num_ellipse_gc_all_rsplit/area_gc_rsplit
	densup_gc_all_rsplit = (num_ellipse_gc_all_rsplit + sqrt(num_ellipse_gc_all_rsplit))/area_gc_rsplit	
	denslow_gc_all_rsplit = (num_ellipse_gc_all_rsplit - sqrt(num_ellipse_gc_all_rsplit))/area_gc_rsplit	
	dens_gc_red_rsplit = num_ellipse_gc_red_rsplit/area_gc_rsplit
	densup_gc_red_rsplit = (num_ellipse_gc_red_rsplit + sqrt(num_ellipse_gc_red_rsplit))/area_gc_rsplit	
	denslow_gc_red_rsplit = (num_ellipse_gc_red_rsplit - sqrt(num_ellipse_gc_red_rsplit))/area_gc_rsplit	
	dens_gc_blue_rsplit = num_ellipse_gc_blue_rsplit/area_gc_rsplit
	densup_gc_blue_rsplit = (num_ellipse_gc_blue_rsplit + sqrt(num_ellipse_gc_blue_rsplit))/area_gc_rsplit	
	denslow_gc_blue_rsplit = (num_ellipse_gc_blue_rsplit - sqrt(num_ellipse_gc_blue_rsplit))/area_gc_rsplit	
	cat("So far so good! Densities of the GCs spatial distribution in elliptical annuli calculated!\n")
	}
	
## Drawing simulated coordinates of sources around the center of the galaxy using
## multiple components, if required.
galaxy1_ellipticity = (galaxy_a - galaxy_b)/galaxy_a
if (sim_flag == "yes") {
	## Creating multiple sets of simulated point distribution of sources.
	# Simulated positions following the elliptical light distribution.
	sim_radeg_all_ellipse_rsplit = list(); sim_radeg_red_ellipse_rsplit = list(); sim_radeg_blue_ellipse_rsplit = list()
	sim_num_radeg_all_ellipse_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_red_ellipse_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_blue_ellipse_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)		
	# Homogeneusly random simulated positions representing the contamination 
	# in the whole covered field and in elliptical annuli.
	sim_radeg_all_cont = list(); sim_radeg_red_cont = list(); sim_radeg_blue_cont = list()
	sim_radeg_all_cont_rsplit = list(); sim_radeg_red_cont_rsplit = list(); sim_radeg_blue_cont_rsplit = list()
	sim_num_radeg_all_cont_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_red_cont_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_blue_cont_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)		
	# Total simulated positions within the elliptical annuli (sum of all the 
	# components considered).
	sim_radec_total_all_rsplit = list(); sim_radec_total_red_rsplit = list(); sim_radec_total_blue_rsplit = list()
	sim_num_radeg_all_total_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_red_total_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)
	sim_num_radeg_blue_total_rsplit = matrix(0, nrow = num_simulations, ncol = length(bins_r_polar_sim) - 1)		
	# Total simulated positions (in the whole footprinf of the observations), sum of all
	# the components considered.
	sim_radec_total_all = list(); sim_radec_total_red = list(); sim_radec_total_blue = list()	
	for (i in seq(num_simulations)) {
		cat(paste("		So far so good! Drawing data for simulation number ", i, "!\n", sep = ""))
		## Elliptical light distribution.
		# All candidate GCs.
		sim_radeg_all_ellipse_rsplit[[i]] = generate_points_ellipse(bins_r_polar_sim, 
														   		 	galaxy1_ellipticity, 
														    		galaxy1_angle, 
														    		sum(num_gc_all_rsplit_sim), 
														    		num_gc_all_rsplit_sim, 
														    		galaxy1_center[1], galaxy1_center[2], 
														    		range_ra, range_dec,
														    		500)
		#random_radeg_sim_all_allrandom[[i]] = random_radeg_sim_all[[i]]$x 
		sim_radeg_all_ellipse_rsplit[[i]] = sim_radeg_all_ellipse_rsplit[[i]]$y		
		# Red candidate GCs
		sim_radeg_red_ellipse_rsplit[[i]] = generate_points_ellipse(bins_r_polar_sim, 
														   		 	galaxy1_ellipticity, 
														    		galaxy1_angle, 
														    		sum(num_gc_red_rsplit_sim), 
														    		num_gc_red_rsplit_sim, 
														    		galaxy1_center[1], galaxy1_center[2], 
														    		range_ra, range_dec,
														    		500)
		#random_radeg_sim_red_allrandom[[i]] = random_radeg_sim_red[[i]]$x 
		sim_radeg_red_ellipse_rsplit[[i]] = sim_radeg_red_ellipse_rsplit[[i]]$y		
		# Blue candidate GCs
		sim_radeg_blue_ellipse_rsplit[[i]] = generate_points_ellipse(bins_r_polar_sim, 
														    		 galaxy1_ellipticity, 
														    		 galaxy1_angle, 
														    		 sum(num_gc_blue_rsplit_sim), 
														    		 num_gc_blue_rsplit_sim, 
														    		 galaxy1_center[1], galaxy1_center[2], 
														    		 range_ra, range_dec,
														    		 500)
		#random_radeg_sim_blue_allrandom[[i]] = random_radeg_sim_blue[[i]]$x 
		sim_radeg_blue_ellipse_rsplit[[i]] = sim_radeg_blue_ellipse_rsplit[[i]]$y
		# Counting simulated sources in radial bins.
		sim_num_radeg_all_ellipse_rsplit[i, ] = unlist(lapply(sim_radeg_all_ellipse_rsplit[[i]], length))/2		
		sim_num_radeg_red_ellipse_rsplit[i, ] = unlist(lapply(sim_radeg_red_ellipse_rsplit[[i]], length))/2
		sim_num_radeg_blue_ellipse_rsplit[i, ] = unlist(lapply(sim_radeg_blue_ellipse_rsplit[[i]], length))/2	
		## Simulating the spatially homogeneous distribution of contaminants.
		# All. 
		sim_radeg_all_cont[[i]] = generate_points_homogeneous(dens_cont_all, 
														   	  range_ra, range_dec)
		sim_radeg_all_cont_rsplit[[i]] = points_in_ellipses(bins_r_polar_sim, 
														    galaxy1_ellipticity, 
														    galaxy1_angle, 
														    unlist(sim_radeg_all_cont[[i]])[, 1], unlist(sim_radeg_all_cont[[i]])[, 2], 	
														    galaxy1_center[1], galaxy1_center[2])$z												    																								   	  
		# Red. 
		sim_radeg_red_cont[[i]] = generate_points_homogeneous(dens_cont_red, 
														   	  range_ra, range_dec)
		sim_radeg_red_cont_rsplit[[i]] = points_in_ellipses(bins_r_polar_sim, 
														    galaxy1_ellipticity, 
														    galaxy1_angle, 
														    unlist(sim_radeg_red_cont[[i]])[, 1], unlist(sim_radeg_red_cont[[i]])[, 2], 	
														    galaxy1_center[1], galaxy1_center[2])$z													    																								   	  
		# Blue. 
		sim_radeg_blue_cont[[i]] = generate_points_homogeneous(dens_cont_blue, 
														   	   range_ra, range_dec)
		sim_radeg_blue_cont_rsplit[[i]] = points_in_ellipses(bins_r_polar_sim, 
														     galaxy1_ellipticity, 
														     galaxy1_angle, 
														     unlist(sim_radeg_blue_cont[[i]])[, 1], unlist(sim_radeg_blue_cont[[i]])[, 2], 	
														     galaxy1_center[1], galaxy1_center[2])$z													    																								   	  
		# Counting simulated sources in radial bins.
		sim_num_radeg_all_cont_rsplit[i, ] = unlist(lapply(sim_radeg_all_cont_rsplit[[i]], length))/2	
		sim_num_radeg_red_cont_rsplit[i, ] = unlist(lapply(sim_radeg_red_cont_rsplit[[i]], length))/2
		sim_num_radeg_blue_cont_rsplit[i, ] = unlist(lapply(sim_radeg_blue_cont_rsplit[[i]], length))/2	
		# Joining the list of simulated positions generated from different components.
		sim_radec_total_all_rsplit[[i]] = mapply(rbind, sim_radeg_all_ellipse_rsplit[[i]], sim_radeg_all_cont_rsplit[[i]])
		sim_radec_total_red_rsplit[[i]] = mapply(rbind, sim_radeg_red_ellipse_rsplit[[i]], sim_radeg_red_cont_rsplit[[i]])
		sim_radec_total_blue_rsplit[[i]] = mapply(rbind, sim_radeg_blue_ellipse_rsplit[[i]], sim_radeg_blue_cont_rsplit[[i]])
		sim_radec_total_all[[i]] = rbind(do.call(rbind, sim_radeg_all_ellipse_rsplit[[i]]), sim_radeg_all_cont[[i]])
		sim_radec_total_red[[i]] = rbind(do.call(rbind, sim_radeg_red_ellipse_rsplit[[i]]), sim_radeg_red_cont[[i]])
		sim_radec_total_blue[[i]] = rbind(do.call(rbind, sim_radeg_blue_ellipse_rsplit[[i]]), sim_radeg_blue_cont[[i]])
		}
	avgsim_num_gc_all_rsplit = colMeans(sim_num_radeg_all_ellipse_rsplit + sim_num_radeg_all_cont_rsplit)
	avgsim_num_gc_red_rsplit = colMeans(sim_num_radeg_red_ellipse_rsplit + sim_num_radeg_red_cont_rsplit)
	avgsim_num_gc_blue_rsplit = colMeans(sim_num_radeg_blue_ellipse_rsplit + sim_num_radeg_blue_cont_rsplit)
	avgsim_dens_gc_all_rsplit = avgsim_num_gc_all_rsplit/area_gc_rsplit
	avgsim_dens_gc_red_rsplit = avgsim_num_gc_red_rsplit/area_gc_rsplit
	avgsim_dens_gc_blue_rsplit = avgsim_num_gc_blue_rsplit/area_gc_rsplit							
	}

## Plot of the distribution of number of background sources (contaminants)
## for each simulation for each radial elliptical annulus.
pathpdf50 = paste(path_plot, "number_background_sources_elliptical_annuli.pdf", sep = "")
pdf(pathpdf50, width = 7, height = 7, paper = "special")
layout(matrix(seq(20), ncol = 4, nrow = 5))
for (y in seq(length(bins_r_polar) - 1)) {
	max_y_histo = max(c(max(hist(sim_num_radeg_all_cont_rsplit[, y], plot = F)$density), 
						max(hist(sim_num_radeg_red_cont_rsplit[, y], plot = F)$density), 
						max(hist(sim_num_radeg_blue_cont_rsplit[, y], plot = F)$density)))	
	par(family = "sans", tcl = 0.5, mai = c(0.3, 0.4, 0.1, 0.1))
	hist(sim_num_radeg_all_cont_rsplit[, y], 
		 freq = F, 
		 xlim = c(0, max(sim_num_radeg_all_cont_rsplit)), 
		 ylim = c(0, max_y_histo),
		 breaks = seq(0, max(sim_num_radeg_all_cont_rsplit), length.out = 20), 
		 border = NA, 
		 col = addalpha("black", 50), 
		 main = "")
	abline(v = area_gc_rsplit[y]*dens_cont_all, 
		  col = "black", 
		  lty = 1, 
		  lwd = 1)
	hist(sim_num_radeg_red_cont_rsplit[, y], 
		 freq = F, 
		 border = NA,
		 col = addalpha(color_red_gc, 50),
		 add = T)
	abline(v = area_gc_rsplit[y]*dens_cont_red, 
		   col = color_red_gc, 
		   lty = 1, 
		   lwd = 1)
	hist(sim_num_radeg_blue_cont_rsplit[, y], 
		 freq = F, 
		 border = NA, 
		 col = addalpha(color_blue_gc, 50), 
		 add = T)
	abline(v = area_gc_rsplit[y]*dens_cont_blue, 
		   col = color_blue_gc, 
		   lty = 1, 
		   lwd = 1)	
	legend("topright", paste(y, "-th annulus"), cex = 0.6, border = NA, pch = NA)   	 
	}
dev.off()
cat("So far so good! Plot of the distribution of random contaminants per elliptical annulus created!\n")

## Plotting the density radial profile of the GCs spatial distribution for
## red and blue sources (no angular slicing).	
# Saving the values of the radial density profiles from GCs to a file
# to be used in the LMXBs program.
save(dens_gc_all_rsplit, densup_gc_all_rsplit, denslow_gc_all_rsplit, num_gc_all_rsplit, 
	 dens_gc_red_rsplit, densup_gc_red_rsplit, denslow_gc_red_rsplit, num_gc_red_rsplit,
	 dens_gc_blue_rsplit, densup_gc_blue_rsplit, denslow_gc_blue_rsplit, num_gc_blue_rsplit,
	 file = paste(path_data, galaxy, "/GCs/radial_density_GCs.RData", 
	  			  sep = ""))
flag_plot_inset = "no"
flag_plot_sim_all = "yes"
flag_plot_sim = "yes"
flag_plot_colors = "yes"
flag_plot_lums = "no"  			  
pathpdf5 = paste(path_plot, "density_radialprofile_ellipses.pdf", sep = "")
pdf(pathpdf5, width = 7, height = 7, paper = "special")
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
	mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
plot(bincenters_r_polar[which(dens_gc_all_rsplit > 0)], 
	 dens_gc_all_rsplit[which(dens_gc_all_rsplit > 0)],
	 xlab = "Radial distance [arcmin]",
	 ylab = expression(paste("Density of sources [", n/arcmin^2, "]", sep = "")), 
	 xlim = c(0, max_radius_fit_arcmin), 
	 ylim = c(0.1, 15),
	 type = "n",
	 log = "y")
# All GCs.
points(bincenters_r_polar[which(dens_gc_all_rsplit > 0)], 
   	   dens_gc_all_rsplit[which(dens_gc_all_rsplit > 0)], 
	   pch = symbol_gc, 
	   col = "black", 
	   cex = cex_gc*2.8, 
	   type = "b") 
segments(bincenters_r_polar[which(dens_gc_all_rsplit > 0)], 
		 densup_gc_all_rsplit[which(dens_gc_all_rsplit > 0)],
		 bincenters_r_polar[which(dens_gc_all_rsplit > 0)],   
		 denslow_gc_all_rsplit[which(dens_gc_all_rsplit > 0)],
		 col = "black", 
		 lwd = lwd_gc, 
	   	 lty = ltype_gc)
# Average simulated densities.
points(bincenters_r_polar[which(avgsim_dens_gc_all_rsplit > 0)], 
	   avgsim_dens_gc_all_rsplit[which(avgsim_dens_gc_all_rsplit > 0)], 
	   pch = 1, 
	   col = "black", 
	   cex = cex_gc*1.8, 
	   type = "b") 		  
# Red GCs.   	 
points(bincenters_r_polar[which(dens_gc_red_rsplit > 0)], 
	   dens_gc_red_rsplit[which(dens_gc_red_rsplit > 0)], 
	   pch = symbol_gc - 2, 
	   col = color_red_gc, 
	   cex = cex_gc*2.6, 
	   type = "b") 
segments(bincenters_r_polar[which(dens_gc_red_rsplit > 0)], 
		 densup_gc_red_rsplit[which(dens_gc_red_rsplit > 0)],
		 bincenters_r_polar[which(dens_gc_red_rsplit > 0)],   
		 denslow_gc_red_rsplit[which(dens_gc_red_rsplit > 0)],
		 col = color_red_gc, 
		 lwd = lwd_gc, 
		 lty = ltype_gc)	   
# Average simulated densities.
points(bincenters_r_polar[which(avgsim_dens_gc_red_rsplit > 0)], 
 	   avgsim_dens_gc_red_rsplit[which(avgsim_dens_gc_red_rsplit > 0)], 
	   pch = 5, 
	   col = color_red_gc, 
	   cex = cex_gc*1.5, 
	   type = "b")   
# Blue GC.
points(bincenters_r_polar[which(dens_gc_blue_rsplit > 0)], 
	   dens_gc_blue_rsplit[which(dens_gc_blue_rsplit > 0)], 
	   pch = symbol_gc - 2, 
	   col = color_blue_gc, 
	   cex = cex_gc*2.6, 
	   type = "b")
segments(bincenters_r_polar[which(dens_gc_blue_rsplit > 0)], 
		 densup_gc_blue_rsplit[which(dens_gc_blue_rsplit > 0)],
		 bincenters_r_polar[which(dens_gc_blue_rsplit > 0)],   
		 denslow_gc_blue_rsplit[which(dens_gc_blue_rsplit > 0)],
		 col = color_blue_gc, 
		 lwd = lwd_gc, 
		 lty = ltype_gc)	   
# Average simulated densities.
points(bincenters_r_polar[which(avgsim_dens_gc_blue_rsplit > 0)], 
	   avgsim_dens_gc_blue_rsplit[which(avgsim_dens_gc_blue_rsplit > 0)], 
	   pch = 5, 
	   col = color_blue_gc, 
	   cex = cex_gc*1.5, 
	   type = "b")	    
## Background level for all, red and blue candidate GCs.
# All
abline(h = dens_cont_all, col = "black", lty = 2, lwd = 3) 
# Red
abline(h = dens_cont_red, col = color_red_gc, lty = 2, lwd = 3) 
# Blue
abline(h = dens_cont_blue, col = color_blue_gc, lty = 2, lwd = 3) 
# Other graphical elements.
abline(v = (0.1*10^galaxy1_isodiameter + (0.1*10^galaxy1_isodiameter)/10^(galaxy1_isodiameter_ratio))/4, 
	   col = "black",
	   lty = 2,
	   lwd = 1)
# Customizing the legend of the plot.
if (flag_plot_colors == "yes" & flag_plot_lums == "yes") {
#	legend_legend = c("All GCs", "Red GCs", "Blue GCs", "High L", "Mid L", "Low L")
#	col_legend = c("black", color_red_gc, "blue", "black", "black", "black")
#	pch_legend = c(rep(symbol_gc, 3), 17, 6, 2)
	legend_legend = c("All GCs", "Red GCs", "Blue GCs", "High L", "Low L")
	col_legend = c("black", color_red_gc, color_blue_gc, "black", "black")
	pch_legend = c(20, 18, 18, 17, 15)
	} else if (flag_plot_colors == "no" & flag_plot_lums == "yes") {
#		legend_legend = c("All GCs", "High L", "Mid L", "Low L")
#		col_legend = c("black", "black", "black")
#		pch_legend = c(symbol_gc, 17, 6, 2)		
		legend_legend = c("All GCs", "High L", "Low L")
		col_legend = c("black", "black")
		pch_legend = c(symbol_gc, 17, 2)		
		} else if (flag_plot_colors == "yes" & flag_plot_lums == "no") {
			legend_legend = c("All GCs", "Red GCs", "Blue GCs")
			col_legend = c("black", color_red_gc, color_blue_gc)
			pch_legend = c(symbol_gc, symbol_gc - 2, symbol_gc - 2)						
			}			     
legend("topright", 
	   legend_legend,
	   col = col_legend,
	   pch = pch_legend, 
	   cex = size_legend,
	   pt.cex = 2,
	   inset = 0.03, 
	   bty = "n")
if (flag_plot_inset == "yes") {	   
	par(fig = c(0.05, 0.5, 0.05, 0.5), new = TRUE)	   
	plot_annuli(radeg_gc, 
				decdeg_gc)	   
	}		
ntx = par()$xaxp
ntx[3] = ntx[3]*5
nty = par()$yaxp
nty[3] = nty[3]*2
par(xaxp = ntx, tcl = 0.2)
axis(1, label = FALSE)
axis(2, label = FALSE)
axis(3, label = FALSE)
axis(4, label = FALSE)
dev.off()	
cat("So far so good! Plot of the radial surface density profiles of the field population and galaxy population created!\n")

## KNN density estimation with different values of N for the distribution of all, blue and
## red GCs, high-luminosity and low-luminosity GCs.
# Determining the grid for density evaluation.
radeg_gridpoint_all = seq(from = range(radeg_gc)[1], to = range(radeg_gc)[2], 
						  length.out = nsteps_smooth)
decdeg_gridpoint_all = seq(from = range(decdeg_gc)[1], to = range(decdeg_gc)[2], 
						   length.out = nsteps_smooth)
grid_dens = expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)
# Size of the pixel along R.A. and Dec. axes, in degrees and 
# arcsecs.
cat(paste("So far so good! Pixel size along R.A. axis [deg]: ", round(radeg_gridpoint_all[2] - radeg_gridpoint_all[1], 2), "\n", sep  = ""))
cat(paste("So far so good! Pixel size along Dec. axis [deg]: ", round(decdeg_gridpoint_all[2] - decdeg_gridpoint_all[1], 2), "\n", sep  = ""))
cat(paste("So far so good! Pixel size along R.A. axis [arcsec]: ", round(3600*(radeg_gridpoint_all[2] - radeg_gridpoint_all[1]), 2), "\n", sep  = ""))
cat(paste("So far so good! Pixel size along Dec. axis [arcsec]: ", round(3600*(decdeg_gridpoint_all[2] - decdeg_gridpoint_all[1]), 2), "\n", sep  = ""))
				  
# Area of each pixel used throughout the program (in arcmin squared).
area_pixel = 3600*((max(radeg_gc) - min(radeg_gc))*(max(decdeg_gc) - min(decdeg_gc))/(nsteps_smooth^2))
width_gridpoint_radeg = radeg_gridpoint_all[2] - radeg_gridpoint_all[1]
width_gridpoint_decdeg = decdeg_gridpoint_all[2] - decdeg_gridpoint_all[1]
## Determining the distances of each grid center from all the points of the 
## distribution.
dist_all = list()
dist_all_rcolor = list()
dist_blue = list(); dist_red = list()
dens_point_all = list()
dens_point_blue = list(); dens_point_red = list()
errdens_point_all = list()
errdens_point_blue = list(); errdens_point_red = list()
names_dens_all = list()
names_dens_red = list(); names_dens_blue = list()
# Average color for cell.
color_all = vector(mode = "numeric", length(grid_dens[, 1]))
numgcs_all = vector(mode = "numeric", length(grid_dens[, 1]))
numgcs_red = vector(mode = "numeric", length(grid_dens[, 1]))
numgcs_blue = vector(mode = "numeric", length(grid_dens[, 1]))
step_radeg = radeg_gridpoint_all[2] - radeg_gridpoint_all[1]
step_decdeg = decdeg_gridpoint_all[2] - decdeg_gridpoint_all[1]
for (i in 1:length(grid_dens[, 1])) {
	# Distances of all points from the central point of each cell used
	# to evaluate the density.
	numgcs_all[i] = length(which(abs(radeg_gc - grid_dens[i, 1]) <= step_radeg/2 & 
							     abs(decdeg_gc - grid_dens[i, 2]) <= step_decdeg/2))
	numgcs_red[i] = length(which(abs(radeg_gc[which(gmi_gc >= color_threshold)] - grid_dens[i, 1]) <= step_radeg/2 & 
							     abs(decdeg_gc[which(gmi_gc >= color_threshold)] - grid_dens[i, 2]) <= step_decdeg/2))
	numgcs_blue[i] = length(which(abs(radeg_gc[which(gmi_gc < color_threshold)] - grid_dens[i, 1]) <= step_radeg/2 & 
					 		      abs(decdeg_gc[which(gmi_gc < color_threshold)] - grid_dens[i, 2]) <= step_decdeg/2))
	color_all[i] = mean(gmi_gc[which(abs(radeg_gc - grid_dens[i, 1]) <= step_radeg & 
									 abs(decdeg_gc - grid_dens[i, 2]) <= step_decdeg)])
	dist_all[[i]] = sqrt((grid_dens[i, 1] - radeg_gc)^2 + (grid_dens[i, 2] - decdeg_gc)^2)
	dist_blue[[i]] = sqrt((grid_dens[i, 1] - radeg_gc_colorclass[[1]])^2 + 
						  (grid_dens[i, 2] - decdeg_gc_colorclass[[1]])^2)
	dist_red[[i]] = sqrt((grid_dens[i, 1] - radeg_gc_colorclass[[2]])^2 + 
						 (grid_dens[i, 2] - decdeg_gc_colorclass[[2]])^2)
	## Names of the sources contained in each cell used to evaluate the 
	## density.
	names_dens_all[[i]] = as.vector(id_gc[which(radeg_gc >= grid_dens[i, 1] & radeg_gc < (grid_dens[i, 1] + width_gridpoint_radeg) &
									  		    decdeg_gc >= grid_dens[i, 2] & decdeg_gc < (grid_dens[i, 2] + width_gridpoint_decdeg))])
	names_dens_red[[i]] = as.vector(id_gc_red[which(radeg_gc_red >= grid_dens[i, 1] & radeg_gc_red < (grid_dens[i, 1] + width_gridpoint_radeg) &
									  		    	decdeg_gc_red >= grid_dens[i, 2] & decdeg_gc_red < (grid_dens[i, 2] + width_gridpoint_decdeg))])
	names_dens_blue[[i]] = as.vector(id_gc_blue[which(radeg_gc_blue >= grid_dens[i, 1] & radeg_gc_blue < (grid_dens[i, 1] + width_gridpoint_radeg) &
									  		    	  decdeg_gc_blue >= grid_dens[i, 2] & decdeg_gc_blue < (grid_dens[i, 2] + width_gridpoint_decdeg))])
	## Density estimation for each point of the grid for each different population.
	dens_point_all[[i]] = vector(mode = "numeric", length = length(k_vector))
	dens_point_blue[[i]] = vector(mode = "numeric", length = length(k_vector))
	dens_point_red[[i]] = vector(mode = "numeric", length = length(k_vector))
	errdens_point_all[[i]] = vector(mode = "numeric", length = length(k_vector))
	errdens_point_blue[[i]] = vector(mode = "numeric", length = length(k_vector))
	errdens_point_red[[i]] = vector(mode = "numeric", length = length(k_vector))
	for (k in 1:length(k_vector)) {
		dens_point_all[[i]][k] = k_vector[k]/(pi*(sort(dist_all[[i]])[k_vector[k]])^2)
		dens_point_blue[[i]][k] = k_vector[k]/(pi*(sort(dist_blue[[i]])[k_vector[k]])^2)
		dens_point_red[[i]][k] = k_vector[k]/(pi*(sort(dist_red[[i]])[k_vector[k]])^2)
		errdens_point_all[[i]][k] = sqrt(k_vector[k])/(pi*(sort(dist_all[[i]])[k_vector[k]])^2)
		errdens_point_blue[[i]][k] = sqrt(k_vector[k])/(pi*(sort(dist_blue[[i]])[k_vector[k]])^2)
		errdens_point_red[[i]][k] = sqrt(k_vector[k])/(pi*(sort(dist_red[[i]])[k_vector[k]])^2)
		}
	}
cat("So far so good! KNN spatial densities for different values of K calculated!\n")

## Determining the pixels within field observed by HST. 
idx_within_footprint = gridknots_inside_polygon(radeg_gridpoint_all, 
											   decdeg_gridpoint_all, 
											   glob_poly_xy)											   
idx_within_footprint = idx_within_footprint
cat("So far so good! Pixels within the HST field determined!\n")

## Plot of the observed KNN densities for the different distibutions of 
## sources (all, red and blue GCs).
idx_densthresh_gc_all = list()
idx_densthresh_gc_red = list(); idx_densthresh_gc_blue = list()
percentile_2d_density_all = vector(mode = "numeric", length(k_vector)) 
percentile_2d_density_red = vector(mode = "numeric", length(k_vector)) 
percentile_2d_density_blue = vector(mode = "numeric", length(k_vector)) 
percentile_2d_normdensity_all = vector(mode = "numeric", length(k_vector)) 
percentile_2d_normdensity_red = vector(mode = "numeric", length(k_vector)) 
percentile_2d_normdensity_blue = vector(mode = "numeric", length(k_vector)) 
for (k in 1:length(k_vector)) {
	# Normalizing densities.
	dens_gc_all = sapply(dens_point_all, "[", k)
	dens_gc_blue = sapply(dens_point_blue, "[", k)
	dens_gc_red = sapply(dens_point_red, "[", k)	
	normdens_gc_all = dens_gc_all/max(dens_gc_all)
	normdens_gc_blue = dens_gc_blue/max(dens_gc_blue)
	normdens_gc_red = dens_gc_red/max(dens_gc_red)	
	percentile_2d_density_all[k] = quantile(dens_gc_all, pixel_percentile)
	percentile_2d_density_red[k] = quantile(dens_gc_red, pixel_percentile)
	percentile_2d_density_blue[k] = quantile(dens_gc_blue, pixel_percentile)	
	percentile_2d_normdensity_all[k] = quantile(normdens_gc_all, pixel_percentile)
	percentile_2d_normdensity_red[k] = quantile(normdens_gc_red, pixel_percentile)
	percentile_2d_normdensity_blue[k] = quantile(normdens_gc_blue, pixel_percentile)				
	# All GCs.
	radeg_normdens_gc_all = radeg_gridpoint_all; decdeg_normdens_gc_all = decdeg_gridpoint_all
	dens_threshold_gc_all = lseq(min(normdens_gc_all), max(normdens_gc_all), length.out = 10)[index_dens_threshold]
	idx_densthresh_gc_all[[k]] = which(normdens_gc_all >= dens_threshold_gc_all)
	matrix_normdens_gc_all = matrix(dens_gc_all, nsteps_smooth, nsteps_smooth)
	matrix_normdens_gc_all = matrix(normdens_gc_all, nsteps_smooth, nsteps_smooth)		
	save(radeg_normdens_gc_all, decdeg_normdens_gc_all, matrix_normdens_gc_all, 
		 file = paste(path_data, galaxy, "/GCs/dens2d_knn_all_K_", k_vector[k], ".RData", 
		 			  sep = ""))
	pathpdfdensknnall = paste(path_plot, "dens2d_all_knn_K_", k_vector[k], ".pdf", sep = "")	
	pdf(pathpdfdensknnall, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
	     xlim = rev(range_ra), 
     	 ylim = range_dec, 
		 main = "", 
		 type = "n", 
		 axes = F)			
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, z = matrix_normdens_gc_all,
          col = rev(heat.colors(100)),
		  add = TRUE)	
	points(radeg_gc, decdeg_gc, 
		   pch = 20,
		   cex = 0.3,
		   col = "darkgray")
	contour(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
			z = matrix_normdens_gc_all,
			add = TRUE, 
			levels = lseq(min(matrix_normdens_gc_all), max(matrix_normdens_gc_all), length.out = 10)[8:10], 
			drawlabels = FALSE, 
			col = "black", 
			lwd = 1)		  
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))			   
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)			
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	## Labels on the plot.			  
	#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[k], 0)	
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)	
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)			
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
	# Blue GCs.
	dens_threshold_gc_blue = lseq(min(normdens_gc_blue), max(normdens_gc_blue), length.out = 10)[index_dens_threshold]
	idx_densthresh_gc_blue[[k]] = which(normdens_gc_blue >= dens_threshold_gc_blue)
	matrix_normdens_gc_blue = matrix(normdens_gc_blue, nsteps_smooth, nsteps_smooth)
	save(radeg_normdens_gc_all, decdeg_normdens_gc_all, matrix_normdens_gc_blue, 
		 file = paste(path_data, galaxy, "/GCs/dens2d_knn_blue_K_", k_vector[k], ".RData", 
		 			  sep = ""))
	pathpdfdensknnblue = paste(path_plot, "dens2d_blue_knn_K_", k_vector[k], ".pdf", sep = "")	
	pdf(pathpdfdensknnblue, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
	     xlim = rev(range_ra),
     	 ylim = range_dec,  
		 main = "", 
		 type = "n", 
		 axes = F)
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, z = matrix_normdens_gc_blue,
          col = rev(heat.colors(100)), 
     	  add = TRUE)	
	points(radeg_gc_colorclass[[1]], decdeg_gc_colorclass[[1]], 
		   pch = 20,
		   cex = 0.3, 
		   col = "darkgray")
	contour(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
			z = matrix_normdens_gc_blue,
			add = TRUE, 
			levels = lseq(min(matrix_normdens_gc_blue), max(matrix_normdens_gc_blue), length.out = 10)[8:10], 
			drawlabels = FALSE, 
			col = "black", 
			lwd = 1) 
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)					 	   
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	## Labels on the plot.			  
	#plot_labels("topleft", "Blue GCs", 0.05, 1.4, k_vector[k], 0)		
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)		
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)					  
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
	# Red GCs.
	dens_threshold_gc_red = lseq(min(normdens_gc_red), max(normdens_gc_red), length.out = 10)[index_dens_threshold]
	idx_densthresh_gc_red[[k]] = which(normdens_gc_red >= dens_threshold_gc_red)
	matrix_normdens_gc_red = matrix(normdens_gc_red, nsteps_smooth, nsteps_smooth)
	save(radeg_normdens_gc_all, decdeg_normdens_gc_all, matrix_normdens_gc_red, 
		 file = paste(path_data, galaxy, "/GCs/dens2d_knn_red_K_", k_vector[k], ".RData", 
		 			  sep = ""))
	pathpdfdensknnred = paste(path_plot, "dens2d_red_knn_K_", k_vector[k], ".pdf", sep = "")	
	pdf(pathpdfdensknnred, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
	     xlim = rev(range_ra), 
     	 ylim = range_dec,  
		 main = "", 
		 type = "n", 
		 axes = F)		
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, z = matrix_normdens_gc_red,
          col = rev(heat.colors(100)), 
     	  add = TRUE)	
	points(radeg_gc_colorclass[[2]], decdeg_gc_colorclass[[2]], 
		   pch = 20,
		   cex = 0.3, 
		   col = "darkgray")
	contour(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
			z = matrix_normdens_gc_red,
			add = TRUE, 
			levels = lseq(min(matrix_normdens_gc_red), max(matrix_normdens_gc_red), length.out = 10)[8:10], 
			drawlabels = FALSE, 
			col = "black", 
			lwd = 1) 	   
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)	
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  				   
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)	
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)			
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)	
	## Labels on the plot.			  
	#plot_labels("topleft", "Red GCs", 0.05, 1.4, k_vector[k], 0)			
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
cat("So far so good! KNN spatial densities maps from observed distribution of GCs created!\n")

if (sim_flag == "yes") {
	## Deriving the densities from the simulated random distribution of all, red, blue, 
	## low-L and high-L sample of GCs. 
	sim_dens_point_all = list()
	sim_dens_point_red = list(); sim_dens_point_blue = list()
	sim_errdens_point_all = list() 
	sim_errdens_point_red = list(); sim_errdens_point_blue = list()
	sim_radial_ellipse_all = list() 
	sim_radial_ellipse_red = list(); sim_radial_ellipse_blue = list()
	sim_num_pixels_above_threshold_all = matrix(0, nrow = num_simulations, ncol = length(k_vector))
	sim_num_pixels_above_threshold_red = matrix(0, nrow = num_simulations, ncol = length(k_vector))
	sim_num_pixels_above_threshold_blue = matrix(0, nrow = num_simulations, ncol = length(k_vector))
	sim_idx_pixels_densthresh_all = list()
	sim_idx_pixels_densthresh_red = list()
	sim_idx_pixels_densthresh_blue = list()
	sim_num_pixels_inside_contour_all = matrix(0, nrow = num_simulations, ncol = length(k_vector))
	sim_num_pixels_inside_contour_red = matrix(0, nrow = num_simulations, ncol = length(k_vector)) 
	sim_num_pixels_inside_contour_blue = matrix(0, nrow = num_simulations, ncol = length(k_vector))
	## Using the simulated coordinates already evaluated to create the 
	## scatterplots of the simulated distributions of GCs for all the 
	## different classes of GCs (all, red, blue, high-L, low-L).
	for (i in seq(num_simulations)) {
		# Plotting the simulated spatial distributions of GCs.
		# All GCs.
		pathpdfsimulatedradec_all = paste(path_plot_simulations, "positions_all_simulation_", i, ".pdf", sep = "")	
		pdf(pathpdfsimulatedradec_all, width = 7, height = 7, paper = "special")
		par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
			mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
		plot(sim_radec_total_all[[i]][, 1], 
			 sim_radec_total_all[[i]][, 2], 
			 xlab = "Right Ascension [deg]", 
			 ylab = "Declination [deg]",
			 xlim = rev(range_ra), 
			 ylim = range_dec,
			 type = "n", 
			 axes = F)
		points(sim_radec_total_all[[i]][, 1], 
			   sim_radec_total_all[[i]][, 2], 
			   pch = symbol_gc, 
			   cex = cex_gc, 
			   col = "black")
		# Drawing the polygons of the HST fields. 
		plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   	   "black", 0)	
		# Shape of the first galaxy.
		plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, galaxy1_angle, cex_center = cex_center_galaxy)
		# Embellishing the axes.	
		axis(1, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.025), label = FALSE, tcl = 0.2)
		axis(2, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.025), label = FALSE, tcl = 0.2)
		axis(3, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.2), label = FALSE, tcl = 0.5)
		axis(3, at = seq(min(range_ra) - 1, max(range_ra) + 1, by = 0.025), label = FALSE, tcl = 0.2)
		axis(4, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.2), label = FALSE, tcl = 0.5)
		axis(4, at = seq(min(range_dec) - 1, max(range_dec) + 1, by = 0.025), label = FALSE, tcl = 0.2)
		dev.off()							 
		# Red GCs.
		pathpdfsimulatedradec_red = paste(path_plot_simulations, "positions_red_simulation_", i, ".pdf", sep = "")	
		pdf(pathpdfsimulatedradec_red, width = 7, height = 7, paper = "special")
		par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
			mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
		plot(sim_radec_total_red[[i]][, 1], 
			 sim_radec_total_red[[i]][, 2], 
			 xlab = "Right Ascension [deg]", 
			 ylab = "Declination [deg]",
			 xlim = rev(range_ra), 
			 ylim = range_dec,
			 type = "n", 
		 	 axes = F)
		points(sim_radec_total_red[[i]][, 1], 
			   sim_radec_total_red[[i]][, 2], 
			   pch = symbol_gc, 
			   cex = cex_gc,			
			   col = color_red_gc)
		# Drawing the polygons of the HST fields. 
		plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   	   "black", 0)
		# Shape of the first galaxy.
		plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, galaxy1_angle, cex_center = cex_center_galaxy)
		# Plotting the angular scale.
		plot_scale("topright", 180, 0)		
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
		# Blue GCs.
		pathpdfsimulatedradec_blue = paste(path_plot_simulations, "positions_blue_simulation_", i, ".pdf", sep = "")	
		pdf(pathpdfsimulatedradec_blue, width = 7, height = 7, paper = "special")
		par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
			mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
		plot(sim_radec_total_blue[[i]][, 1], 
			 sim_radec_total_blue[[i]][, 2], 
			 xlab = "Right Ascension [deg]", 
			 ylab = "Declination [deg]",
			 xlim = rev(range_ra), 
			 ylim = range_dec,
			 type = "n", 
		 	 axes = F)
		points(sim_radec_total_blue[[i]][, 1], 
			   sim_radec_total_blue[[i]][, 2], 
			   pch = symbol_gc, 
			   cex = cex_gc,
			   col = color_blue_gc)
		# Drawing the polygons of the HST fields. 
		plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   	   "black", 0)			   
		# Shape of the first galaxy.
		plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, galaxy1_angle, cex_center = cex_center_galaxy)
		# Plotting the angular scale.
		plot_scale("topright", 180, 0)		
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
		# Evaluating the density with the KNN method on the simulated 
		# distribution of all, blue and red GCs.
		sim_dist_all = list() 
		sim_dist_red = list() 
		sim_dist_blue = list()
		sim_dens_point_all[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		sim_dens_point_blue[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		sim_dens_point_red[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		sim_errdens_point_all[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		sim_errdens_point_blue[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		sim_errdens_point_red[[i]] = matrix(0, nrow = length(grid_dens[, 1]), ncol = length(k_vector))
		for (k in 1:length(grid_dens[, 1])) {
			sim_radeg_all = sim_radec_total_all[[i]][, 1]
			sim_decdeg_all = sim_radec_total_all[[i]][, 2]
			sim_radeg_red = sim_radec_total_red[[i]][, 1]
			sim_decdeg_red = sim_radec_total_red[[i]][, 2]
			sim_radeg_blue = sim_radec_total_blue[[i]][, 1]
			sim_decdeg_blue = sim_radec_total_blue[[i]][, 2]
			sim_dist_all[[k]] = sqrt((grid_dens[k, 1] - sim_radeg_all)^2 + 
									 (grid_dens[k, 2] - sim_decdeg_all)^2)
			sim_dist_blue[[k]] = sqrt((grid_dens[k, 1] - sim_radeg_blue)^2 + 
									  (grid_dens[k, 2] - sim_decdeg_blue)^2)
			sim_dist_red[[k]] = sqrt((grid_dens[k, 1] - sim_radeg_red)^2 + 
									 (grid_dens[k, 2] - sim_decdeg_red)^2)
			# Density estimation for each point of the grid for each different population.
			for (f in 1:length(k_vector)) {
				sim_dens_point_all[[i]][k, f] = k_vector[f]/(pi*(sort(sim_dist_all[[k]])[k_vector[f]])^2)
				sim_dens_point_blue[[i]][k, f] = k_vector[f]/(pi*(sort(sim_dist_blue[[k]])[k_vector[f]])^2)
				sim_dens_point_red[[i]][k, f] = k_vector[f]/(pi*(sort(sim_dist_red[[k]])[k_vector[f]])^2)
				sim_errdens_point_all[[i]][k, f] = sqrt(k_vector[f])/(pi*(sort(sim_dist_all[[k]])[k_vector[f]])^2)
				sim_errdens_point_blue[[i]][k, f] = sqrt(k_vector[f])/(pi*(sort(sim_dist_blue[[k]])[k_vector[f]])^2)
				sim_errdens_point_red[[i]][k, f] = sqrt(k_vector[f])/(pi*(sort(sim_dist_red[[k]])[k_vector[f]])^2)			
				}
			}
		cat(paste("		So far so good! KNN densities of the ", i, "-th simulated coordinates of all GCs classes calculated!\n", sep = ""))
		}
	cat("So far so good! Simulations of the positions of all, red and blue, high-L and low-L GCs created!\n")							   						   							   
	
	## Plotting the density from the simulated spatial distributions and 
	## calculating the residual for each pixel of the 2D densities.
	# Defining the observed spatial distribution for the three classes
	# of the sources.
	res_all = list()
	res_red = list(); res_blue = list()
	errres_all = list()
	errres_red = list(); errres_blue = list()
	num_sigmas_all = list()
	num_sigmas_red = list(); num_sigmas_blue = list()
	perc_tail_all = list()
	perc_tail_red = list(); perc_tail_blue = list()
	num_more1sigmas_simulation_all = list()
	num_more1sigmas_simulation_blue = list()
	num_more1sigmas_simulation_red = list()
	num_more2sigmas_simulation_all = list()
	num_more2sigmas_simulation_blue = list()
	num_more2sigmas_simulation_red = list()
	num_more3sigmas_simulation_all = list()
	num_more3sigmas_simulation_blue = list()
	num_more3sigmas_simulation_red = list()		
	perc_more1sigmas_all = vector(mode = "numeric", length = length(k_vector))
	perc_more1sigmas_red = vector(mode = "numeric", length = length(k_vector))
	perc_more1sigmas_blue = vector(mode = "numeric", length = length(k_vector))		
	for (f in 1:length(k_vector)) {
		res_all[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))
		res_red[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth)) 
		res_blue[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth)) 
		errres_all[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))
		errres_red[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth)) 
		errres_blue[[f]] = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth)) 
		num_sigmas_all[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		num_sigmas_red[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		num_sigmas_blue[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		perc_tail_all[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		perc_tail_red[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		perc_tail_blue[[f]] = array(0, dim = c(nsteps_smooth, nsteps_smooth))
		num_sigmas_simulation_all = list()
		num_sigmas_simulation_blue = list()
		num_sigmas_simulation_red = list()
		num_more1sigmas_simulation_all[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more1sigmas_simulation_red[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more1sigmas_simulation_blue[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more2sigmas_simulation_all[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more2sigmas_simulation_red[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more2sigmas_simulation_blue[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more3sigmas_simulation_all[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more3sigmas_simulation_red[[f]] = vector(mode = "numeric", length = num_simulations)
		num_more3sigmas_simulation_blue[[f]] = vector(mode = "numeric", length = num_simulations)
		for (i in seq(num_simulations)) {
			# Deriving the normalized observed densities of all, red and blue, high-L and low-L
			# GCs.
			dens_all = sapply(dens_point_all, "[", f)
			dens_blue = sapply(dens_point_blue, "[", f)
			dens_red = sapply(dens_point_red, "[", f)
			#normdens_all = dens_all/sum(dens_all)
			#normdens_blue = dens_blue/sum(dens_blue)
			#normdens_red = dens_red/sum(dens_red)
			normdens_all = dens_all/max(dens_all)
			normdens_blue = dens_blue/max(dens_blue)
			normdens_red = dens_red/max(dens_red)
			errdens_all = sapply(errdens_point_all, "[", f)
			errdens_blue = sapply(errdens_point_blue, "[", f)
			errdens_red = sapply(errdens_point_red, "[", f)	
			#errnormdens_all = errdens_all/sum(dens_all)
			#errnormdens_red = errdens_red/sum(dens_red)
			#errnormdens_blue = errdens_blue/sum(dens_blue)
			errnormdens_all = errdens_all/max(dens_all)
			errnormdens_red = errdens_red/max(dens_red)
			errnormdens_blue = errdens_blue/max(dens_blue)
			# Deriving the normalized simulated densities of all, red and blue
			# GCs for the i-th simulations.
			# Normalizing densities.
			sim_dens_all = sim_dens_point_all[[i]][, f]
			sim_dens_blue = sim_dens_point_blue[[i]][, f]
			sim_dens_red = sim_dens_point_red[[i]][, f]
			sim_errdens_all = sim_errdens_point_all[[i]][, f]
			sim_errdens_blue = sim_errdens_point_blue[[i]][, f]
			sim_errdens_red = sim_errdens_point_red[[i]][, f]
			sim_normdens_all = sim_dens_all/max(sim_dens_all)
			sim_normdens_blue = sim_dens_blue/max(sim_dens_blue)
			sim_normdens_red = sim_dens_red/max(sim_dens_red)
			sim_errnormdens_all = sim_errdens_all/max(sim_dens_all)
			sim_errnormdens_blue = sim_errdens_blue/max(sim_dens_blue)
			sim_errnormdens_red = sim_errdens_red/max(sim_dens_red)
			# Counting the number of pixels above the threshold derived from the 
			# percentile and within the density contour derived from the 
			# 2D density maps.
			sim_idx_pixels_densthresh_all[[i]] = which(sim_dens_all >= percentile_2d_density_all[f])
			sim_idx_pixels_densthresh_red[[i]] = which(sim_dens_red >= percentile_2d_density_red[f])
			sim_idx_pixels_densthresh_blue[[i]] = which(sim_dens_blue >= percentile_2d_density_blue[f])
			sim_num_pixels_above_threshold_all[i, f] = length(sim_idx_pixels_densthresh_all[[i]])
			sim_num_pixels_above_threshold_red[i, f] = length(sim_idx_pixels_densthresh_red[[i]])
			sim_num_pixels_above_threshold_blue[i, f] = length(sim_idx_pixels_densthresh_blue[[i]])
			sim_num_pixels_inside_contour_all[i, f] = length(which(sim_idx_pixels_densthresh_all[[i]] %in% idx_densthresh_gc_all[[f]]))
			sim_num_pixels_inside_contour_red[i, f] = length(which(sim_idx_pixels_densthresh_red[[i]] %in% idx_densthresh_gc_red[[f]]))
			sim_num_pixels_inside_contour_blue[i, f] = length(which(sim_idx_pixels_densthresh_blue[[i]] %in% idx_densthresh_gc_blue[[f]]))						
			cat(paste("		So far so good! Number of pixels above the threshold for the ", i, "-th simulated coordinates obtained with the ", 
					  f, "-th K values of all GC classes counted!\n", sep = ""))
			# Calculating the residual map for all, red and blue
			# GCs for the i-th simulation and f-th value of K.
			res_all[[f]][i, , ] = normdens_all - sim_normdens_all
			res_red[[f]][i, , ] = normdens_red - sim_normdens_red
			res_blue[[f]][i, , ] = normdens_blue - sim_normdens_blue
			errres_all[[f]][i, , ] = sqrt(errnormdens_all^2 + sim_errnormdens_all^2)
			errres_red[[f]][i, , ] = sqrt(errnormdens_red^2 + sim_errnormdens_red^2)
			errres_blue[[f]][i, , ] = sqrt(errnormdens_blue^2 + sim_errnormdens_blue^2)
			## Plots of the simulated spatial distributions.
			# All (simulated) GCs.
			pathpdfdensknnall_sim = paste(path_plot_simulations, "sim_dens2d_knn_all_K_", 
										  k_vector[f], "_simulation_", i,".pdf", sep = "")	
			pdf(pathpdfdensknnall_sim, width = 7, height = 7, paper = "special")
			par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
				mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
			sim_matrix_normdens_gc_all = matrix(sim_normdens_all, nsteps_smooth, nsteps_smooth)
			plot(rev(range_ra), range_dec,
				 xlim = rev(range_ra), 
     			 ylim = range_dec,  
 			     xlab = "Right Ascension [deg]", 
				 ylab = "Declination [deg]",
 				 type = "n",
 				 main = paste("2D density, KNN, K=", k_vector[f], 
				 			  ", all (simulated) GCs, i=", i,  sep = ""), 
		 		axes = F)
			image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
				  z = sim_matrix_normdens_gc_all,
				  col = rev(heat.colors(100)), 
				  add = T)
		    points(sim_radec_total_all[[i]][, 1], sim_radec_total_all[[i]][, 2], 
		   		   pch = 20,
		   		   cex = 0.3,
		   		   col = "darkgray")				  
			# Drawing the polygons of the HST fields. 
			plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   		   "black", 0)		
			# Drawing arrows indicating the North and West directions. 	 
			plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)		   		   						   
			# Shape of the first galaxy.
			plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
						  galaxy1_angle, cex_center = cex_center_galaxy)
			## Labels on the plot.			  
			#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[f], 0)	
			# Plotting the angular scale.
			plot_scale("topright", 180, 0)							  
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
			# Red (simulated) GCs.																	
			pathpdfdensknnred_sim = paste(path_plot_simulations, "sim_dens2d_knn_red_K_", 
										  k_vector[f], "_simulation_", i,".pdf", sep = "")	
			pdf(pathpdfdensknnred_sim, width = 7, height = 7, paper = "special")
			par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
				mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
			sim_matrix_normdens_gc_red = matrix(sim_normdens_red, nsteps_smooth, nsteps_smooth)
			plot(rev(range_ra), range_dec,
				 xlim = rev(range_ra), 
     			 ylim = range_dec,  
 			     xlab = "Right Ascension [deg]", 
				 ylab = "Declination [deg]",
 				 type = "n",
 				 main = paste("2D density, KNN, K=", k_vector[f], 
				 			  ", all (simulated) GCs, i=", i,  sep = ""), 
		 		 axes = F)
			image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
				  z = sim_matrix_normdens_gc_red,
				  col = rev(heat.colors(100)), 
				  add = T)
		    points(sim_radec_total_red[[i]][, 1], sim_radec_total_red[[i]][, 2], 
		   		   pch = 20,
		   		   cex = 0.3,
		   		   col = "darkgray")
			# Drawing the polygons of the HST fields. 
			plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   		   "black", 0)							   
			# Drawing arrows indicating the North and West directions. 	 
			plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)		   		   						   
			# Shape of the first galaxy.
			plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
						  galaxy1_angle, cex_center = cex_center_galaxy)
			## Labels on the plot.			  
			#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[f], 0)	
			# Plotting the angular scale.
			plot_scale("topright", 180, 0)	
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
			# Blue (simulated) GCs.																	
			pathpdfdensknnblue_sim = paste(path_plot_simulations, "sim_dens2d_knn_blue_K_", 
										   k_vector[f], "_simulation_", i,".pdf", sep = "")	
			pdf(pathpdfdensknnblue_sim, width = 7, height = 7, paper = "special")
			par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
				mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
			sim_matrix_normdens_gc_blue = matrix(sim_normdens_blue, nsteps_smooth, nsteps_smooth)
			plot(rev(range_ra), range_dec,
				 xlim = rev(range_ra), 
     			 ylim = range_dec,  
 			     xlab = "Right Ascension [deg]", 
				 ylab = "Declination [deg]",
 				 type = "n",
 				 main = paste("2D density, KNN, K=", k_vector[f], 
				 			  ", all (simulated) GCs, i=", i,  sep = ""), 
		 		 axes = F)
			image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
				  z = sim_matrix_normdens_gc_blue,
				  col = rev(heat.colors(100)), 
				  add = T)
		    points(sim_radec_total_blue[[i]][, 1], sim_radec_total_blue[[i]][, 2], 
		   		   pch = 20,
		   		   cex = 0.3,
		   		   col = "darkgray")
			# Drawing the polygons of the HST fields. 
			plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   		   "black", 0)								   
			# Drawing arrows indicating the North and West directions. 	 
			plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)		   		   						   
			# Shape of the first galaxy.
			plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
						  galaxy1_angle, cex_center = cex_center_galaxy)
			## Labels on the plot.			  
			#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[f], 0)	
			# Plotting the angular scale.
			plot_scale("topright", 180, 0)			
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
			# Calculating the number of pixels of the i-th density maps with residuals 
			# larger than 1, 2 and 3 sigmas (where sigmas is the standard deviation of
			# the total distributions of simulated density maps).
			num_sigmas_simulation_all[[i]] = vector(mode = "numeric", length = nsteps_smooth*nsteps_smooth)
			num_sigmas_simulation_red[[i]] = vector(mode = "numeric", length = nsteps_smooth*nsteps_smooth)
			num_sigmas_simulation_blue[[i]] = vector(mode = "numeric", length = nsteps_smooth*nsteps_smooth)
			for (k in seq(nsteps_smooth*nsteps_smooth)) {
				num_sigmas_simulation_all[[i]][k] = (sim_dens_all[k] - mean(sapply(sim_dens_point_all, "[", i = k, j = f)))/
										  		  	sd(sapply(sim_dens_point_all, "[", i = k, j = f))
				num_sigmas_simulation_red[[i]][k] = (sim_dens_red[k] - mean(sapply(sim_dens_point_red, "[", i = k, j = f)))/
												  	sd(sapply(sim_dens_point_red, "[", i = k, j = f))						  		  
				num_sigmas_simulation_blue[[i]][k] = (sim_dens_blue[k] - mean(sapply(sim_dens_point_blue, "[", i = k, j = f)))/
												  	 sd(sapply(sim_dens_point_blue, "[", i = k, j = f))
				}
			num_more1sigmas_simulation_all[[f]][i] = 100*length(which(num_sigmas_simulation_all[[i]] > 1 | 
																  	  num_sigmas_simulation_all[[i]] < -1))/(nsteps_smooth*nsteps_smooth)
			num_more1sigmas_simulation_red[[f]][i] = 100*length(which(num_sigmas_simulation_red[[i]] > 1 | 
																	  num_sigmas_simulation_red[[i]] < -1))/(nsteps_smooth*nsteps_smooth)
			num_more1sigmas_simulation_blue[[f]][i] = 100*length(which(num_sigmas_simulation_blue[[i]] > 1 | 
																	   num_sigmas_simulation_blue[[i]] < -1))/(nsteps_smooth*nsteps_smooth)									  	 				 
			num_more2sigmas_simulation_all[[f]][i] = 100*length(which(num_sigmas_simulation_all[[i]] > 2 | 
																	  num_sigmas_simulation_all[[i]] < -2))/(nsteps_smooth*nsteps_smooth)
			num_more2sigmas_simulation_red[[f]][i] = 100*length(which(num_sigmas_simulation_red[[i]] > 2 | 
																	  num_sigmas_simulation_red[[i]] < -2))/(nsteps_smooth*nsteps_smooth)
			num_more2sigmas_simulation_blue[[f]][i] = 100*length(which(num_sigmas_simulation_blue[[i]] > 2 | 
																	   num_sigmas_simulation_blue[[i]] < -2))/(nsteps_smooth*nsteps_smooth)			
			num_more3sigmas_simulation_all[[f]][i] = 100*length(which(num_sigmas_simulation_all[[i]] > 3 | 
																	  num_sigmas_simulation_all[[i]] < -3))/(nsteps_smooth*nsteps_smooth)
			num_more3sigmas_simulation_red[[f]][i] = 100*length(which(num_sigmas_simulation_red[[i]] > 3 | 
																	  num_sigmas_simulation_red[[i]] < -3))/(nsteps_smooth*nsteps_smooth)
			num_more3sigmas_simulation_blue[[f]][i] = 100*length(which(num_sigmas_simulation_blue[[i]] > 3 | 
																	   num_sigmas_simulation_blue[[i]] < -3))/(nsteps_smooth*nsteps_smooth)
			}
		cat(paste("		So far so good! KNN densities of the ", i, "-th simulated coordinates plotted!\n", sep = ""))
		# Calculating the number of sigmas of difference between
		# the observed value of the 2D densities of GCs and the
		# distribution of simulated 2D densities for each pixel
		# of the 2D map.
		for (k in seq(nsteps_smooth*nsteps_smooth)) {
			num_sigmas_all[[f]][k] = (dens_all[k] - mean(sapply(sim_dens_point_all, "[", i = k, j = f)))/
									  sd(sapply(sim_dens_point_all, "[", i = k, j = f))
			num_sigmas_red[[f]][k] = (dens_red[k] - mean(sapply(sim_dens_point_red, "[", i = k, j = f)))/
									  sd(sapply(sim_dens_point_red, "[", i = k, j = f))
			num_sigmas_blue[[f]][k] = (dens_blue[k] - mean(sapply(sim_dens_point_blue, "[", i = k, j = f)))/
									  sd(sapply(sim_dens_point_blue, "[", i = k, j = f))
			perc_tail_all[[f]][k] = pnorm(dens_all[k], mean = mean(sapply(sim_dens_point_all, "[", i = k, j = f)), 
									  	  sd = sd(sapply(sim_dens_point_all, "[", i = k, j = f)), lower.tail = TRUE)
			perc_tail_red[[f]][k] = pnorm(dens_red[k], mean = mean(sapply(sim_dens_point_red, "[", i = k, j = f)),
									  	  sd = sd(sapply(sim_dens_point_red, "[", i = k, j = f)), lower.tail = TRUE)
			perc_tail_blue[[f]][k] = pnorm(dens_blue[k], mean = mean(sapply(sim_dens_point_blue, "[", i = k, j = f)),
									       sd = sd(sapply(sim_dens_point_blue, "[", i = k, j = f)), lower.tail = TRUE)
			} 	
		perc_more1sigmas_all[f] = 100*length(which(num_sigmas_all[[f]] > 1 | 
												   num_sigmas_all[[f]] < -1))/(nsteps_smooth*nsteps_smooth)
		perc_more1sigmas_red[f] = 100*length(which(num_sigmas_red[[f]] > 1 | 
												   num_sigmas_red[[f]] < -1))/(nsteps_smooth*nsteps_smooth)
		perc_more1sigmas_blue[f] = 100*length(which(num_sigmas_blue[[f]] > 1 | 
											    	num_sigmas_blue[[f]] < -1))/(nsteps_smooth*nsteps_smooth)		
		}
	cat("So far so good! KNN densities of the simulations plotted!\n")
	}
	
## Plots of the averaged residual maps of the spatial
## distributions of GCs.
idx_resthresh_gc_all = list() 
idx_resthresh_gc_red = list(); idx_resthresh_gc_blue = list()
numgcs_all_pos = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
numgcs_red_pos = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
numgcs_blue_pos = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
numgcs_all_neg = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
numgcs_red_neg = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
numgcs_blue_neg = array(0, dim = c(length(k_vector), length(num_sigmas_threshold)))
idx_residuals_matrix_all_pos = list(); idx_residuals_matrix_all_neg = list()
idx_residuals_matrix_red_pos = list(); idx_residuals_matrix_red_neg = list()
idx_residuals_matrix_blue_pos = list(); idx_residuals_matrix_blue_neg = list()
structure_image_res2d_all = list()
structure_image_res2d_red = list()
structure_image_res2d_blue = list()
blurred_structure_image_res2d_all = list()
blurred_structure_image_res2d_red = list()
blurred_structure_image_res2d_blue = list()
for (f in 1:length(k_vector)) {
	residuals_matrix_all = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_red = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	residuals_matrix_blue = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))
	errresiduals_matrix_all_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_red_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_blue_aux = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_all = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_red = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))	
	errresiduals_matrix_blue = array(0, dim = c(num_simulations, nsteps_smooth, nsteps_smooth))			
	residuals_matrix_all = apply(res_all[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_red = apply(res_red[[f]], c(2, 3), sum)/num_simulations
	residuals_matrix_blue = apply(res_blue[[f]], c(2, 3), sum)/num_simulations
	errresiduals_matrix_all_aux = apply(errres_all[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_red_aux = apply(errres_red[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_blue_aux = apply(errres_blue[[f]]^2, c(2, 3), sum)
	errresiduals_matrix_all = sqrt(errresiduals_matrix_all_aux)/num_simulations
	errresiduals_matrix_red = sqrt(errresiduals_matrix_red_aux)/num_simulations
	errresiduals_matrix_blue = sqrt(errresiduals_matrix_blue_aux)/num_simulations
	## Plots of the average residuals of the simulated 2D densities evaluated with 
	## the KNN method.
	radeg_lims_maps = rev(range(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, 1]))
	decdeg_lims_maps = range(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, 2])
	radeg_lims_maps[2] = radeg_lims_maps[2] - 0.01
	# Average residuals of all (simulated) GCs 2D density relative to the observed
	# 2D density of all GCs.
	res_threshold_gc_all = lseq(max(residuals_matrix_all)/100, max(residuals_matrix_all), length.out = 10)[index_res_threshold]
	idx_resthresh_gc_all[[f]] = which(residuals_matrix_all >= res_threshold_gc_all)
	idx_residuals_matrix_all_pos[[f]] = list(); idx_residuals_matrix_all_neg[[f]] = list()
	for (z in seq(length(num_sigmas_threshold))) {
		idx_residuals_matrix_all_pos[[f]][[z]] = which(num_sigmas_all[[f]] >= num_sigmas_threshold[z])
		idx_residuals_matrix_all_neg[[f]][[z]] = which(num_sigmas_all[[f]] <= -num_sigmas_threshold[z])
		}
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_all_K_", 
											 		   k_vector[f], "_thresholds_sigmas_largestructures.pdf", sep = "")	
			} else if (typechull_flag == "large_intermediate") {
				pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_all_K_", 
												 		   k_vector[f], "_thresholds_sigmas_largeintermediatestructures.pdf", sep = "")
				} 
		} else {
			pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_all_K_", 
													   k_vector[f], "_thresholds_sigmas.pdf", sep = "")
			}				
	pdf(pathpdfdensknnallres_simthresholds, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, ], 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
		 xlim = rev(range_ra), 
		 ylim = range_dec, 
		 col = "black",
		 pch = 15, 
		 cex = 1,
		 main = "", 
		 type = "n",
		 axes = F)	
	structure_image_res2d_all[[f]] = image_residual_maps(num_sigmas_all[[f]], 
														 radeg_gridpoint_all, decdeg_gridpoint_all, 
														 num_sigmas_threshold, 
														 basic_colors_residual)
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		  z = structure_image_res2d_all[[f]]$image,
		  levels = structure_image_res2d_all[[f]]$levels,
		  col = structure_image_res2d_all[[f]]$colors,
		  add = TRUE)		
	text(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 1], 
		 expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 2],
		 round(num_sigmas_all[[f]], 1), 
		 cex = 0.5)		 		  	 	 	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[3]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[3]], idx_within_footprint), 2], 
	#	   pch = "-", 
	#	   cex = 1.8, 
	#	   col = "black")
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[3]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[3]], idx_within_footprint), 2], 
	#	   pch = "-", 
	#	   cex = 1.8, 
	#	   col = "white")	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[4]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[4]], idx_within_footprint), 2], 
	#	   pch = "I", 
	#	   cex = 1.5, 
	#	   col = "black")	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[4]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[4]], idx_within_footprint), 2], 
	#	   pch = "I", 
	#	   cex = 1.5, 
	#	   col = "white")	
	# Graying-out the area outside of the convex hull(s)
	# of the large structure(s), if requested.	    
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			chull_coordinates_external = externalify(chull_coordinates[[1]][[1]]$x, 
													 chull_coordinates[[1]][[1]]$y)
			polygon(chull_coordinates_external[, 1], 
					chull_coordinates_external[, 2], 
					col = addalpha("darkgray", 180), 
					lwd = 3, 
					border = NA)	
			} else if (typechull_flag == "large_intermediate") {
				for (k in 1:length(chull_coordinates$xy_coords_contour)) {
					polygon(chull_coordinates$xy_coords_contour[[k]][[1]]$x, 
							chull_coordinates$xy_coords_contour[[k]][[1]]$y, 
							border = "black", 
							lwd = 6, 
							col = NA)
			 		}
				}	 		 				
		}	
	# Drawing the polygons of the footprint of the observations. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)	
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  				   							  			 
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	## Labels on the plot.			  
	#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[f], 0)
	# Plotting points.	
	if (points_average_residuals == "yes") {
		points(radeg_gc, decdeg_gc,      
			   pch = luminclass_gc, 
			   cex = cex_gc, 
			   col = colorclass_gc)      
		}	
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)	
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)				
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)		
	# Graying-out the area outside of the convex hull(s)
	# of the large structure(s), if requested.	    
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			text(187.38, 8.07, 
				 "A1", 
				 cex = 3, 
				 col = "darkgreen")	
			} else if (typechull_flag == "large_intermediate") {
				for (k in 1:length(chull_coordinates$xy_coords_contour)) {
					text(mean(chull_coordinates$xy_coords_contour[[k]][[1]]$x), 
						 mean(chull_coordinates$xy_coords_contour[[k]][[1]]$y),  
			 			 paste("A", k, sep = ""), 
			 			 cex = 3, 
			 			 col = "darkgreen")
			 		print(paste("Number of GCs in excess in structure ", "A", k, ": ", 
			 					chull_coordinates$number_gc_excess[[k]], sep = ""))	 
			 		}
				}	 		 				
		}	
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
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			pathpdfdensknnallres_simthresholds_blurred = paste(path_plot, "res2d_all_K_", 
											 		   		   k_vector[f], "_thresholds_sigmas_largestructures_blurred.pdf", sep = "")	
			} else if (typechull_flag == "large_intermediate") {
				pathpdfdensknnallres_simthresholds_blurred = paste(path_plot, "res2d_all_K_", 
												 		   		   k_vector[f], "_thresholds_sigmas_largeintermediatestructures_blurred.pdf", sep = "")
				} 
		} else {
			pathpdfdensknnallres_simthresholds_blurred = paste(path_plot, "res2d_all_K_", 
													   		   k_vector[f], "_thresholds_sigmas_blurred.pdf", sep = "")
			}
	pdf(pathpdfdensknnallres_simthresholds_blurred, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, ], 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
		 xlim = rev(range_ra), 
		 ylim = range_dec, 
		 col = "black",
		 pch = 15, 
		 cex = 1,
		 main = "", 
		 type = "n",
		 axes = F)	
	blurred_structure_image_res2d_all[[f]] = image_residual_maps_blurred(num_sigmas_all[[f]], 
														 				 radeg_gridpoint_all, decdeg_gridpoint_all, 
														 				 num_sigmas_threshold, 
														 				 basic_colors_residual, 
														 				 1)
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		  z = blurred_structure_image_res2d_all[[f]]$image,
		  levels = blurred_structure_image_res2d_all[[f]]$levels,
		  col = blurred_structure_image_res2d_all[[f]]$colors,
		  add = TRUE)		
	text(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 1], 
		 expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 2],
		 #round(num_sigmas_all[[f]], 1), 
		 as.vector(blurred_structure_image_res2d_all[[f]]$blurred_sigmas),
		 cex = 0.5)		 		  	 	 	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[3]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[3]], idx_within_footprint), 2], 
	#	   pch = "-", 
	#	   cex = 1.8, 
	#	   col = "black")
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[3]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[3]], idx_within_footprint), 2], 
	#	   pch = "-", 
	#	   cex = 1.8, 
	#	   col = "white")	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[4]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_pos[[f]][[4]], idx_within_footprint), 2], 
	#	   pch = "I", 
	#	   cex = 1.5, 
	#	   col = "black")	 
	#points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[4]], idx_within_footprint), 1], 
	#	   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_all_neg[[f]][[4]], idx_within_footprint), 2], 
	#	   pch = "I", 
	#	   cex = 1.5, 
	#	   col = "white")	
	# Graying-out the area outside of the convex hull(s)
	# of the large structure(s), if requested.	    
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			chull_coordinates_external = externalify(chull_coordinates[[1]][[1]]$x, 
													 chull_coordinates[[1]][[1]]$y)
			polygon(chull_coordinates_external[, 1], 
					chull_coordinates_external[, 2], 
					col = addalpha("darkgray", 180), 
					lwd = 3, 
					border = NA)	
			} else if (typechull_flag == "large_intermediate") {
				for (k in 1:length(chull_coordinates$xy_coords_contour)) {
					polygon(chull_coordinates$xy_coords_contour[[k]][[1]]$x, 
							chull_coordinates$xy_coords_contour[[k]][[1]]$y, 
							border = "black", 
							lwd = 6, 
							col = NA)
			 		}
				}	 		 				
		}	
	# Drawing the polygons of the footprint of the observations. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)	
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  				   							  			 
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	## Labels on the plot.			  
	#plot_labels("topleft", "GCs", 0.05, 1.4, k_vector[f], 0)
	# Plotting points.	
	if (points_average_residuals == "yes") {
		points(radeg_gc, decdeg_gc,      
			   pch = luminclass_gc, 
			   cex = cex_gc, 
			   col = colorclass_gc)      
		}	
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)	
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)				
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)		
	# Graying-out the area outside of the convex hull(s)
	# of the large structure(s), if requested.	    
	if (chull_flag == "yes") {
		if (typechull_flag == "large") {
			text(187.38, 8.07, 
				 "A1", 
				 cex = 3, 
				 col = "darkgreen")	
			} else if (typechull_flag == "large_intermediate") {
				for (k in 1:length(chull_coordinates$xy_coords_contour)) {
					text(mean(chull_coordinates$xy_coords_contour[[k]][[1]]$x), 
						 mean(chull_coordinates$xy_coords_contour[[k]][[1]]$y),  
			 			 paste("A", k, sep = ""), 
			 			 cex = 3, 
			 			 col = "darkgreen")
			 		print(paste("Number of GCs in excess in structure ", "A", k, ": ", 
			 					chull_coordinates$number_gc_excess[[k]], sep = ""))	 
			 		}
				}	 		 				
		}	
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
	# Plot of the average residuals of the red (simulated) GCs 2D density relative to the
	# observed red GCs 2D density.		
	res_threshold_gc_red = lseq(max(residuals_matrix_red)/100, max(residuals_matrix_red), length.out = 10)[index_res_threshold]
	idx_resthresh_gc_red[[f]] = which(residuals_matrix_red >= res_threshold_gc_red)
	idx_residuals_matrix_red_pos[[f]] = list(); idx_residuals_matrix_red_neg[[f]] = list()
	for (z in seq(length(num_sigmas_threshold))) {
		idx_residuals_matrix_red_pos[[f]][[z]] = which(num_sigmas_red[[f]] >= num_sigmas_threshold[z])
		idx_residuals_matrix_red_neg[[f]][[z]] = which(num_sigmas_red[[f]] <= -num_sigmas_threshold[z])
		}
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	if (chull_flag != "yes") {
		pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_red_K_", 
										 		   k_vector[f], "_thresholds_sigmas.pdf", sep = "")	
		} else {
			pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_red_K_", 
											 		   k_vector[f], "_thresholds_sigmas_largestructures.pdf", sep = "")
			}
	pdf(pathpdfdensknnallres_simthresholds, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, ], 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
		 xlim = rev(range_ra), 
		 ylim = range_dec, 
		 col = "black",
		 pch = 15, 
		 cex = 1,
		 main = "", 
		 type = "n",
		 axes = F)	
	structure_image_res2d_red[[f]] = image_residual_maps(num_sigmas_red[[f]], 
														 radeg_gridpoint_all, decdeg_gridpoint_all, 
														 num_sigmas_threshold, 
														 basic_colors_residual)
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		  z = structure_image_res2d_red[[f]]$image,
		  levels = structure_image_res2d_red[[f]]$levels,
		  col = structure_image_res2d_red[[f]]$colors,
		  add = TRUE)		
	#text(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 1], 
	#	 expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 2],
	#	 round(num_sigmas_red[[f]], 1), 
	#	 cex = 0.5)		 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_pos[[f]][[3]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_pos[[f]][[3]], idx_within_footprint), 2], 
		   pch = "-", 
		   cex = 1.8, 
		   col = "black")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_neg[[f]][[3]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_neg[[f]][[3]], idx_within_footprint), 2], 
		   pch = "-", 
		   cex = 1.8, 
		   col = "white")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_pos[[f]][[4]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_pos[[f]][[4]], idx_within_footprint), 2], 
		   pch = "I", 
		   cex = 1.5, 
		   col = "black")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_neg[[f]][[4]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_red_neg[[f]][[4]], idx_within_footprint), 2], 
		   pch = "I", 
		   cex = 1.5, 
		   col = "white")	 
	# Drawing the polygons of the HST fields. 
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)	
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  				   
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)				
	## Labels on the plot.			  
	#plot_labels("topleft", "Red GCs", 0.05, 1.4, k_vector[f], 0)	
	# Plotting points.
	if (points_average_residuals == "yes") {
		points(radeg_gc, decdeg_gc,      
			   pch = luminclass_gc, 
			   cex = cex_gc, 
			   col = colorclass_gc)      
		}			
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)			
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)			
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
	# Plot of the average residuals of the blue (simulated) GCs 2D density relative to the
	# observed red GCs 2D density.		
	res_threshold_gc_blue = lseq(max(residuals_matrix_blue)/100, max(residuals_matrix_blue), length.out = 10)[index_res_threshold]
	idx_resthresh_gc_blue[[f]] = which(residuals_matrix_blue >= res_threshold_gc_blue)
	idx_residuals_matrix_blue_pos[[f]] = list(); idx_residuals_matrix_blue_neg[[f]] = list()
	for (z in seq(length(num_sigmas_threshold))) {
		idx_residuals_matrix_blue_pos[[f]][[z]] = which(num_sigmas_blue[[f]] >= num_sigmas_threshold[z])
		idx_residuals_matrix_blue_neg[[f]][[z]] = which(num_sigmas_blue[[f]] <= -num_sigmas_threshold[z])
		}
	# Maps of the pixels >1sigma and other symbols showing the pixels with 2sigmas
	# and 3sigmas significance.
	if (chull_flag != "yes") {
		pathpdfdensknnallres_simthresholds = paste(path_plot, "res2d_blue_K_", 
										 		   k_vector[f], "_thresholds_sigmas.pdf", sep = "")	
		} else {
			pathpdfdensknnallres_simthresholds = paste(path_plot,"res2d_blue_K_", 
											 		   k_vector[f], "_thresholds_sigmas_largestructures.pdf", sep = "")
			}
	pdf(pathpdfdensknnallres_simthresholds, width = 7, height = 7, paper = "special")
	par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
		mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1, 0, 0))
	plot(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[idx_within_footprint, ], 
		 xlab = "Right Ascension [deg]", 
		 ylab = "Declination [deg]",
		 xlim = rev(range_ra), 
		 ylim = range_dec, 
		 col = "black",
		 pch = 15, 
		 cex = 1,
		 main = "", 
		 type = "n",
		 axes = F)	
	structure_image_res2d_blue[[f]] = image_residual_maps(num_sigmas_blue[[f]], 
														 radeg_gridpoint_all, decdeg_gridpoint_all, 
														 num_sigmas_threshold, 
														 basic_colors_residual)
	image(x = radeg_gridpoint_all, y = decdeg_gridpoint_all, 
		  z = structure_image_res2d_blue[[f]]$image,
		  levels = structure_image_res2d_blue[[f]]$levels,
		  col = structure_image_res2d_blue[[f]]$colors,
		  add = TRUE)		
	#text(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 1], 
	#	 expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[, 2],
	#	 round(num_sigmas_blue[[f]], 2), 
	#	 cex = 0.5)		 		  	 	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_pos[[f]][[3]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_pos[[f]][[3]], idx_within_footprint), 2], 
		   pch = "-", 
		   cex = 1.8, 
		   col = "black")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_neg[[f]][[3]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_neg[[f]][[3]], idx_within_footprint), 2], 
		   pch = "-", 
		   cex = 1.8, 
		   col = "white")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_pos[[f]][[4]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_pos[[f]][[4]], idx_within_footprint), 2], 
		   pch = "I", 
		   cex = 1.5, 
		   col = "black")	 
	points(expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_neg[[f]][[4]], idx_within_footprint), 1], 
		   expand.grid(radeg_gridpoint_all, decdeg_gridpoint_all)[intersect(idx_residuals_matrix_blue_neg[[f]][[4]], idx_within_footprint), 2], 
		   pch = "I", 
		   cex = 1.5, 
		   col = "white")	 
	# Drawing the polygons of the HST fields. 
	plot_footprint(polygons, polygons_external, polygons_external_avoid, 1, 
				   "white", 255)	
	# Background.
	polygon(glob_poly_back_ra, glob_poly_back_dec, col = addalpha("darkgray", 155))  
	# Shape of the first galaxy.
	plot_galaxies(galaxy1_center, galaxy1_isodiameter, galaxy1_isodiameter_ratio, 
				  galaxy1_angle, cex_center = cex_center_galaxy)
	## Labels on the plot.			  
	#plot_labels("topleft", "Blue GCs", 0.05, 1.4, k_vector[f], 0)	
	# Plotting points.
	if (points_average_residuals == "yes") {
		points(radeg_gc, decdeg_gc,      
			   pch = luminclass_gc, 
			   cex = cex_gc, 
			   col = colorclass_gc)      
		}	
	# Drawing arrows indicating the North and West directions. 	 
	plot_arrows_new("bottomleft", ra_length_arrow, dec_length_arrow)	
	# Grid of lines.
	plot_line_grid(range_ra, range_dec, 0.05)		
	# Plotting the angular scale.
	plot_scale("topright", 180, 0)		
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
cat("So far so good! Plots of the 2D residual maps of simulations for all, red, blue, high-L and low-L GCs created!\n")

## Saving the whole session into a file with a name understandable enough 
## to be read in into R later for plots etc.
if (sim_flag == "yes") {
	path_data_sessions = paste(path_data, galaxy, "/Sessions/", sep = "")
	dir.create(path_data_sessions, showWarnings = FALSE)	
	name_file_session = paste(path_data_sessions, "gc_nsims_", 
							  num_simulations, "_nsteps_", nsteps_smooth, 
							  "_simtype_", sim_type, ".Rdata", sep = "")
	save.image(name_file_session, compress = TRUE)
	cat("So far so good! Image of the session saved!\n")
	}