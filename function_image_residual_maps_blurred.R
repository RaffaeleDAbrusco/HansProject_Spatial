## Function that produces an image of values that can be used to plot
## the blurred map of the 2D residuals of the spatial distribution of GCs/LMXBs. 
##
## Written by R. D'Abrusco
##
## Last modified: 31/8/16.

image_residual_maps_blurred <- function(vector_sigmas, 
										ra_grid, dec_grid, 
										thresholds_sigmas,  
										basic_colors_residual,
										sigma_blur = 1) {
										
	# Loading the package required to blur the image.
	suppressMessages(library(spatstat, quietly = TRUE))										
  	
  	# Creating an image from the vector of sigmas values, 
  	# with size equal to the size of the RA vs Dec grid.
 	residual_map = matrix(vector_sigmas, nrow = length(ra_grid), ncol = length(dec_grid))

	# Applying the "blur", i.e. applying a gaussian smoothing with 
	# sigma equal to sigma_blur. 
	image_blurred_residual_map = blur(as.im(residual_map), sigma_blur)
	blurred_residual_map = as.matrix(image_blurred_residual_map)
 	print(range(blurred_residual_map))
 	# Defining a new matrix that will contain the values required to plot
	# the residual map as an image.
	image_residual_map = matrix(0, nrow = length(ra_grid), ncol = length(dec_grid))
	
	# Populating the image_residual_map matrix.
	for (f in 1:(length(thresholds_sigmas) - 1)) {
		print(f)
		image_residual_map[which(blurred_residual_map >= thresholds_sigmas[f] & blurred_residual_map < thresholds_sigmas[f + 1])] = f
		image_residual_map[which(blurred_residual_map <= -thresholds_sigmas[f] & blurred_residual_map > -thresholds_sigmas[f + 1])] = -f
		}
	image_residual_map[which(blurred_residual_map >= thresholds_sigmas[length(thresholds_sigmas)])] = length(thresholds_sigmas)
	image_residual_map[which(blurred_residual_map <= -thresholds_sigmas[length(thresholds_sigmas)])] = -length(thresholds_sigmas)	
	
	# Determining the levels that are visible in the image of the residual map.
	all_levels = seq(-length(thresholds_sigmas), length(thresholds_sigmas), by = 1)
	levels_image = all_levels[sort(unique(c(image_residual_map))) + length(thresholds_sigmas) + 1]
	
  	# Creating the vector of colors that will be used to plot
  	# the residual images.
  	num_levels_color = 2*length(thresholds_sigmas) + 1
  	levels_color = vector(mode = "character", length = num_levels_color) 
  	levels_color[length(thresholds_sigmas) + 1] = basic_colors_residual[2]
  	levels_color[1] = basic_colors_residual[1]; levels_color[num_levels_color] = basic_colors_residual[3]
  	levels_color[2:length(thresholds_sigmas)] = addalpha(basic_colors_residual[1], seq(255, 0, length.out = length(thresholds_sigmas) + 1)[2:(length(thresholds_sigmas))])
  	levels_color[(length(thresholds_sigmas) + 2):(num_levels_color - 1)] = addalpha(basic_colors_residual[3], seq(0, 255, length.out = length(thresholds_sigmas) + 1)[2:(length(thresholds_sigmas))])	
	
	# Selecting only the colors associated to the levels that will 
	# appear in the residual image.
	levels_color_image = levels_color[levels_image + length(thresholds_sigmas) + 1]
	
	# Defining a structure containing everything needed to create the residual image
	structure_image_residual = list("image" = image_residual_map, 
									"levels" = levels_image, 
									"colors" = levels_color_image, 
									"blurred_sigmas" = round(blurred_residual_map, 1))
	return(structure_image_residual)
	}
	
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red >= -0.5 & matrix_num_sigmas_red <= 0.5)] = 0
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red < -0.5 & matrix_num_sigmas_red >= -1)] = -1
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red < -1 & matrix_num_sigmas_red >= -2)] = -2
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red < -2 & matrix_num_sigmas_red >= -3)] = -3
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red < -3)] = -4		
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red > 0.5 & matrix_num_sigmas_red < 1)] = 1
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red >= 1 & matrix_num_sigmas_red < 2)] = 2
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red >= 2 & matrix_num_sigmas_red < 3)] = 3
#	image_matrix_num_sigmas_red[which(matrix_num_sigmas_red > 3)] = 4			
	