## This program defines functions that will be used
## by external programs to evaluate whether a given
## set of point is located within or ourside of a
## generic polygon.
##
## Written by R. D'Abrusco
##
## Last modified: 27/8/15.
#
# 	27/8/15			->		Added the function "area_intersections_variable_isophotes"
#							that evaluates the area of the region between two ellipses
#							representing two isophotes at different radial distances, 
#							with variable position angle.


## Defining functions.
point_in_polygon <- function(polygon, p) {
	# Return whether a given point is outside or
	# inside a polygon.
	count <- 0
	for(side in polygon) {
		if (ray_intersect_segment(p, side)) {
			count <- count + 1
			}
		}
	  	if (count %% 2 == 1)
			1
	  	else
			0
	}
 
ray_intersect_segment <- function(p, side) {
	# Determine if the half-line originating 
	# in a generic point p intersects one of the 
	# sides of the polygon. It is used by the 
	# previous function to determine if the point
	# is within the polygon or not.
  	eps <- 0.0001
  	a <- side$A
  	b <- side$B
  	if (a$y > b$y) {
    	a <- side$B
    	b <- side$A
  		}
  	if ((p$y == a$y) || (p$y == b$y)) {
    	p$y <- p$y + eps
  		}
  	if ((p$y < a$y) || (p$y > b$y))
    	return(FALSE)
  	else if (p$x > max(a$x, b$x))
    	return(FALSE)
  	else {
    	if (p$x < min(a$x, b$x))
      		return(TRUE)
    	else {
      		if (a$x != b$x)
        		m_red <- (b$y - a$y)/(b$x - a$x)
      		else
        		m_red <- Inf
      		if (a$x != p$x)
       	 		m_blue <- (p$y - a$y)/(p$x - a$x)
      		else
        		m_blue <- Inf
      		return(m_blue >= m_red)
    		}
  		}
	}
	
## Defining a new function that takes x and y coordinates
## of a point and convert them to the elements of a list, 
# named "x" and "y" respectively. 
point <- function(x, y) list(x = x, y = y)
pointnew <- function(xy) list(x = xy[1], y = xy[2])	
	
gridknots_inside_polygon <- function(grid_knots_x, grid_knots_y, xy_polygon) {
	# Return the indices of the points of a
	# regular grid contained within a polygon, 
	# using the other functions defined for the 
	# ray casting.
	# grid_knots_x 		<-	vector containing the x coordinates
	#						of the knots of the grid.
	# grid_knots_y 		<-	vector containing the y coordinates
	#						of the knots of the grid.
	# xy_polygon		<-  array (n rows x 2 columns) containing
	#						the x and y coordinates of the vertices
	#						of the polygon.
	xy_grid = expand.grid(grid_knots_x, grid_knots_y)
	xy_grid_poly = apply(xy_grid, 1, pointnew)
	flag_insideoutside = vector(mode = "numeric", length = length(xy_grid[, 1]))
	points_polygon = apply(xy_polygon, 1, pointnew)
	polygon = createPolygonNew(points_polygon)	
	for(p in 1:length(xy_grid[, 1])) {
		flag_insideoutside[p] = point_in_polygon(polygon, point(xy_grid[p, 1], xy_grid[p, 2]))
   		}
	idx_inside = which(flag_insideoutside == 1)
	return(idx_inside)
	}

points_inside_polygon <- function(points_x, 
								  points_y, 
								  xy_polygon) {
	# Return the indices of the points contained within an irregular polygon, 
	# using the other functions defined for the ray casting.
	# points_x 			<-	vector containing the x coordinates
	#						of the points
	# points_y 			<-	vector containing the y coordinates
	#						of the points.
	# xy_polygon		<-  array (n rows x 2 columns) containing
	#						the x and y coordinates of the vertices
	#						of the polygon.
	if (length(points_x) == 0) {
		idx_inside = vector(mode = "numeric", length = 0)
		return(idx_inside)
		}
	xy_points = cbind(points_x, points_y)
	xy_points_poly = apply(xy_points, 1, pointnew)
	flag_insideoutside = vector(mode = "numeric", length = length(xy_points[, 1]))
	points_polygon = apply(xy_polygon, 1, pointnew)
	polygon = createPolygonNew(points_polygon)	
	for(p in 1:length(xy_points[, 1])) {
		flag_insideoutside[p] = point_in_polygon(polygon, point(xy_points[p, 1], xy_points[p, 2]))
   		}
	idx_inside = which(flag_insideoutside == 1)
	return(idx_inside)
	}

createPolygonNew <- function(pts) {
	# Takes list of points coordinates and list of 
	# indices couples. 
	pol <- list()
  	#for(pseg in segs) {
    # 	pol <- c(pol, list(list(A = pts[[pseg[1]]], B = pts[[pseg[2]]])))
  	#	}
	for(i in 1:(length(pts) - 1)) {
    	pol <- c(pol, list(list(A = pts[[i]], B = pts[[i + 1]])))
  		}  	
  	pol
	}
	
area_intersections <- function(x_span, 				# Interval along the x direction (R.A.)
							   y_span, 				# Interval along the y direction (Dec.)
							   num_knots_grid, 		# Number of knots of the grid.
							   polygon,				# Irregular polygon representing the footpring of the observations.
							   galaxy_center1, 		# Coordinates of the center of the primary galaxy.
							   galaxy_angle1, 		# Inclination of the main axis of the primary galaxy.
							   axes_larger, 		# Major and minor axes of the outside ellipse of the 
							   						# main galaxy.
							   axes_smaller,		# Major and minor axes of the inside ellipse of the 
							   						# main galaxy.
							   galaxy_center2, 		# Coordinates of the center of the secondary galaxy.
							   galaxy_angle2, 		# Inclination of the main axis of the secondary galaxy.
							   axes_secondary) {	# Major and minor axes of the ellipse of the 
							   						# secondary ellipse.
	# Defining the center of the cells of the grid that will be 
	# used to evaluate the area of the intersections.
	x_seeds = seq(x_span[1], x_span[2], length.out = num_knots_grid)	
	y_seeds = seq(y_span[1], y_span[2], length.out = num_knots_grid)						
	# Defining the grid.
	xy_grid = expand.grid(x_seeds, y_seeds)							
	# Using functions already existing to extract the points of the 
	# grid within the irregular polygon.						
	idx_inside_polygon = gridknots_inside_polygon(x_seeds, y_seeds, polygon) 		
	# Preliminary transformations.
	shifted_x1 = xy_grid[, 1] - galaxy_center1[1]
	shifted_y1 = xy_grid[, 2] - galaxy_center1[2]
	major_axis_larger = axes_larger[1]; minor_axis_larger = axes_larger[2]
	major_axis_smaller = axes_smaller[1]; minor_axis_smaller = axes_smaller[2]
	if (missing(galaxy_center2) == FALSE) {
		shifted_x2 = xy_grid[, 1] - galaxy_center2[1]
		shifted_y2 = xy_grid[, 2] - galaxy_center2[2]
		major_axis_secondary = axes_secondary[1]; minor_axis_secondary = axes_secondary[2]
		}	
	locus_ellipse1 = vector(mode = "numeric", length(xy_grid[, 1]))
	locus_ellipse2 = vector(mode = "numeric", length(xy_grid[, 1]))
	# Extracting indices of pixels within the larger primary ellipse.
	locus_ellipse1 = (shifted_x1*cos(pi*galaxy_angle1/180) + shifted_y1*sin(pi*galaxy_angle1/180))^2/major_axis_larger^2 + 
					 (shifted_x1*sin(pi*galaxy_angle1/180) - shifted_y1*cos(pi*galaxy_angle1/180))^2/minor_axis_larger^2
	idx_inside_ellipse_primary = which(locus_ellipse1 <= 1)
	# Extracting indices of pixels outside of the smaller primary ellipse.
	locus_ellipse2 = (shifted_x1*cos(pi*galaxy_angle1/180) + shifted_y1*sin(pi*galaxy_angle1/180))^2/major_axis_smaller^2 + 
					 (shifted_x1*sin(pi*galaxy_angle1/180) - shifted_y1*cos(pi*galaxy_angle1/180))^2/minor_axis_smaller^2
	idx_outside_ellipse_primary = which(locus_ellipse2 >= 1)
	if (missing(galaxy_center2) == FALSE) {
		# Extracting indices of pixels outside of the secondary ellipse.
		locus_ellipse3 = vector(mode = "numeric", length(xy_grid[, 1]))
		locus_ellipse3 = (shifted_x2*cos(pi*galaxy_angle2/180) + shifted_y2*sin(pi*galaxy_angle2/180))^2/major_axis_secondary^2 + 
						 (shifted_x2*sin(pi*galaxy_angle2/180) - shifted_y2*cos(pi*galaxy_angle2/180))^2/minor_axis_secondary^2
		idx_outside_ellipse_secondary = which(locus_ellipse3 >= 1)		
		}				
	# Combining the indices to find the intersection of all regions.
	idx_aux1 = intersect(idx_inside_polygon, idx_inside_ellipse_primary)
	idx = intersect(idx_aux1, idx_outside_ellipse_primary)
	if (missing(galaxy_center2) == FALSE) {
		idx = intersect(idx, idx_outside_ellipse_secondary)
		}
	# Area of the single pixel.
	area_pixel = abs(x_span[2] - x_span[1])*(y_span[2] - y_span[1])/length(xy_grid[, 1])
	area_intersection = area_pixel*length(idx)	
	stuff = list(indices = idx, area = area_intersection)
	return(stuff)
	}	

area_intersections_variable_isophotes <- function(x_span, 			# Interval along the x direction (R.A.)
							   					  y_span, 			# Interval along the y direction (Dec.)
							   					  num_knots_grid, 	# Number of knots of the grid.
							   					  polygon,			# Irregular polygon representing the HST field.
							   					  galaxy_center, 	# Coordinates of the center of the primary galaxy.
							   					  galaxy_angles, 	# Position angles of the isophotes
							   					  axes_larger, 		# Major axes of the isophotes
							   					  axes_smaller,		# Minor axes of the isophotes
							   					  galaxy_center2, 	# Coordinates of the center of the secondary galaxy.
							   					  galaxy_angle2, 	# Inclination of the main axis of the secondary galaxy.
							   					  axes_secondary) {	# Major and minor axes of the ellipse of the 
							   										# secondary ellipse.
	# Defining the center of the cells of the grid that will be 
	# used to evaluate the area of the intersections.
	x_seeds = seq(x_span[1], x_span[2], length.out = num_knots_grid)	
	y_seeds = seq(y_span[1], y_span[2], length.out = num_knots_grid)						
	# Defining the grid.
	xy_grid = expand.grid(x_seeds, y_seeds)							
	# Using functions already existing to extract the points of the 
	# grid within the irregular polygon describing the general 
	# footprint of the observations.						
	idx_inside_polygon = gridknots_inside_polygon(x_seeds, y_seeds, polygon) 		
	# Preliminary transformations.
	shifted_x1 = xy_grid[, 1] - galaxy_center[1]
	shifted_y1 = xy_grid[, 2] - galaxy_center[2]
	major_axis_larger = axes_larger[1]
	minor_axis_larger = axes_larger[2]
	pa_angle_larger = galaxy_angles[2]
	major_axis_smaller = axes_smaller[1]
	minor_axis_smaller = axes_smaller[2]
	pa_angle_smaller = galaxy_angles[1]
	if (missing(galaxy_center2) == FALSE) {
		shifted_x2 = xy_grid[, 1] - galaxy_center2[1]
		shifted_y2 = xy_grid[, 2] - galaxy_center2[2]
		major_axis_secondary = axes_secondary[1]; minor_axis_secondary = axes_secondary[2]
		}	
	locus_ellipse1 = vector(mode = "numeric", length(xy_grid[, 1]))
	locus_ellipse2 = vector(mode = "numeric", length(xy_grid[, 1]))
	# Extracting indices of pixels within the larger primary ellipse.
	locus_ellipse1 = (shifted_x1*cos(pi*pa_angle_larger/180) + shifted_y1*sin(pi*pa_angle_larger/180))^2/major_axis_larger^2 + 
					 (shifted_x1*sin(pi*pa_angle_larger/180) - shifted_y1*cos(pi*pa_angle_larger/180))^2/minor_axis_larger^2
	idx_inside_ellipse_primary = which(locus_ellipse1 <= 1)
	# Extracting indices of pixels outside of the smaller primary ellipse.
	locus_ellipse2 = (shifted_x1*cos(pi*pa_angle_smaller/180) + shifted_y1*sin(pi*pa_angle_smaller/180))^2/major_axis_smaller^2 + 
					 (shifted_x1*sin(pi*pa_angle_smaller/180) - shifted_y1*cos(pi*pa_angle_smaller/180))^2/minor_axis_smaller^2
	idx_outside_ellipse_primary = which(locus_ellipse2 >= 1)
	if (missing(galaxy_center2) == FALSE) {
		# Extracting indices of pixels outside of the secondary ellipse.
		locus_ellipse3 = vector(mode = "numeric", length(xy_grid[, 1]))
		locus_ellipse3 = (shifted_x2*cos(pi*galaxy_angle2/180) + shifted_y2*sin(pi*galaxy_angle2/180))^2/major_axis_secondary^2 + 
						 (shifted_x2*sin(pi*galaxy_angle2/180) - shifted_y2*cos(pi*galaxy_angle2/180))^2/minor_axis_secondary^2
		idx_outside_ellipse_secondary = which(locus_ellipse3 >= 1)		
		}				
	# Combining the indices to find the intersection of all regions.
	idx_aux1 = intersect(idx_inside_polygon, idx_inside_ellipse_primary)
	idx = intersect(idx_aux1, idx_outside_ellipse_primary)
	if (missing(galaxy_center2) == FALSE) {
		idx = intersect(idx, idx_outside_ellipse_secondary)
		}
	# Area of the single pixel.
	area_pixel = abs(x_span[2] - x_span[1])*(y_span[2] - y_span[1])/length(xy_grid[, 1])
	area_intersection = area_pixel*length(idx)	
	stuff = list(indices = idx, area = area_intersection)
	return(stuff)
	}	

if (1 == 0) {
## Testing the code.
# Defining a regular grid and plotting it.
x_grid = seq(-2, 2, by = 0.1)
y_grid = x_grid
# Defining the points that will be connected
# in different ways to create the polygons.
n_vertices = 4
xy_vertices = array(0, dim = c(n_vertices, 2))
#vert1 = c(-0.93, 1.4)
#vert2 = c(0.1, 0)
#vert3 = c(-0.13, -1.23)
#vert45 = c(0.9, 1)
#vert4 = c(2, 1.9)
#vert5 = c(0.5, 0.8)
vert1 = c(-1.9, 1.9)
vert2 = c(-1.1, -1.1)
vert3 = c(0.95, 1.05)
vert4 = c(0.87, 1)
#xy_vertices = rbind(vert1, vert2, vert3, vert4, vert45, vert5, vert1)
xy_vertices = rbind(vert1, vert2, vert3, vert4, vert1)
rownames(xy_vertices) <- NULL
x_polygon = xy_vertices[, 1]
y_polygon = xy_vertices[, 2]
plot(expand.grid(x_grid, y_grid)[, 1], 
     expand.grid(x_grid, y_grid)[, 2],
     pch = 20, 
     cex = 0.5, 
     col = "darkgray", 
     xlab = "", ylab = "") 
for (i in 1:length(y_polygon)) {
	points(x_polygon[i], y_polygon[i], pch = 17, col = "red")  
	}
polygon(x_polygon, y_polygon, lwd = 2, col = NA, dens = 5)

## Creating the polygon vertices in the right format.
points_polygon = apply(xy_vertices, 1, pointnew)
pts <- list(point(0,0), point(10,0), point(10,10), point(0,10),
           point(2.5,2.5), point(7.5,2.5), point(7.5,7.5), point(2.5,7.5), 
            point(0,5), point(10,5), 
            point(3,0), point(7,0), point(7,10), point(3,10))
 
## Creating the polygons using the "createPolygon" function.
natural_indices_polygon = seq(length(xy_vertices[, 1]) + 1)
polygon = createPolygonNew(points_polygon)

## Preparing the testpoints, i.e. all points of the grid, to be 
## tested relative to the polygon built previously. 

#testpoints <-
#  list(
#       point(5,5), point(5, 8), point(-10, 5), point(0,5), point(10,5),
#       point(8,5), point(9.9,9.9)
#      )
testpoints = expand.grid(x_grid, y_grid)
testpoints_poly = apply(testpoints, 1, pointnew)
flag_insideoutside = vector(mode = "numeric", length = length(testpoints[, 1]))
for(p in 1:length(testpoints[, 1])) {
	flag_insideoutside[p] = point_in_polygon(polygon, point(testpoints[p, 1], testpoints[p, 2]))
   	}
colors_points1 = rep("darkgreen", length(testpoints[, 1]))
colors_points2 = rep("red", length(testpoints[, 1]))
colors_points = cbind(colors_points1, colors_points2)
colors_points = c("darkgreen", "red")
points(expand.grid(x_grid, y_grid)[, 1], 
       expand.grid(x_grid, y_grid)[, 2],	
       pch = 13, 
       cex = 0.9, 
       col = colors_points[flag_insideoutside + 1])

x_poly1a = c(190.9246377, 190.8951063, 190.9204987, 190.9500331, 190.9246377)
y_poly1a = c(11.5790871, 11.5310763, 11.5189782, 11.5669864, 11.5790871)       
x_poly1b = c(190.8684378, 190.8946098, 190.9238245, 190.8976493, 190.8684378)        
y_poly1b = c(11.5436762, 11.5309681, 11.5799322, 11.5926429, 11.5436762)

x_poly2a = c(190.92425, 190.946877, 190.986129, 190.963341, 190.92425) 
y_poly2a = c(11.562748, 11.545317, 11.586431, 11.604756, 11.562748)
x_poly2b = c(190.969499, 191.00876, 190.986686, 190.947399, 190.969499)
y_poly2b = c(11.528252, 11.568284, 11.585941, 11.544851, 11.528252)

x_poly3a = c(190.906144, 190.877867, 190.858626, 190.886399, 190.906144)
y_poly3a = c(11.538152, 11.546298, 11.492662, 11.485278, 11.538152)
x_poly3b = c(190.887065, 190.914021, 190.934205, 190.906852, 190.887065)
y_poly3b = c(11.485048, 11.478229, 11.53011, 11.537915, 11.485048)

x_poly4a = c(190.964119, 190.935841, 190.9166, 190.944374, 190.964119)
y_poly4a = c(11.546269, 11.554415, 11.500779, 11.493395, 11.546269)
x_poly4b = c(190.971997, 190.992181, 190.964827, 190.94504, 190.971997)
y_poly4b = c(11.486346, 11.538226, 11.546031, 11.493165, 11.486346)

x_poly5a = c(190.921006, 190.885941, 190.910086, 190.945393, 190.921006)
y_poly5a = c(11.638185, 11.592888, 11.577532, 11.621952, 11.638185)
x_poly5b = c(190.910648, 190.934191, 190.969607, 190.945993, 190.910648)
y_poly5b = c(11.577113, 11.562539, 11.605882, 11.621514, 11.577113)

x_poly6a = c(190.835584, 190.816349, 190.844127, 190.863866, 190.835584)
y_poly6a = c(11.586495, 11.532857, 11.525477, 11.578353, 11.586495)
x_poly6b = c(190.891933, 190.864574, 190.844793, 190.871754, 190.891933)
y_poly6b = c(11.570315, 11.578116, 11.525247, 11.518431, 11.570315)

glob_poly_x = c(191.0085, 190.9668, 190.9693, 190.9209, 190.8860, 190.8942, 190.8854, 190.8354,
 			    190.8162, 190.8678, 190.8587, 190.9143, 190.9225, 190.9718, 190.9924, 190.9822,
				191.0085)
glob_poly_y = c(11.56832, 11.60238, 11.60608, 11.63807, 11.59314, 11.58727, 11.57248, 11.58677,
 				11.53264, 11.51923, 11.49292, 11.47800, 11.49914, 11.48682, 11.53835, 11.54114,
				11.56832)
glob_poly_xy = cbind(glob_poly_x, glob_poly_y)

#points_polygon = apply(glob_poly_xy, 1, pointnew)
#natural_indices_polygon = seq(length(xy_vertices[, 1]))
#polygon = createPolygonNew(points_polygon)

ra_seq = seq(190.81, 191.01, length.out = 36) 
dec_seq = seq(11.48, 11.64, length.out = 36)
radec_grid = expand.grid(ra_seq, dec_seq)
#radec_grid_poly = apply(radec_grid, 1, pointnew)
#flag_insideoutside = vector(mode = "numeric", length = length(radec_grid[, 1]))
#for(p in 1:length(radec_grid[, 1])) {
#	flag_insideoutside[p] = point_in_polygon(polygon, point(radec_grid[p, 1], radec_grid[p, 2]))
#   	}

idx_within = gridknots_inside_polygon(ra_seq, dec_seq, glob_poly_xy)
#colors_points1 = rep("darkgreen", length(radec_grid[, 1]))
#colors_points2 = rep("red", length(radec_grid[, 1]))
#colors_points = cbind(colors_points1, colors_points2)
#colors_points = c("darkgreen", "red")
rand_x = runif(1000, 190.82, 191.01)
rand_y = runif(1000, 11.48, 11.64)
idx_within_rand = points_inside_polygon(rand_x, rand_y, glob_poly_xy)

# Plot.
pathpdf1 = paste("~/Desktop/footprint.pdf", sep = "")
pdf(pathpdf1, width = 7, height = 7, paper = "special")
par(family = "sans", tcl = 0.5, cex.lab = 1.2, cex.axis = 0.8, 
	mai = c(0.75, 0.75, 0.20, 0.20), mgp = c(1.8, 0.2, 0), las = 1)
plot(ra_seq, dec_seq, 
	 xlab = "Right Ascension [deg]", 
     ylab = "Declination [deg]",
     xlim = c(191.01, 190.82), 
     ylim = c(11.48, 11.64),
     type = "n")
points(radec_grid[, 1], 
       radec_grid[, 2],	
       pch = 13, 
       cex = 0.9, 
       col = "darkgreen")
points(radec_grid[idx_within, 1], 
       radec_grid[idx_within, 2],	
       pch = 20, 
       cex = 1.9, 
       col = "midnightblue")
points(rand_x, 
	   rand_y, 
	   pch = 20, 
	   cex = 1, 
	   col = "black")
points(rand_x[idx_within_rand], 
	   rand_y[idx_within_rand], 
	   pch = 20, 
	   cex = 1, 
	   col = "goldenrod")
polygon(glob_poly_x, glob_poly_y, border = "violetred", lwd = 2.5, col = NA)
# Drawing the polygons of the HST fields. 
polygon(x_poly1a, y_poly1a, border = "red", lwd = 0.4, col = NA)	 
polygon(x_poly1b, y_poly1b, border = "red", lwd = 0.4, col = NA)	 
polygon(x_poly2a, y_poly2a, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly2b, y_poly2b, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly3a, y_poly3a, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly3b, y_poly3b, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly4a, y_poly4a, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly4b, y_poly4b, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly5a, y_poly5a, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly5b, y_poly5b, border = "black", lwd = 0.4, col = NA)
polygon(x_poly6a, y_poly6a, border = "black", lwd = 0.4, col = NA)	 
polygon(x_poly6b, y_poly6b, border = "black", lwd = 0.4, col = NA)
nty = par()$yaxp
nty[3] = nty[3]*5
ntx = par()$xaxp
ntx[3] = ntx[3]*5
par(xaxp = ntx, yaxp = nty, tcl = 0.2)
axis(1, label = FALSE)
axis(2, label = FALSE)
axis(3, label = FALSE)
axis(4, label = FALSE)
dev.off()
}