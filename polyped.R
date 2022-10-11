#Modified function
#covs	raster stack of environmental covariates
#pts	sampling points, object of class SpatialPointsDataframe. The first and only column must be the names of the soil class at each observation
#buffer	Radius of the disk around each point that similarity will be calculated. The smaller the buffer the larger the similarity of all cells. Standardizing does not affect this. 
#fac	numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present
#metric	character string specifying the similarity metric to be used. The currently available options are "euclidean", "manhattan" and "gower" (the default). See daisy from the cluster package for more details
#stand	logical flag: if TRUE, then the measurements in x are standardized before calculating the dissimilarities.
#thresh numeric, a threshold similarity value. All cells with similarity > than this threshold will be retained. Values with similarity = 1 will not be kept as this is often the same cell as the observation. In practice this has been found to be useful when set to 0.95.  

# The error "Error in data.frame(class.name = pts@data[i, 1], res_filt_buff[, -1], : arguments imply differing number of rows: 1, 0", means that no cells with a threshold similarity value were found so you should lower the threshold value.

polyped <- function (covs, pts, buffer, fac = NA, metric = "gower", 
                     stand = FALSE, thresh = 0.99, ...) 
{
  res_l <- plyr::llply(1:nrow(pts), function(i) {
    coords <- pts[i, ]
    buff_data <- raster::extract(x = covs, y = coords, buffer = buffer, 
                                 cellnumbers = TRUE, method = "simple", df = TRUE)
    cellnum <- cellFromXY(covs, coords)
    buff_data <- data.frame(buff_data[complete.cases(buff_data), 
    ], stringsAsFactors = TRUE)
    if (!any(is.na(fac))) {
      buff_data[, fac + 1] <- lapply(buff_data[fac + 1], 
                                     factor)
    }
    gower_dissim <- daisy(x = buff_data[, names(covs)], metric = metric, 
                          stand = stand, warnBin = FALSE)
    gower_dissim <- cbind(buff_data$cell, as.matrix(gower_dissim))
    gower_sim <- 1 - gower_dissim[gower_dissim[, 1] == cellnum, ]
    
    res_df <- data.frame(cells = buff_data$cell, similarity = gower_sim[-1], stringsAsFactors = TRUE)
    
    # filter for all that are >= a similarity threshold
    res_df2 <- res_df[res_df$similarity >= thresh & res_df$similarity < 1,]
    # Join, which will keep only the covariate values that meet the threshold criteria
    res_filt_buff <- dplyr::inner_join(res_df2, buff_data[,-1], by = 'cells')
    # Get the xy coordinates
    resxy <- xyFromCell(covs, res_filt_buff$cell)
    # Join the xy to the coordinates
    resall <- data.frame(class.name = pts@data[i,1], res_filt_buff[, -1], resxy, stringsAsFactors = TRUE)
    resall
  }, ...)
  
  res_s <- dplyr::bind_rows(res_l)
  res_s
}