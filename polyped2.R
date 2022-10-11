#' polypedon 
#' @description Calculates Gower's similarity index for every pixel within an given radius buffer of each sampling point and returns the values of those cells that meet or exceed a user specified threshold
#' @param covs raster stack of environmental covariates
#' @param pts sampling points, object of class SpatialPointsDataframe. The first and only column must be the names of the soil class at each sampling location.
#' @param buffer Radius in meters of the disk around each point that similarity will be calculated
#' @param fac numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present
#' @param metric character string specifying the similarity metric to be used. The currently available options are "euclidean", "manhattan" and "gower" (the default).  See \code{daisy} from the \code{cluster} package for more details
#' @param stand logical flag: if TRUE, then the measurements in x are standardized before calculating the dissimilarities. 
#' @param thresh numeric, a threshold similarity value. All cells with similarity > than this threshold will be retained. Values with similarity = 1 will not be kept as this is often the same cell as the observation. In practice this has been found to be useful when set to > 0.95.
#' @param ... passed to plyr::llply
#' @return A dataframe with columns class.name, similarity value, covariate values, and x and y coordinates.
#' 
#' @author Colby Brungard  \email{cbrung@@nmsu.edu}
#' 
#' 
#' @importFrom stats complete.cases
#' @importFrom raster extract cellFromXY xyFromCell
#' @importFrom cluster daisy
#' @importFrom plyr llply
#' @importFrom dplyr inner_join bind_rows
#' 
#' @export
#' 
#' @examples 
#' library(raster)
#' library(sp)
#' 
#' data(meuse.grid)
#' coordinates(meuse.grid) = ~x+y
#' proj4string(meuse.grid) <- CRS("+init=epsg:28992")
#' gridded(meuse.grid) = TRUE
#' ms <- stack(meuse.grid)
#' 
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(1)
#' pts <- clhs(ms, size = 3, iter = 100, progress = FALSE, simple = FALSE)
#' 
#' pts <- pts$sampled_data 
#' pts$soilclassname <- c('A', 'B', 'C')
#' pts <- pts[,6]
#' 
#' synthobs <- polyped(covs=ms, pts=pts, fac = c(4,5), buffer = 50, thresh = 0.8)
#' 
similarity_buffer <- function(covs, pts, buffer, fac = NA, metric = "gower", stand = FALSE, ...) {
  
  # Iterate over every point
  res_l <- plyr::llply(1:nrow(pts), function(i) {
    
    coords <- pts[i, ]
    
    # 2. Extract all cells within x m of the sampling points. 
    buff_data <- raster::extract(
      x = covs, 
      y = coords, 
      buffer = buffer, 
      cellnumbers = TRUE, 
      method = 'simple', 
      df = TRUE
    )
    
    # 3. Apply Gower's similarity index to each element of list of extracted raster values
    
    # Get the cell numbers from each sample point to identity the right column in the similarity matrix. 
    cellnum <- cellFromXY(covs, coords)
    
    # 3.b Calculate Gower's similarity index around each point. 
    #   I used Gower's because it can handle categorical covariates, 
    #   but there could be other options. 
    
    # Only retain cases without NA values
    buff_data <- data.frame(buff_data[complete.cases(buff_data), ], stringsAsFactors = TRUE)
    
    # If there are some factor data
    if (!any(is.na(fac))) {
      buff_data[, fac + 1] <- lapply(buff_data[fac + 1], factor)
    }
    
    # Calculate gowers similarity index 
    gower_dissim <- daisy(x = buff_data[, names(covs)], metric = metric, stand = stand, warnBin = FALSE)
    # turn dissimilarity object to matrix
    gower_dissim <- cbind(buff_data$cell, as.matrix(gower_dissim)) 
    
    # Select the row of similarity indices with cell number equal to the cell number of the 
    # sample point and convert dissimilarity to similarity by subtracting from 1.  
    gower_sim <- 1 - gower_dissim[gower_dissim[, 1] == cellnum, ] 
    
    # Combine the cellnumbers of the raster to the similarity index. 
    res_df <- data.frame(cells = buff_data$cell, similarity = gower_sim[-1], stringsAsFactors = TRUE)
    
    # filter for all that are >= a similarity threshold
    res_df2 <- res_df[res_df$similarity >= thresh & res_df$similarity < 1,]
    
    # Join, which will keep only the covariate values that meet the threshold criteria
    res_filt_buff <- dplyr::inner_join(res_df2, buff_data[,-1], by = 'cells')
    
    # Get the xy coordinates
    resxy <- xyFromCell(covs, res_filt_buff$cell)
    
    # Join the observed class name, the covariate values and the xy coordinates 
    resall <- data.frame(class.name = pts@data[i,1], res_filt_buff[, -1], resxy, stringsAsFactors = TRUE)
    resall
  }, ...)
  
  res_s <- dplyr::bind_rows(res_l)
  res_s

} 