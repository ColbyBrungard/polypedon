#' polypedon 
#' @description Calculates a similarity index for every pixel within an given radius buffer of each sampling point and returns the values of those cells.   
#' @param covs raster stack of environmental covariates.
#' @param pts sampling points. A two column object of class SpatialPointsDataframe. The first column must be a unique ID for each sampling location. The second column must be the class label at each sampling location.
#' @param buffer Radius in map units (usually meters) of the disk around each point that similarity will be calculated
#' @param fac numeric, can be > 1, (e.g., fac = c(2,3)). Raster layer(s) which are categorical variables. Set to NA if no factor is present
#' @param metric character string specifying the similarity metric to be used. The currently available options are "euclidean", "manhattan" and "gower" (the default). Use gower if there are categorical variables in the raster stack. See \code{daisy} from the \code{cluster} package for more details
#' @param stand logical flag: if TRUE, then the covariate values are standardized before calculating the dissimilarities.   
#' @param cores numeric, the number of cores (i.e., threads on a local workstation) to use for parallel computing.
#' @param ... passed to plyr::llply
#' @return A dataframe with columns uid, the observation, the class name, the similarity value, all covariate values, and x and y coordinates. uid is a unique id, the class name is the label of the polypedon points, and similarity is the calculated similarity value.
#' 
#' @author Colby Brungard  \email{cbrung@@nmsu.edu}
#' 
#' @importFrom doFuture registerDoFuture
#' @importFrom stats complete.cases
#' @importFrom raster extract cellFromXY xyFromCell
#' @importFrom cluster daisy
#' @importFrom plyr llply
#' @importFrom dplyr inner_join bind_rows
#' @importFrom progressr handlers
#' 
#' @export
#' 
#' @examples 
#' library(raster)
#' library(sp)
#' library(clhs)
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
#' pts$ID <- seq(1:nrow(pts))
#' pts$soilclassname <- letters[1:nrow(pts)]
#' pts <- pts[,c(6,7)]
#' 
#' synthobs <- polyped(covs=ms, pts=pts, fac = c(4,5), buffer = 50, cores = 2)
#' 

# If I want a threshold
#polyped <- function (covs, pts, buffer, fac = NA, metric = "gower", stand = FALSE, thresh = 0.95, cores=NULL, ...)  {

polyped <- function (covs, pts, buffer, fac = NA, metric = "gower", stand = FALSE, cores=NULL, ...)  {
  
  start.time <- Sys.time()
  
  doFuture::registerDoFuture()
  future::plan("multisession", workers = cores)
  
  handlers(global = TRUE)
  
  # Iterate over every point
  # add .progress='text' if not running in parallel
  res_l <- plyr::llply(1:nrow(pts), function(i) {
    
    coords <- pts[i, ]
    
    # 2. Extract all cells within 'buffer' m of the sampling points. 
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
    # rather than a threshold I believe it is better to calculate everything as it provides a quicker way to choose a threshold.
    #res_df2 <- res_df[res_df$similarity >= thresh & res_df$similarity < 1,]
    
    # Join, which will keep only the covariate values that meet the threshold criteria
    # uncomment the following if I want to integrate a threshold into the code
    #res_filt_buff <- dplyr::inner_join(res_df2, buff_data[,-1], by = 'cells')
    res_filt_buff <- dplyr::inner_join(res_df, buff_data[,-1], by = 'cells')
    
    # if there there are no similarity values, skip appending appending 
    # the results. I would rather do this in the buff_data step to 
    # avoid calculating the distance matrix, but putting it here 
    # makes the function work 
    if(nrow(res_filt_buff) > 0) {
      
      # Get the xy coordinates
      resxy <- xyFromCell(covs, res_filt_buff$cell)
      
      # Join the observed class name, the covariate values and the xy coordinates 
       resall <- data.frame(uid=res_filt_buff[, 1], pts@data[i,1:2], res_filt_buff[, -1], resxy, stringsAsFactors = TRUE, row.names = NULL)
      names(resall)[2] <- names(pts)[1]
      names(resall)[3] <- names(pts)[2]
      resall
      
    } else {                       
      
      cat(paste("Location", pts@data[i, 1], "had 0 cells >= the threshold", "\n"))
      
    }
    
  }, .progress = "progressr", .parallel = TRUE, ...)
  end.time <- Sys.time()
  print(end.time - start.time)
  
  res_s <- dplyr::bind_rows(res_l)
  res_s
  
}