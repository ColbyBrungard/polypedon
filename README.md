# polypedon

A parallelized function to calculate a similarity index for every pixel within an given radius buffer of each sampling point and returns the values of those cells. This is intended to supplement training data for digital soil class mapping based on the polypedon idea. 

Thoughts on covariates:
- Initial testing suggests that all covariates which will be used for prediction should be included in this function. Choosing only two or three covariates, then using the resulting polypedon locations to exctract covariates values to be used in the modeling results in an accuracy decrease and a duplication of effort. 

Observations of buffer size:
- Smaller buffers generally result in lower similarity values.
- The buffer size should generally be chosen based on soil spatial variation. I have generally found the buffers between 20 and 100 m are sufficient, but this also depends on the grid size of input variables. I recomend that the users choose a few observations then compare multiple buffer sizes by plotting these in a gis. It has been my observation that this is the easiset way to determine an appropriate buffer size. 
- I generally favor more conservative (smaller) buffers because it is easy to very quickly generate a (probably too) large number of observations. 

Observations on chosing threshold values:
- This function calucates a similarity value between the central point and the surrounding raster cells. The user must then choose a similarity threshold above which the cells should be similar enough to the observation to be confidently labled as the same class. This is subjective and should be chosen based on knowledge of the area. 

Misc notes
- The number of observations is more influential than is the size of the radius (e.g., 50 vs 20) on the speed of the algorithm. 
- This function does NOT retain the center cell where the observation was located. 
