# polypedon

A parallelized function to calculate a similarity index for every pixel within an given radius buffer of each sampling point and returns the values of those cells that meet or exceed a user specified threshold. This is intended to supplement training data for digital soil class mapping based on the polypedon idea. 

Thoughts on covariates:
- It generally apperas that choosing fewer covariates is better than many covariates. Covariates should be calculated over local neighborhoods because the intention of this function is to find surrounding areas that are very similar to a location where there is a real observation.  

Observations of buffer size:
- Smaller buffers generally result in lower similarity values.
- The buffer size should generally be chosen based on soil spatial variation. I have generally found the buffers between 20 and 100 m are sufficient, but this also depends on the grid size of input variables. I recomend that the users choose a few observations then compare multiple buffer sizes by plotting these in a gis. It has been my observation that this is the easiset way to determine an appropriate buffer size. 
- I generally favor more conservative (smaller) buffers because it is easy to very quickly generate a (probably too) large number of observations. 

Observations on threshold values:
- I find that it is most beneficial to adjust the radius first, then modify the threshold. I find that a threshold of 0.9 seems to consistently work in both arid and recently glaciated humid environments. 

Notes on speed:
- The number of observations is more influential than is the size of the radius (e.g., 50 vs 20). 
