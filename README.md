

# Dimensionality Reduction Metrics

This repository contains R functions to evaluate the quality of projections obtained after using dimensionality reduction techniques. A nextjournal notebook is associated to this repository and uses the functions described in this README file to evaluate the quality of a molecular map of lung neuroendocrine tumors produced using the UMAP algorithm.

Sequence difference view (SD) metric
========================

SD metric calculation for one sample `compute_SD`
-------------------------------------------------

### Description

This function computes the sequence difference (*SD*) view metric value for a single given sample (*i*), following the equation 3 described by [Martins *et al.* in 2015](https://bdpi.usp.br/bitstream/handle/BDPI/49441/2722542.pdf;sequence=1). This dissimilarity metric compares the *k*-neighborhood of a given sample in two different dimensional spaces. The lower is the SD value, the better is the neighborhood preservation.

### Usage

`compute_SD(dist_space1,dist_space2,k)`

### Arguments
  
-   **dist_space1**: vector containing the distances of sample *i* to all samples in space1

-   **dist_space2**: vector containing the distances of sample *i* to all samples in space2

-   **k**: number of neighbors considered

### Value

A numeric value corresponding to the *SD* value is returned.

SD metric calculation for all samples `compute_SD_allSamples`
-----------------------------------------------------

### Description

This function computes the *SD* metric for all samples included in the dimensionality reduction. The metric is computed to compare one or multiple comparison reduced spaces to a the reference space. The *SD* values are computed for several *k* values (number of neighbors to consider).

### Usage
`compute_SD_allSamples(distRef,List_projection,k_values,colnames_res_df, threads=2)`

### Arguments 

-   **distRef**: vector containing the distances of sample *i* to all samples in the reference space
    
-   **List_projection**: list of data frames where each data frame contains the coordinates of all samples in each reduced space for which the SD metric needs to be calculated.
    
-   **k_values**: vector listing the *k* values corresponding to the number of neighbors considered

-   **colnames_res_df**: vector specifying the colnames associated to the computed SD values in the returned data frame. The vector should have the same length as *List_projection*


### Value

-   Data frame containing a column with the samples IDs, a column correspoding to the *k* values, and *n* colunms containing
    the *SD* values, *n* corresponding to the number of data frames listed in *List_projection*. 

Visualizing the SD metric in a two dimensional map `SD_map_f`
----------------------------------

### Description

This function allows to display, on a two dimensional projection, the samples *SD* values averaged over different values of *k* (number of neighbors considered to compute the *SD* metric).

### Usage

`SD_map_f(SD_df, Coords_df, legend_pos = "right")` 

### Arguments

-   **SD_df**: a data frame resulting from the call to the function *compute_SD_allSamples*. The data frame contains the following columns: i) the samples IDs, ii) *k* values, the number of neighbors considered to compute the *SD* metric, and iii) the *SD* values 

-   **Coords_df**: data frame containing the coordinates of each sample in the projection to use for the representation of the samples

-   **legend_pos**: Optional argument to define the position of the legend 

### Value

A list containing:
-   A data frame containing the same columns as *Coords_df* and a column corresponding to the averaged *SD* values over *k*.
-   The plot representing all samples in a two dimensional space. A color gradient is used to represent the *SD* values averaged over the *k* levels.

Spatial autocorrelation
=======================

Moran's Index (MI) computation `moran_I_knn`
----------------------------------------

### Description

This function allows to compute the Moran’s Index autocorrelation coefficient for a given feature used in the dimensionality reduction technique, for different levels of the parameter *k* which corresponds to the number of samples to consider for the samples neighborhood definition. The MI values are computed using the *Moran.I* function from the R package *ape*.

### Usage

`moran_I_knn(expr_data , spatial_data, listK)`

### Arguments

-   **expr_data**: matrix containing, for each sample (in rows), the values of the features (in columns) for which the *MI* values will be calculated
  
-   **spatial_data**: matrix containing the coordinates of each sample in the projection used to define the samples neighborhood 

-   **listK**: vector listing the *k* values corresponding to the number of samples considered to define samples neighborhood


### Value

-   **MI_array**: 3D array containing the MI values and their associated *p*-values for each feature (in columns), and each *k* level (in rows).





