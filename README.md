

# Dimensionality Reduction Metrics

This repository contains scripts to compute several metrics that evaluate dimensionality reduction techniques:

-   [Sequence difference view](#sequence-difference-view)
    -   [Sequence difference view calculation
        `compute_SD`](#sequence-difference-view-calculation-compute_SD)
    -   [Main function for sequence difference view
        `Seq_main`](#main-function-for-sequence-difference-view-seq_main)
    -   [Sequence difference map
        `SD_map_f`](#sequence-difference-map-sd_map_f)
-   [Spatial autocorrelation](#spatial-autocorrelation)
    -   [Moran index main function
        `moran_I_main`](#moran-index-main-function-moran_i_main)
    -   [Calcul of Moran Indexes for high dimensional data
        `moran_index_HD`](#calcul-of-moran-indexes-for-high-dimensional-data-moran_index_hd)

Sequence difference view
========================

Sequence difference view calculation `compute_SD`
-------------------------------------------------

### Description

This function computes the sequence diffeence (*SD*) view metric values for a single given point. For a point *i* the formula is :
\(D\_k(i) = \\frac{1}{2} \\sum\_{j \\in V^l\_k(i)}\[k-\\rho^l\_i(j)\].\|\\rho^l\_i(j)-\\rho^h\_i(j)\|+ \\frac{1}{2} \\sum\_{j \\in V^h\_k(i)}\[k-\\rho^h\_i(j)\].\|\\rho^l\_i(j)-\\rho^h\_i(j)\|, \\label{EqSD}\)
 where *V*<sub>*k*</sub><sup>*d*</sup>(*i*) is the *k*−neighborhood of
*i* in the dimension *d*, and *ρ*<sub>*i*</sub><sup>*d*</sup>(*j*) is
the rank of *j* in the *k*−neighborhood of *i*

### Usage

`compute_SD(dist_space1,dist_space2,k)`

### Arguments
  
-   **dist_space1** : vector containing the distances of sample *i* to all samples in space1

-   **dist_space2** : vector containing the distances of sample *i* to all samples in space2

-   **k** : number of neighbors considered

### Value

A numeric value corresponding to the *SD* value is returned.

Main function for sequence difference view `compute_SD_allSamples`
-----------------------------------------------------

### Description

This function computes for all samples the *SD* metrics comparing each multiple chosen reduced space to a the reference space. The *SD* values are computed for several different *k* values (number of neighbors to consider).

### Usage
`compute_SD_allSamples(distRef,List_projection,k_values,colnames_res_df, threads=2)`

### Arguments 

-   **distRef** : vector containing the distances of sample *i* to all samples in the reference space
    
-   **List_projection** : list of data frames where each data frame contains the coordinates of all samples in each reduced space for which the SD metric needs to be calculated.
    
-   **k_values** : vector listing the k values corresponding to the number of neighbors considered

-   **colnames_res_df** : vector specifying the colnames of the returned data frame. 


### Value

-   Data frame containing a column with the samples IDs, a column correspoding to the *k* values, and *n* colunms containing
    the *SD* values, *n* corresponding to the number of data frames listed in *List_projection*. 

Sequence difference map `SD_map_f`
----------------------------------

### Description

This function allows display the samples mean $\\overline{SD}\_k$ values on a two dimensional projection.

### Usage

`SD_map_f(SD_df, Coords_df, legend_pos = "right")` 

### Arguments

-   **SD\_df** : data frame defined such as :

| Sample\_ID | k       | SD      | 
|------------|---------|---------|
| ID1        | levelk1 | SD\_1,k |

-   **Coords_df** : data frame containing the coordinates of each sample in the projection to use for the representation of the samples

-   **legend\_pos**: Optional argument to define the position of the legend 

### Value

This return a map of the mean $\\overline{SD}\_k$ per sample.

Spatial autocorrelation
=======================

Moran index main function `moran_I_main`
----------------------------------------

### Description

This function allows to calculate Moran’s Index, spatial autocorrelation
index such as:
$$ I = \\frac{N \\sum\_{i=1}^N \\sum\_{j=1}^N W\_{ij}(x\_i - \\bar{x})(x\_j - \\bar{x})}{\\sum\_{i=1}^N \\sum\_{j=1}^N (W\_{ij}) \\sum\_{i=1}^N (x\_i - \\bar{x})^2}$$

where *W* is a binary spatial weight matrix, defining through the
*K*−nearest neighbors method (KNN) such as *W*<sub>*i**j*</sub> equals
one if *i* belongs to the first *k* neighbors of *j*, and zero
otherwise, and where *x* is the value of the variable associated to the
sample *i*, and reciprocally for *x**j*, and *x̄* corresponds to the
general mean of *x*. The results values are calculated for several
variables according several projection and for differents k levels.
Graphics of Moran Indexes distribution for each variable, could be
computed. Finally significance tests according the Monte Carlo
procedure, could be computed. \#\#\# Usage

`moran_I_main <-function(l_coords_data , spatial_att, listK, nsim = 500, Stat=FALSE, Graph = FALSE, methods_name = NULL),`

### Arguments

-   **l\_coords\_data** : list of coordinates data frames whose
    structure is :

| Sample\_ID | x        | y         | …   |
|------------|----------|-----------|-----|
| ID1        | x\_coord | y\_coords | …   |

These data frames contain samples’ coordinates which could be defined in
ℝ<sup>𝕟</sup>.

-   **spatial\_att** : data frame containing variables values.

| Sample\_ID | Variable1 | Variable2 | …   |
|------------|-----------|-----------|-----|
| ID1        | V1\_id1   | V2\_id1   | …   |

-   **listK** : list k values

-   **nsim** : number of simulations for the significance test.

-   **Stat** : optional boolean argument, if this argument is set to
    TRUE, then the significance test will be calculated.

-   **Graph** : optional boolean argument, if this argument is set to
    TRUE, then the graphic of Moran Index distributions is drawn. This
    graphic depicts Moran Indexes distributions for each variable.

-   **methods\_name** : optional parameter allowing to specify the named
    of the space included in `l_coords_data` argument.

### Details

A inner join on samples’ ID is effected if those differs between the
different data frames. Moran Indexes and statistics are computed
according `moran_index_HD` and `moran_stat_HD`functions.

### Value

According options activated the return list contains the following
elements :

-   **MI\_array** : 3D array containing Moran Index for each projection
    in row *i*, each variable in colunm *j* and each *k* level.

-   **MS\_array** : 3D array containing Moran Significance tests’
    p.value for each projection in row *i*, each variable in colunm *j*
    and each *k* level.

-   **Graph** : GGplots are printed if the option is activated. The
    plots correspond to the Moran’s index values for each variables in
    function of the k levels, for each spaces.

### See also

`moran_index_HD`, `moran_stat_HD` and `moran_I_scatter_plot`

Calcul of Moran Indexes for high dimensional data `moran_index_HD`
------------------------------------------------------------------

### Description

This function allows to calculate Moran Indexes for high dimensional
data by generalizing the process effected in 2D. In order to get Moran
Indexes the *k*−nearest neighbors are defined for each sample according
the brute method of knn algorithm. This *k*−nearest neighbors is use to
define the spatial weights matrix. Then Moran Indexes are computed
classically.

### Usage

`moran_index_HD <- function(data, spatial_att, K, merge = TRUE)`

### Arguments

-   **data** : data frame defining such as :

| Sample\_ID | x         | y         | z         | …   |
|------------|-----------|-----------|-----------|-----|
| MYID       | x\_coords | y\_coords | z\_coords | …   |

-   **spatial\_att** : data frame which contains variables values.

| Sample\_ID | Variable1 | Variable2 | …   |
|------------|-----------|-----------|-----|
| MYID       | V1\_myid  | V2\_myid  | …   |

-   **K** : numeric argument defining k level.

-   **merge** : optional boolean argument that allows to checked if
    `spatial_att` and `datta` contains the same samples’ ID. If samples’
    ID differs then an inner join will be done.

### Value

Moran Index (numeric value).


