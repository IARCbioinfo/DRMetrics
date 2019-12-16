# Dimensionality Reduction Metrics

This repository contains scripts to compute several metrics that evaluate dimensionality reduction techniques:

-   [Sequence difference view](#sequence-difference-view)
    -   [Sequence difference view calculation
        `Seq_calcul`](#sequence-difference-view-calculation-seq_calcul)
    -   [Main function for sequence difference view
        `Seq_main`](#main-function-for-sequence-difference-view-seq_main)
    -   [Graphic of means of sequence difference values by k values
        `Seq_graph_by_k`](#graphic-of-means-of-sequence-difference-values-by-k-values-seq_graph_by_k)
    -   [Sequence difference values permutation test
        `seq_permutation_test`](#sequence-difference-values-permutation-test-seq_permutation_test)
    -   [Sequence difference map
        `SD_map_f`](#sequence-difference-map-sd_map_f)
-   [Spatial autocorrelation](#spatial-autocorrelation)
    -   [Moran index main function
        `moran_I_main`](#moran-index-main-function-moran_i_main)
    -   [Calcul of Moran Indexes for high dimensional data
        `moran_index_HD`](#calcul-of-moran-indexes-for-high-dimensional-data-moran_index_hd)
    -   [Moran significance test for high dimensional data
        `moran_stat_HD`](#moran-significance-test-for-high-dimensional-data-moran_stat_hd)
    -   [Graphic of Moran Indexes for each variable and each method
        `moran_I_scatter_plot_by_k`](#graphic-of-moran-indexes-for-each-variable-and-each-method-moran_i_scatter_plot_by_k)

Sequence difference view
========================

Sequence difference view calculation `Seq_calcul`
-------------------------------------------------

### Description

This function allows to calculate the sequence diffeence (*S**D*) view
metrics. For a point *i* the formula is :
\(D\_k(i) = \\frac{1}{2} \\sum\_{j \\in V^l\_k(i)}\[k-\\rho^l\_i(j)\].\|\\rho^l\_i(j)-\\rho^h\_i(j)\|+ \\frac{1}{2} \\sum\_{j \\in V^h\_k(i)}\[k-\\rho^h\_i(j)\].\|\\rho^l\_i(j)-\\rho^h\_i(j)\|, \\label{EqSD}\)
 where *V*<sub>*k*</sub><sup>*d*</sup>(*i*) is the *k*−neighborhood of
*i* in the dimension *d*, and *ρ*<sub>*i*</sub><sup>*d*</sup>(*j*) is
the rank of *j* in the *k*−neighborhood of *i*

### Usage

`Seq_calcul <- function(l_data, dataRef, listK)`

### Arguments

-   **l\_data** : list of data frame whose structure is :

| Sample\_ID | x        | y         | …   |
|------------|----------|-----------|-----|
| ID1        | x\_coord | y\_coords | …   |

These data frames contain samples’ coordinates which could be defined in
ℝ<sup>𝕟</sup>. **warning : ** It must be a list of dates frames and not
a list of data tables.

-   **dataRef** : reference data frame whose structure is defined
    above.  
    **warning : ** It must be a list of dates frames and not a list of
    data tables.

-   **listK** : list k levels

### Details

-   A inner join on samples’ ID is effected if they differs between the
    different data frames.
-   Calculations use a parallel computing according the levels *k*.

### Value

A list of containing *l* elements is returned (where *l* corresponds to
the number of data frames containing in `l_data`). Each element contains
*n* SQ values, where *n* is the number of common samles’ ID between the
reference data frame and the data frames in `l_data`.

Main function for sequence difference view `Seq_main`
-----------------------------------------------------

### Description

This function allows to calculate the *S**D* values for several data
frames and for differents *k* levels. Distributions of means *S**D*
values by levels *k*, *i.e* $\\overline{SD}\_k$ could be plot. Finally
statistic tests could be computed if at least two low dimensionals
projections are given in input. \#\#\# Usage

`Seq_main <- function(l_data, dataRef, listK, colnames_res_df = NULL , filename = NULL , graphics = FALSE, stats = FALSE)`

### Arguments

-   **l\_data** : list of data frames whose the respective structure
    must be :

| Sample\_ID | x        | y         | …   |
|------------|----------|-----------|-----|
| ID1        | x\_coord | y\_coords | …   |

These data frames contain samples’ coordinates which could be defined in
ℝ<sup>𝕟</sup>.

-   **dataref** : data frame of reference whose strucuture is the same
    as define above.

-   **listK** : list *k* levels.

-   **colnames\_res\_df** : This optional argument allows to specify the
    colnames of the returned data frame and also the plot’s legend if it
    was computed. If this argument is unsecified then the default values
    will be set to : *V*1, *V*2, ..., *V**n* (where *n* is the length of
    `l_data`).

-   **filename** : This optional arguement allows to defined the
    filename on which results will be written. If this argument is
    unspecified then results will be returned and not written. If users
    choose a filename that ever exits in the current directory a
    incrementation in the filename will be done.

-   **graphics** : This boolean argument allows to computes plot. This
    plot will represent means of *S**D* values for the different *k*
    levels and for the different data frames in `l_data`.

-   **stats** : This option allows to run statistic tests, it is
    available only if the number of defines method is higher at least
    equal to two, (*i.e* `l_data`’s length is  ≥ 2). If only two data
    frames were given as input via the `l_data` then a test will be
    computed to compare the distribution of the the means by k levels of
    absolute differences between *S**D* values. If more than two methods
    were defined then paired tests are done. If more than 30
    $\\overline{SD}\_k$ values have been computed Student tests are
    done, otherwise Wilcoxon tests are preferred.

### Details

-   A inner join on samples’ ID is effected if those differs between the
    different data frames.

### Value

According options activated the return list contains the following
elements :

-   **Seq\_df** : data frame containing a column with the samples’ Id, a
    column correspoding to the levels *k*, and *n* colunms corresponding
    to the *S**D* values. This data frame could be written in a file if
    `filename` is defined.

-   **Seq\_mean\_by\_k** : data frame containing the
    $\\overline{SD}\_k$, for each data frame contained in `l_data`.

-   **paired\_test** : Student or Wilcoxon paired test’s results,
-   **pairwise\_tests** : Matrix of Student or Wilcoxon pairwise tests’
    p.value.

-   **graphics** : GGplot of the $\\overline{SD}\_k$ in function of the
    levels *k* if the graphic option was activated.

Graphic of means of sequence difference values by k values `Seq_graph_by_k`
---------------------------------------------------------------------------

### Description

This function displays the graphic of means of sequence difference
values by k values.

### Usage

`Seq_graph_by_k  <-function (data_Seq, Names=NULL, data_diff_mean_K = NULL, log=False)`

### Arguments

-   **data\_Seq** : data frame of sequence difference values structured
    such as :

| Sample\_ID | K           | CP1      |
|------------|-------------|----------|
| ID1        | 1Kst\_level | CP1\_id1 |

-   **Names** : optional argument allowing to precise legend labels. If
    this argument is unprecised lengend labels are equal to `data_Seq`’s
    colnames.

-   **data\_diff\_mean\_K** : optional data frame contining means of SQ
    values by K level. If this argument is precised then means are not
    calculated.

-   **log**:Boolean optional argument, if it set to true then a
    logarithmic scale will be used.

### Value

A ggplot object is returned.

### See also

`Seq_main`

Sequence difference values permutation test `seq_permutation_test`
------------------------------------------------------------------

### Description

Then this function test the random hypothesis *i.e*.: Does *S**D* values
calculated on real data set are equivalent to those expected on random
data ? In order to do this *n* simulations are realized. According these
simulations the $\overline{SD}\_k$ are calculated. Finally wilcoxon
test is effected to compare the mean random distribution and the real
one.

### Usage

`seq_permutation_test <- function(data, data_ref, list_K, n=30, graph = TRUE)`

### Arguments

-   **data** : data frame defined such as :

| Sample\_ID | x        | y         | …   |
|------------|----------|-----------|-----|
| ID1        | x\_coord | y\_coords | …   |

-   **data\_ref** : reference data frame whose structure is equivalent
    to the one defined above.

-   **listK** : list k levels.

-   **n** : number of simulations.

-   **graph** : optional boolean argument, if this argument is TRUE,
    simulations resulting graphic is computed.

### Value

This function returns the Wilcoxon test’s results.

### Details

According the `n` the simulation could be long.

Sequence difference map `SD_map_f`
----------------------------------

### Description

This function allows display the sample’s mean $\\overline{SD}\_k$ on a
two dimensional projection.

### Usage

`SD_map_f <- function(SD_df, Coords_df, legend_pos = "right")` \#\#\#
Arguments

-   **SD\_df** : data frame defined such as :

| Sample\_ID | k       | SD      | …   |
|------------|---------|---------|-----|
| ID1        | levelk1 | SD\_1,k | …   |

-   **Coords\_df** : data frame defined such as : \| Sample\_ID \| X \|
    Y \| … \| \|———–\|———\|———-\|—–\| \| ID1 \| x1 \| y1 \| … \|

-   **legend\_pos**: Optional argument to define legend’s position

### Value

This return a map of the ean $\\overline{SD}\_k$ per sample.

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

Moran significance test for high dimensional data `moran_stat_HD`
-----------------------------------------------------------------

### Description

Singnificance test are computed according Monte Carlo procedure. Like
this *n* simulations are done, at each iteration the vector of the
variable of interest is shuffle, and then Moran Indexes are clculatated
using `moran_index_HD` function. Finally the rank of the observed Moran
Index in the resulting vector is computed to infer the p.value. This
p.value is the proportion of Moran Indexes obtained with random data
that are greater then the observed Moran Index.

### Usage

`moran_stat_HD <- function(data, K, spatial_att, obs_moran_I, nsim = 99)`

### Arguments

-   **data** : data frame defining such as :

| Sample\_ID | x         | y         | z         | …   |
|------------|-----------|-----------|-----------|-----|
| MYID       | x\_coords | y\_coords | z\_coords | …   |

-   **K** : numeric argument defining a k level.

-   **spatial\_att** : data frame which contains variables values.

| Sample\_ID | Variable1 |
|------------|-----------|
| MYID       | V1\_myid  |

-   **obs\_moran\_I** : observed moran Index computed according the real
    spatial distribution of the variable.

-   **nsim** : number of simulations.

### Value

Moran Significance test p.value.

### See also

`moran_I_main`

Graphic of Moran Indexes for each variable and each method `moran_I_scatter_plot_by_k`
--------------------------------------------------------------------------------------

### Description

This function allows to displays the plot of Moran Index values for each
variable and and each method,either a scatter or a boxplot is display if
the Moran Index values have been calculated for several k levels.

### Usage

`moran_I_scatter_plot <- function(data, Xlab = NULL, Ylab=NULL, Title= NULL)`

### Arguments

-   **data** : 3D array containing Moran Index values whose the
    structure is the following :

k = i

|         | Varaible1  | Varaible2  | Varaible3  | …   |
|---------|------------|------------|------------|-----|
| Method1 | MI\_v1\_m1 | MI\_v2\_m1 | MI\_v3\_m1 | …   |
| Method2 | MI\_v1\_m2 | …          | …          | …   |
| …       | …          | …          | …          | …   |

------------------------------------------------------------------------

k = j

|         | Varaible1  | Varaible2  | Varaible3  | …   |
|---------|------------|------------|------------|-----|
| Method1 | MI\_v1\_m1 | MI\_v2\_m1 | MI\_v3\_m1 | …   |
| Method2 | MI\_v1\_m2 | …          | …          | …   |
| …       | …          | …          | …          | …   |

-   **Xlab** : this optional argument is used to define the x-axis
    label.

-   **Ylab** : this optional argument is used to define the y-axis
    label.

-   **Title** : This optional argument is used to define the plot title.

### Value

This function return a ggplot object.

### See also

`moran_I_main`
