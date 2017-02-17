# Gps Integration
Model for integrating GPS tracks and traditional survey data, based on methods from *Talluto et al. 2016. Glob. Ecol. Biogeogr.*

## Data
*Data are not included in repo; please see maintainers for access to data.*

`dat/shw_scale.Rdata`

R data file for Mediterranean shearwater data. Column descriptions:

* `count` = aerial survey counts
* `nlocs` = number of GPS locations (0 are cells within the range of the birds, which where not visited, I used the pseudo tracks for the moment but could use all cells with a specific radius later on).
* `ninds` = number of unique individuals visiting the cell
* `nvisits` = number of visits in a cell
* environmental variables are scaled.

I removed the count data from May and June to avoid temporal mismatch which could mess up the environment relationships since tracks data were collected in August. The count data are then from July and August, and the environmental variables are monthly averages from August.
