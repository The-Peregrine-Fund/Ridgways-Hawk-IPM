# Ridgway's Hawk Integrated Population Model (IPM) and Population Viability Analysis (PVA)
Supplemental materials for: Rolek, B.W., McClure, CJW, Dunn, L., Curti, M., ... Ridgway's Hawk IPM and PVA

Contact: Brian Rolek  
Email: rolek.brian@peregrinefund.org

Metadata, data, and scripts used in analyses can be found at <https://github.com/The-Peregrine-Fund/XXXXX>.

The full workflow can be accessed as an html website at:
<https://the-peregrine-fund.github.io/XXXXX/>.

A permanent archive and DOI is available at: https://zenodo.org/doi/XXXXX

# Metadata and Data Dictionary
All data are contained in the file "data/data.Rdata" which is readable in R.
Once loaded into R, for example load("data/data.Rdata"), the data are included in two lists to be read into NIMBLE for Bayesian analysis. These lists are saved as the objects "constl" and "datl". Other objects loaded are related to initial values.

## Metadata for Data (object "datl")
* prod- long format data of the number of successfully fledged young by each pair each year.
* y- Mark-recapture- A matrix of dimensions nind by nyr. Detections of individuals and their life stage. Coded as 1 = seen first-year fledgling, 2 = seen nonbreeder adult, 3 = seen breeder adult, and 4 = not seen.
* countsAdults- A matrix of dimensions nyr by nsite. The number of breeder adults seen, divided by two and rounded.
* countsFY- A matrix of dimensions nyr by nsite. The number of first-year fledglings seen, divided by two and rounded.
* mu.zeroes and mu.zeroes2 - A vector of length p. Zeroes input to specify zero-centered means for random effects from site and time.
* constraint_data- A matrix of dimensions nyr by nsites. An index that indicated which data to constrain to be > zero. This constraint prevents negative numbers when hacking occurs. 

## Metadata for Constants (object "constl")
* nind- The number of individuals included in the mark-recapture-resight dataset.
* nyr- The duration of datasets as the number of years.
* nsite- The number of sites included in analysis.
* first- The first year of capture for each individual in the mark-recapture-resight dataset.
* site- a nind x nyr matrix specifying the site an individual occupied (site=1 is Los Haitises, site=2 is Punta Cana)
* yrind.surv- A matrix to help calculate annual averages of survival, recruitment, and detection for the population submodel. It's necessary because NIMBLE can't loop over non-sequntial indices otherwise. 
* surv.end- An index indicating the end for yrind.surv. 
* treat.pair- Whether or not a pair's nest was treated. treated=1, not treated=0.
* npairobs- Total number of observed pairs attempting breeding summed over all years.
* year.pair- Index of year. year.pair=1 is 2011, year.pair=2 is 2012, and so forth. 
* yrind.pair- A matrix to calculate annual averages of fecundity for the population submodel. It's necessary because NIMBLE can't loop over non-sequntial indices otherwise.
* pair.end- An index indicating the end for yrind.pair.
* site.pair- An index indicating site for fecundity submodel. site.pair=1 is Los Haitises and site.pair=2 is Punta Cana.
* pPrior- priors for abundance of each population segment when t=1.
* p and p2- the number of random effects for multivariate normal distributions.
* s.end- The end of each vector for priors for pPrior, that is the priors for abundance when t=1.
* hacked- A matrix of dimensions nind x nyr. A categorical covariate indicating whether an individual was hacked (hacked=1) or not (hacked=0).
* hacked.counts- The number of translocated individuals from (negative values) and to (positive values) each site.
* effort2- Number of surveyor days, centered on the mean, and scaled. Punta Cana was considered a census so we imputed zeroes for that site. We also imputed the mean survey effort (effort=0) for early years at Los Haitises because these count data were removed and appeared to be biased.
