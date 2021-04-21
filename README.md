# c-Q_HESS_paper

This repository contains data and code used in the paper "Hydrologic regimes drive nitrate export behavior in human-impacted watersheds" (https://hess.copernicus.org/articles/25/1333/2021/hess-25-1333-2021.html).

Gorski, G. and Zimmer, M. A.: Hydrologic regimes drive nitrate export behavior in human-impacted watersheds, Hydrol. Earth Syst. Sci., 25, 1333–1345, https://doi.org/10.5194/hess-25-1333-2021, 2021.

The repo contains three sets of data:
  L0-Raw daily discharge and nitrate data for the five watersheds downloaded from USGS through the NWIS portal
  L1-Daily data with "event-periods" and "non-event periods" delieated for later analysis
  L2-Event and non-event data grouped into data stuctures for easier manipulation in R
  
 The repo contains several codes for reproducing the analysis and plots found in the paper. The idea is that each script should be able to stand alone to regenerate plots from the paper. All plots can be regenerated with the data from this repository. The only data from the paper that is not in this repository is the spatial data used to make Figure 1. Due to its size, it is not stored here, but is available upon request. All the data used for this project were publicly available.
 
 Questions? Email: ggorski@ucsc.edu
