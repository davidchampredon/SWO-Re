# SWO-Re
**Estimates of the effective reproduction number for South-West Ontario**

***--WARNING-- this is a work in progress!***


This repo host R scripts to estimate the effective reproduction number ("Re") of reported COVID-19 cases for multiple Public Health Units (PHUs) in Ontario. Currently, there is a focus on the South-West regions. The estimation uses the R package `EpiEstim`. 

## Files Description

 * `swo-Re.R`: This is the main script that retrieves and filter epidemiological data, estimates Re and plot its evolution over time.
 * `utils.R` : Implements utility functions.
 * `phus.csv`: Name of the PHUs selected. The names must be consistent with the data source (which is the [line list from MoH](https://data.ontario.ca/dataset/confirmed-positive-cases-of-covid-19-in-ontario/resource/455fd63b-603d-4608-8216-7d8647f43350)).
 * `phu-group.csv`: Grouping of multiple PHUs (reported cases are aggregated within the group).



