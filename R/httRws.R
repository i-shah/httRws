#' httRws: A package to access the EPA/ORD/NCCT high-throughput transcriptomics (HTTr) data.  
#'
#' The httRws package provides functions to search for chemicals that have been run in the HTTr 
#' platform and to retreive the raw and processed data on a probe / gene level. These functions
#' utilize the NCCT internal (development) RESTful webservices for accessing HTTr data on 
#' ~2,200 chemicals, 8 concentrations, 1 cell type (MCF7), ~58,000 samples, ~21,000 probes, ...
#' Currently, these webservices are only available on the EPA intranet at (httr-dev.epa.gov) and are
#' for internal development and testing purposes only.   
#' 
#' @section Study and treatment
#'
#' Use searchHTTrChem to find chemicals by searching the name/partial name/regex. Enter a zero space  
#' string will retrieve all chemicals. This function returns a dataframe with all the information
#' needed for subsequent queries. Chemicals are identified by: the blind sample identifier
#' (epa_sample_id) or the dsstox_sid (linked to the chemistry dashboard http://comptox.epa.gov).
#' Currently, the raw count data are retrieved by the epa_sample_id and the differential expression 
#' data by the dsstox_sid.  
#'
#' Use getHTTrChemPlates, getHTTrPlateGroups and getHTTrPlateInfo to obtain the complete study
#' design information.
#'
#'
#' @section Raw HTTr data.
#'
#' Use getHTTrProbeCounts to get the raw HTTr data in terms of probe counts using the epa_sample_id. Currently, there are
#' ~21,111 probes and the raw count associated with these probes can be obtained at a sample level.
#' Some of the probe have very low probe counts across samples so use the various options for
#' filtering the data (e.g. average count > 10, median count > 50, etc.)
#'
#'
#' @section Processed HTTr data.
#'
#' Use getHTTrDEG to get the differentially expressed probes/genes for each chemical by dsstox_sid
#' calculated by DESeq2 by comparing each treatment (chemical,concentration) vs DMSO control counts for
#' each probe across multiple replicates, accounting for plate effects. 
#' The default skinny format gives one row per chemical, concentration, probe (gene included), and the #'  log2fc (L2FC) and associated statistics. It is possible to filter these data by p-value,
#' min average probe count (across reps), and |L2FC|. It is also possible to pivot the data into
#' a wide format where rows are chemical treatments and columns are probes/genes. 
#'
#'
