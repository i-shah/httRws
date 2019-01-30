                    
#' Finds HTTr chemicals
#' 
#' `searchHTTrChem` returns the HTTr chemicals than match the input string pattern. 
#' @param txt A text/regex string 
#' @return a dataframe containing the chemical name, dsstox_sid and epa_sample_id
#' @examples
#' searchHTTrChem('phtha')                                  
#' @export                        
searchHTTrChem <- function(txt){
    D1 = getDataFromWS('searchChem',
                       list("txt"=txt))

    if(!is.null(D1$hits)){
        list2df(D1$hits)
    }
}   

#' Get the raw HTTr count data by epa_sample_id
#' 
#' @param epa_sample_id A valid epa_sample_id from searchHTTrChem
#' @param probe_filter removes probes that don't meet minimum requirements
#'    mx0_5:  probe maximum greater than 5
#'    av0_10: probe average greater than 10
#'    md0_20: probe median greater than 20
#'                    
#' @return a named list with two dataframes
#'    treatments: all treatment factors (columns) x samples (rows)
#'    counts: all probe counts (columns) x samples (rows)         
#' @examples
#' Get count data for TP0001718N03 and keep probes with average > 100
#'  getHTTrCountData('TP0001718N03',probe_filter='av0_100')         
#' @export                                                
getHTTrCountData <- function(epa_sample_id,probe_filter='av0_10'){
    D1 = getDataFromWS('getChemData',
                       list("chem_id"=epa_sample_id,
                            "probe_filter",probe_filter))
      
    if(!is.null(D1$treatments) & !is.null(D1$treatments)){
        list(treatments=list2df(D1$treatments),
             counts=list2df(D1$counts))
    }
}

#' Get the raw HTTr count data by epa_sample_id
#' 
#' @param epa_sample_id A valid epa_sample_id (e.g. from searchHTTrChem)
#' @param probe_filter removes probes that don't meet minimum requirements
#'    mx0_5:  probe maximum greater than 5
#'    av0_10: probe average greater than 10
#'    md0_20: probe median greater than 20
#'                    
#' @return a named list with two dataframes
#'    treatments: all treatment factors (columns) x samples (rows)
#'    counts: all probe counts (columns) x samples (rows)                   
#' @examples
#' Get count data for TP0001718N03 and keep probes with average > 100
#'  getChemProbeCounts('TP0001718N03',probe_filter='av0_100')          
#' @export                                                                        
getChemProbeCounts <- function(epa_sample_id,
                               probe_filter='av0_10'){
    D1 = getDataFromWS('getChemProbeCounts',
                       list(chem_id=epa_sample_id,
                            probe_filter=probe_filter))
      
    if(!is.null(D1$treatments) & !is.null(D1$treatments)){
        list(treatments=list2df(D1$treatments),
             counts=list2df(D1$counts))
    }
}


#' Get the  HTTr plate_ids for a chemical
#' 
#' @param epa_sample_id A valid epa_sample_id (e.g. from searchHTTrChem)
#' @return A named list 
#'    plates: all plate_ids with chemical
#' @examples
#' Get plates for TP0001718A17
#'  getHTTrChemPlates('TP0001718A17')      
#' @export                                                                             
getHTTrChemPlates <- function(epa_sample_id){

    D1 <- getDataFromWS('getChemPlates',
                       list("chem_id"=epa_sample_id))
    D1
}

#' Get the  HTTr plate information for a plate_id
#' 
#' @param plate_id A valid plate_id (e.g. from getHTTrChemPlates)
#' @return A dataframe with 
#'    
#' @examples
#' Get plates for plate_id=TC00284691
#'  getHTTrPlateInfo('TC00284691')      
#' @export                                       
getHTTrPlateInfo <- function(plate_id){
    D1 <- getDataFromWS('getPlateInfo',
                       list("plate_id"=plate_id))
    if (!is.null(D1$plateinfo)){
        list2df(D1$plateinfo)
    }
}

#' @export                          
getHTTrPlateGroups <- function(){
    D1 <- getDataFromWS('getPlateGroups',list())
    if (!is.null(D1$plate_groups)){
        list2df(D1$plate_groups)
    }
}

#'  Gets the effect of each chemical treatment on probes  
#' 
#'  This function returns the results from DESeq2 analysis of the raw                 
#'  count data by using the replicates for each sample across multiple plates
#'  as well as the DMSO conrols for each treatment group (conc).
#'  Just calling this with the sid wil return a long skinny format df. The skinny
#'  df has all the information for filtering the results. 
#'  Optionally, these results can be filtered by p-values, minimum average count
#'  l2fc magnitude and also pivoted to create wide df (more suitable for analysis) 
#'
#'  Note:
#'  a) If the chemical was run mulitple times then there will be multiple sample ids!
#'  b) If filters are used to exclude probes by any criterion then the missing value is replaced with 0.0 (naught) 
#'  c) Keep the skinny format without pivoting if you want to see all the data 
#'
#' @param sid A valid dsstox_sid (e.g. from searchHTTrChem)
#' @param mincnt0 the min average raw count of the probe across reps (optional)
#' @param p0 the the unadjusted p-value cut-off for each probe (optional)
#' @param l2fc0 the min abs(L2FC) threshold (optional)
#' @param pivot_on whether to pivot and how. Values
#'      NULL: no pivot provide skinny results
#'      gene: pivot on gene symbol 
#'      probe_id: pivot on probe_id 
#' @param pivot_func:  if pivot_on=='gene' then how to aggregate l2fc
#'      max: take the max +ve or min -ve
#'      mean: take the average 
#'                    
#' @return a named list with two dataframes
#'    l2fc:  all L2FC data by probe_id/genen (columns) x treatment groups (rows)
#'    index: the index columns in the L2FC dataframe 
#' @examples
#'   Get all the probes with differential expression information for DTXSID7020182
#'     X <- getHTTrDEG('DTXSID7020182')
#'
#'   ... above chemical and subset of probes with average count > 100 per replicate 
#'     X <- getHTTrDEG('DTXSID7020182',mncnt0=100)
#'
#'   ... above and just subset of probes unadjusted p-value <= 0.05   
#'     X <- getHTTrDEG('DTXSID7020182',mncnt0=100,p0=0.05)
#'
#'   ... above and just subset of probes unadjusted |l2fc0| >= 2   
#'     X <- getHTTrDEG('DTXSID7020182',mncnt0=100,p0=0.05,l2fc0=2)
#'
#'   ... above and pivot on gene symbol and take the max of l2fc values    
#'     X <- getHTTrDEG('DTXSID7020182',mncnt0=100,p0=0.1,l2fc0=0.6,
#'                      pivot_on='gene', pivot_func='max')
#'
#' @export                                                      
getHTTrDEG <- function(sid,mncnt0=NULL,p0=NULL,l2fc0=NULL,
                        pivot_on=NULL,pivot_func=NULL,
                        svc=SVC){
    Q = list(dsstox_sid=sid)
    
    if (!is.null(pivot_on)){
        Q$pivot_on = pivot_on
        
        if (!is.null(pivot_func)){
            Q$pivot_func=pivot_func
        } else{
            Q$pivot_func='max'
        }
    }
        
    if (!is.null(p0)){
        Q$p0 = p0
    }
    
    if (!is.null(l2fc0)){
        Q$l2fc0=l2fc0
    }
    
    if (!is.null(mncnt0)){
        Q$mncnt0=mncnt0
    }
      
    D1 <- getDataFromWS('getChemDEG',Q)

    if (!is.null(D1$l2fc) & !is.null(D1$cols)){
        list(l2fc=list2df(D1$l2fc),index=D1$cols)
    }
}
                    