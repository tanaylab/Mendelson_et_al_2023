#this file contains methods for selecting patients with different criteria and extracts requested data for these patients.


#create global usefull filters to select patients that were already born, didn't die, were registered in the healthcare system and did not leave if for good.
emr_filter.create('born', 'patients.dob', time.shift=c(-120,0)*year())
emr_filter.create('dead', 'patients.dod', time.shift=c(-120,0)*year())
emr_filter.create('registered', 'patients.status.register', time.shift=c(-120,0)*year())
emr_filter.create('left_for_good', 'patients.status.lfg', time.shift=c(-120,0)*year())
emr_filter.create('absent', 'patients.status.absent', time.shift=c(-1,0)*month())

#create vtrack for sex
emr_vtrack.create('sex', 'patients.dob', time.shift=c(-120,0)*year(), func='earliest') #function doesn't matter, only a single entry per patient

#create vtrack for age
emr_vtrack.create('age', 'patients.dob', time.shift=c(-120,0)*year(), func='dt2.earliest') #function computes the difference in time (in hours) between the iterator time and birth time

#create vtrack for time_to_death
emr_vtrack.create('time_to_death', 'patients.dod', time.shift=c(0,120)*year(), func='dt1.earliest') #function computes the difference in time (in hours) between the iterator time and time of death

#create vtrack for time_to_left_for_good
emr_vtrack.create('time_to_left_for_good', 'patients.status.lfg', time.shift=c(0,120)*year(), func='dt1.earliest') #function computes the difference in time (in hours) between the iterator time and time the patient left for good


#' Select all patients that are alive at a given age in the naryn database 
#' @param target_age - the age of all patients to identify
#' @param start_time - the minimum time in which to start looking for patients
#' @param end_time - the maximum time in which to start looking for patients
#' @param filters - string containing any additional filters to apply to the list of patients
get_all_patients_at_age <- function(target_age, start_time, end_time, end_db, filters=NULL, additional_data=NULL, additional_names=NULL) 
{
    age_filter <- emr_filter.create('f_age', 'patients.dob', time.shift=c(-target_age-1, -target_age+1)*year())
    patients_filter=paste(c(filters, '(born & !dead & registered & !left_for_good & !absent & f_age)'), collapse = " & ")
    patients <- emr_extract(c('age/year()', 'sex', 'time_to_death/year()', 'time_to_left_for_good/year()', additional_data),
                            iterator=list(target_age*year(), 'patients.dob'), 
                            filter=patients_filter,
                            stime=start_time,
                            etime=end_time,
                            names=c('age', 'sex', 'survival', 'left_for_good', additional_names)
                            ) %>% 
    mutate(time_in_system=ifelse(is.na(left_for_good), (end_db - time)/year(), left_for_good))
    ### BUG WORKAROUND ###
    patients <- patients %>% inner_join(emr_extract('age', iterator=patients %>% select(id, time), filter=patients_filter) %>% select(id, time), by=c("id", "time"))
    ### BUG WORKAROUND ###
    return(patients)
}

#' Extract EHR data for a set of patients (id/time)
#' @param patients - dataframe containing patient id and time reference point for lab and diesease extraction
#' @param labs - a vector of tracks representing lab data
#' @param labs_window - a vector of two numbers (e.g. c(-3*year(), 0) representing the time window relative to patient time, in which to look for lab values
#' @param diseases - a vector of tracks representing diseases
#' @param diseases_window - a vector of two numbers (e.g. c(-3*year(), 0) representing the time window relative to patient time, in which to look for diseases

get_patients_features <- function(patients, labs, labs_window, diseases, diseases_window, max_missing_percent=1) {
    #create a vtrack for each lab, representing the time window from which the lab value will be sampled
    purrr::walk(labs, ~ emr_vtrack.create(paste0('v_', .x), .x, time.shift=labs_window, func='avg'))
    
    #create a vtrack for each disease, representing the time window from which the disease value will be sampled
    purrr::walk(diseases, ~ emr_vtrack.create(paste0('v_', .x), .x, time.shift=diseases_window, func='size'))
    
    data <- emr_extract(c(paste0('v_', c(labs, diseases))), iterator=patients, names=c(labs, diseases)) %>% 
        as_tibble() %>% 
        mutate(across(all_of(diseases), ~+as.logical(.x)))
    
    missing_per_feature <- apply(t(data), 1, function(x) { sum(is.na(x)) })
    data <- data[, missing_per_feature/nrow(data) < max_missing_percent] #requiring common features, ignoring sporadic ones
    return(patients %>% 
           select(id, time, age, sex) %>% 
           left_join(data %>% select(-ref), by=c("id", "time")) %>% 
           select(-time) %>% 
           mutate(sex=sex-1)
    )
}