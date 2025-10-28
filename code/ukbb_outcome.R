#computes the patient survival from time of experimental encounter.
get_patients_survival <- function(patients, demog, death_censor_date=DEATH_CENSOR_DATE) {
    patients <- patients %>% select(id, sex, age)
    patients_age <- patients %>%
        left_join(demog %>% select(id, matches("age_.$"), dod), by="id") %>%
        pivot_longer(cols = starts_with("age_"), names_to = "exp_id", values_to = "exp_age") %>% 
        mutate(exp_id = gsub("age_", "", exp_id)) %>% 
        filter(!is.na(exp_age))
    patients_date <- patients %>%
        left_join(demog %>% select(id, matches("date_.$"), dod), by="id") %>%
        pivot_longer(cols = starts_with("date_"), names_to = "exp_id", values_to = "exp_date") %>% 
        mutate(exp_id = gsub("date_", "", exp_id)) %>% 
        filter(!is.na(exp_date))
    patients_survival <- patients_age %>%
        data.table::as.data.table() %>%
        left_join(patients_date,  by = c("id", "sex", "age", "dod", "exp_id")) %>%
        mutate(diff = age - exp_age) %>%
        filter(diff >= 0, diff <= 5) %>% #age model was applied to this experimental encounter (0-3)
        group_by(id, sex, age, dod) %>%
        summarise(exp_age = mean(exp_age), exp_date = mean(exp_date), .groups = "drop") %>%
        as_tibble() %>%
        mutate_at(vars(dod, exp_date), lubridate::as_datetime) %>%
        mutate(dead = !is.na(dod), follow_time = if_else(
            is.na(dod), 
            difftime(death_censor_date, exp_date, units = "days"), 
            difftime(dod, exp_date, units = "days")
        ))
    return(patients_survival)
}

get_patients_disease_outcomes <- function(patients, 
                                          diseases, 
                                          selected_diseases=c('diabetes', 'ckd', 'copd', 'cvd', 'liver'),
                                          censor_date=DISEASE_CENSOR_DATE)
{
    sdiseases <- diseases %>% filter(cohort %in% selected_diseases)
    pd <- plyr::alply(selected_diseases, 1, function(selected_disease) {
        patients_disease <- patients %>% 
            left_join(sdiseases %>% filter(cohort == selected_disease) %>% arrange(age) %>% distinct(id, .keep_all=TRUE) %>% 
                rename(disease_age=age, disease_date=date, disease=cohort), by="id") %>% 
            mutate_at(vars(exp_date, disease_date), lubridate::as_datetime) %>%
            mutate(sick_at_exp=(!is.na(disease_date) & disease_date < exp_date) | (disease_age < exp_age), 
                sick=!is.na(disease_date),
                disease_follow_time=if_else(is.na(disease_date), as.numeric(difftime(censor_date, exp_date, units='days')), 
                    if_else(sick_at_exp, 0, as.numeric(difftime(disease_date, exp_date, units='days')))))
        return(patients_disease %>% mutate(disease=selected_disease))
    }, .parallel=TRUE)
    names(pd) <- selected_diseases
    return(pd)
}