
#' load entire annotation dataset of ukbiobank
load_data <- function() {
    wd <- getwd()
    setwd(here::here('data'))
    data <- ukbtools::ukb_df('full')
    exclusion <- data.table::fread('exclusion.csv', data.table=FALSE) %>% rename
    setwd(wd)
    ukbb_data <- data %>% filter(!(eid %in% exclusion[,1]))
}

#' extract demographic data from dataset: sex, dob, race, date of interactions, dod
#' @param data full data returned from load_data function 
get_demog_data <- function(data) {
    demog <- data %>% 
        select(
            id=eid, 
            sex=sex_f31_0_0, 
            month_of_birth=month_of_birth_f52_0_0, 
            year_of_birth=year_of_birth_f34_0_0,
            date_of_death_0=date_of_death_f40000_0_0,
            date_of_death_1=date_of_death_f40000_1_0,
            race_0 = ethnic_background_f21000_0_0,
            race_1 = ethnic_background_f21000_1_0,
            race_2 = ethnic_background_f21000_2_0,
            date_0=date_of_attending_assessment_centre_f53_0_0,
            date_1=date_of_attending_assessment_centre_f53_1_0,
            date_2=date_of_attending_assessment_centre_f53_2_0,
            date_3=date_of_attending_assessment_centre_f53_3_0
        ) %>% 
        filter(!is.na(sex), !is.na(month_of_birth), !is.na(year_of_birth)) %>% 
        mutate(
            dob=as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"),
            sex=2-as.numeric(sex),
            dod_showcase=pmin(date_of_death_0, date_of_death_1, na.rm=TRUE),
            age_0 = as.numeric(date_0 - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365,
            age_1 = as.numeric(date_1 - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365,
            age_2 = as.numeric(date_2 - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365,
            age_3 = as.numeric(date_3 - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365
        )

    #adding death information from ukbb data portal (downloaded jan 2022, censor date 30/9/2021)
    dod_portal <- data.table::fread(here::here('data/portal_death_jan_2022.txt')) %>%
        mutate(dod_portal=dmy(date_of_death)) %>% 
        select(id=eid, dod_portal) %>% 
        distinct
    demog <- demog %>% left_join(dod_portal, by="id") %>% 
        mutate(dod=if_else(is.na(dod_showcase), dod_portal, dod_showcase)) %>% 
        select(-dod_showcase, -dod_portal)
    return(demog)
}

#' extract icd diagnosis data (both from hospitals and self reported)
#' @param data full data returned from load_data function
#' @param demog demographic data returned from get_demog_data function
get_diagnosis_data <- function(data, demog) {
    ###########################
    #   HOSPITAL
    ###########################
    #showcase fields
    icd10_value <- data %>% select(all_of(c("eid", grep('f41270', colnames(data), value=TRUE)))) %>% 
        pivot_longer(cols=!eid, names_to="variable", values_to="value")
    icd10_date <- data %>% select(all_of(c("eid", grep('f41280', colnames(data), value=TRUE)))) %>% 
        pivot_longer(cols=!eid, names_to="variable", values_to="value")
    icd10_showcase <- cbind(icd10_value %>% select(id=eid, icd=value), icd10_date %>% select(date=value)) %>% filter(!is.na(icd)) 
    icd9_value <- data %>% select(all_of(c("eid", grep('f41271', colnames(data), value=TRUE)))) %>% 
        pivot_longer(cols=!eid, names_to="variable", values_to="value")
    icd9_date <- data %>% select(all_of(c("eid", grep('f41281', colnames(data), value=TRUE)))) %>% 
        pivot_longer(cols=!eid, names_to="variable", values_to="value")
    icd9_showcase <- cbind(icd9_value %>% select(id=eid, icd=value), icd9_date %>% select(date=value)) %>% filter(!is.na(icd)) 

    #hesin followup period
    hesin <- data.table::fread(here::here('data/portal_hesin_jan_2022.txt'))
    hesin_dx <- data.table::fread(here::here('data/portal_hesin_diag_jan_2022.txt')) %>% 
        left_join(hesin, by = c("eid", "ins_index")) %>% 
        select(eid, diag_icd10, diag_icd9, epistart, admidate) %>% 
        mutate(date=ifelse(is.na(epistart) | epistart == "", admidate, epistart)) %>% 
        mutate(date=dmy(date))
    #icd 10 codings
    icd10_hesin <- hesin_dx %>% 
        select(id=eid, icd=diag_icd10, date) %>% 
        filter(!is.na(icd), icd != "")  %>% 
        arrange(id, date) %>% 
        distinct(id, icd, .keep_all=T) #looking at onset
    #icd9 codings
    icd9_hesin <- hesin_dx %>% 
        select(id=eid, icd=diag_icd9, date) %>% 
        filter(!is.na(icd), icd != "") %>% 
        arrange(id, date) %>% 
        distinct(id, icd, .keep_all=T) #looking at onset

    icd10_codes <- data.table::fread(here::here('data/icd10_coding19.tsv'))
    icd9_codes <- data.table::fread(here::here('data/icd9_coding87.tsv'))
    hospital_icd <- icd10_showcase %>% 
        bind_rows(icd10_hesin) %>% 
        arrange(id, date, icd) %>% 
        distinct(id, icd, .keep_all=T) %>% 
        filter(!is.na(icd))  %>% 
        left_join(icd10_codes %>% select(icd=coding, meaning), by="icd") %>% 
        mutate(icd_version=10) %>% 
        bind_rows(
            icd9_showcase %>% 
            bind_rows(icd9_hesin) %>% 
            arrange(id, date, icd) %>% 
            distinct(id, icd, .keep_all=T) %>% 
            filter(!is.na(icd))  %>% 
            left_join(icd9_codes %>% select(icd=coding, meaning), by="icd") %>% 
            mutate(icd_version=9)
        )

    ###########################
    # first occurrences
    ###########################
    major <- icd10_codes %>% 
        tidyr::separate(coding, into=c('major', 'minor'), 3) %>% 
        distinct(major) %>% 
        left_join(icd10_codes %>% select(major=coding, meaning), by="major")

    fields <- purrr::map_df(major$major, function(x) { 
        field=grep(paste0("^date_",tolower(x)), colnames(data), value=TRUE)
        return(data.frame(major=x, field=ifelse(length(field)==0, NA, field))) 
    }) %>% filter(!is.na(field))
    columns <- c("eid", fields$field)
    occurence_data <- data %>% select(all_of(columns)) %>% pivot_longer(cols=!eid, names_to="icd", values_to="date")
    first_occurrence_icd <- plyr::adply(major, 1, function(x) { 
        field <- grep(paste0("^date_",tolower(x$major)), colnames(data), value=TRUE)
        if (length(field) == 0) {
            return()
        }
        return(data %>% select(all_of(c("eid", field))) %>% 
               rename(id=eid) %>%
               rename(date = !!field) %>%
            #    rename_at(all_of(c("eid", field)), ~ c("id", "date")) %>% 
                # rename_with(~ c("id", "date"), where(all_of(c("eid", field)))) %>%
               filter(!is.na(date)) %>% 
               mutate(icd=x$major, meaning=x$meaning, icd_version=10))
    }, .parallel=FALSE)

    ##########################
    # self rep (nurse)
    ##########################
    self_rep_to_icd_codes <- data.table::fread(here::here('data/self_report_medical_coding_f20002_coding609.tsv')) %>% 
        rename(code=coding, icd=meaning) %>% mutate(code=as.character(code))
    #source for self rep is field 20002
    self_rep_data <- data %>% select(all_of(c("eid", grep('f20002', colnames(data), value=TRUE)))) %>% reshape2::melt(id.var='eid')
    self_rep_age <- data %>% select(all_of(c("eid", grep('f20009', colnames(data), value=TRUE)))) %>% reshape2::melt(id.var='eid')
    self_rep_showcase <- cbind(self_rep_data %>% select(id=eid, code=value), self_rep_age %>% select(age=value)) %>% 
        filter(!is.na(code), age > 0) %>% left_join(self_rep_to_icd_codes, by="code") %>% 
        filter(!is.na(icd)) %>% mutate(icd_version=10) %>% 
        select(-code) %>% 
        left_join(icd10_codes %>% select(icd=coding, meaning), by="icd")

    ##########################
    # self rep (touchscreen)
    ##########################
    self_rep_touchscreen_to_icd_codes <- data.table::fread(here::here('data/self_rep_touchscreen_to_icd10_category100044.csv'))
    #source for self rep is field category 100044
    self_rep_touch <- plyr::adply(self_rep_touchscreen_to_icd_codes, 1, function(rep) {
        return(
            data %>% select(all_of(c("eid", grep(rep$field, colnames(data), value=TRUE)))) %>% 
            reshape2::melt(id.var='eid') %>% 
            filter(!is.na(value), value > 1) %>% 
            select(id=eid, age=value) %>% 
            mutate(icd=rep$icd, icd_version=10, age=as.numeric(age)) %>% 
            left_join(icd10_codes %>% select(icd=coding, meaning), by="icd")
        )
    }, .parallel=FALSE)
    ##########################
    # clinic
    ##########################
    gp_clinic_read2_icd10 <- data.table::fread(here::here('data/gp_clinic_read2_to_icd10_coding1834.tsv'))
    gp_clinic_read3_icd10 <- data.table::fread(here::here('data/gp_clinic_read3_to_icd10_coding1835.tsv'))

    clinic <- data.table::fread(here::here('data/gp_clinic.txt'))
    clinic <- clinic %>% 
        filter(!is.na(read_2), read_2 != "") %>% 
        select(id=eid, date=event_dt, coding=read_2) %>% 
        inner_join(gp_clinic_read2_icd10, by="coding") %>% 
        bind_rows(
            clinic %>% filter(!is.na(read_3), read_3 != "") %>% 
            select(id=eid, date=event_dt, coding=read_3) %>% 
            inner_join(gp_clinic_read3_icd10, by="coding")
        ) %>% 
        rename(icd=meaning) %>% 
        mutate(icd_version=10, date=as.Date(date, format="%d/%m/%Y")) %>% 
        filter(!is.na(date)) %>% 
        left_join(icd10_codes %>% select(icd=coding, meaning), by="icd") %>% 
        select(-coding)

    #combine it all together
    diagnosis  <- hospital_icd %>% mutate(rep="hospital") %>% 
        bind_rows(first_occurrence_icd %>% select(-major) %>% mutate(rep="first_occurrence")) %>% 
        bind_rows(clinic %>% mutate(rep="clinic")) %>% 
        left_join(demog %>% select(id, month_of_birth, year_of_birth), by="id") %>% 
        mutate(age = as.numeric(date - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365) %>% 
        select(-month_of_birth, -year_of_birth) %>% 
        bind_rows(
            self_rep_showcase %>% mutate(rep="self_rep_nurse") %>% 
            bind_rows(self_rep_touch %>% mutate(rep="self_rep_touchscreen")) %>% 
            left_join(demog %>% select(id, month_of_birth, year_of_birth), by="id") %>% 
            mutate(date = as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y") + floor(age*365)) %>% 
            select(-month_of_birth, -year_of_birth)
        ) %>% 
        arrange(id, age) %>% 
        filter(age>0)
    return(diagnosis)
}

#' associate the age of the participant at each visit to the center
#' @param demog demographic data returned from get_demog_data function
get_visit_data <- function(demog) {
    demog_tidy <- demog %>% 
        select(id, sex, age_0, age_1, age_2, age_3) %>% 
        reshape2::melt(id.var= c("id", "sex"))
    return(demog_tidy %>% 
        left_join(
            demog_tidy %>% 
                count(variable) %>% 
                tidyr::separate(variable, into=c("age", "timepoint"), sep="_", remove=FALSE) %>% 
                select(-n), 
            by="variable"
        ) %>% 
        select(id, sex, value, timepoint) %>% 
        rename(age=value) %>% 
        mutate(age = as.numeric(age), timepoint=as.numeric(timepoint))
    )
}


#' extract all lab test measurements
#' @param data full data returned from load_data function 
#' @param visits the age of the patients at each visit to the center, as returned from get_visit_data

get_labs_data <- function(data, visits) {
    lab_codes <- data.table::fread(here::here('data/blood_assays_count_and_biochemistry_category_100080.csv'))
    labs_tidy <- plyr::adply(lab_codes, 1, function(lc) {
        columns <- grep(paste0('_f', lc[1,1], '_'), colnames(data), value=TRUE)
        return(data %>% select(all_of(c('eid', columns))) %>% 
            pivot_longer(cols=columns, names_to='feature_name', values_to='value') %>% 
            rename(id=eid) %>%
            filter(!is.na(value)))
    }, .parallel=TRUE) %>% mutate(feature_name=as.character(feature_name))
    feature_timepoint <- labs_tidy %>% 
        distinct(feature_name) %>% 
        mutate(timepoint = as.numeric(unlist(purrr::map(strsplit(feature_name, "_"), ~ .x[length(.x)-1]))))
    labs_tidy <- labs_tidy %>% 
        left_join(feature_timepoint, by="feature_name") %>% 
        left_join(visits, by=c("id", "timepoint")) %>% 
        # filter(!is.na(as.numeric(value))) %>% #removing values that were textual (e.g. urin result flag)
        filter(!is.na(readr::parse_number(value))) %>% #removing values that were textual (e.g. urin result flag)
        mutate(value=as.numeric(value)) 
    return(labs_tidy)
}

#' identify the onset of any of the chronic diseases from either diagnosis icd or cancer report
#' @param data full data returned from load_data function 
#' @param diagnosis contains all reported diagnosis from various sources as returned from get_diagnosis
#' @param demog demographic data returned from get_demog_data function
get_chronic_onset <- function(data, diagnosis, demog) {
    chronic_icd_codes <- data.table::fread(here::here('data/icd9_icd10_non_healthy_codes_longevity_S2.csv'), header=TRUE)
    chronic_icd10 <- diagnosis %>% filter(icd_version == 10) %>% 
        filter(icd %in% chronic_icd_codes$icd10 | substring(icd, 1, 3) %in% chronic_icd_codes$icd10) %>% 
        left_join(demog %>% select(id, year_of_birth, month_of_birth), by="id") %>% 
        mutate(age = as.numeric(date - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365)

    chronic_icd9 <- diagnosis %>% filter(icd_version == 9) %>% 
        filter(icd %in% chronic_icd_codes$icd9 | substring(icd, 1, 3) %in% chronic_icd_codes$icd9) %>% 
        left_join(demog %>% select(id, year_of_birth, month_of_birth), by="id") %>% 
        mutate(age = as.numeric(date - as.Date(paste0("1-", month_of_birth, "-", year_of_birth), format="%d-%m-%Y"))/365)

    cancer <- data %>% select("eid", grep('interpolated_age_of_participant_when_cancer_first_diagnosed_f20007', colnames(data), value=TRUE))  %>% 
        rename(id=eid) %>% 
        reshape2::melt(id.var='id') %>% 
        filter(!is.na(value)) %>% 
        group_by(id) %>% summarize(age=min(value))

    chronic <- chronic_icd10 %>% select(id, age) %>% 
        bind_rows(chronic_icd9) %>% 
        bind_rows(cancer) %>% 
        group_by(id) %>% 
        summarize(age = min(age))
    return(chronic)
}

#' identify the diagnosis codes of each of the target diseases (diabetes, ckd, copd, ncvd, liver disease)
#' @param diagnosis contains all reported diagnosis from various sources as returned from get_diagnosis
get_diseases <- function(diagnosis, cancer_codes) {
    disease_codes <- tgutil::fread(here::here('data/common_disease_inclusion_icd9_icd10.csv')) %>% 
        bind_rows(cancer_codes %>% mutate(subs=0, cohort='disease.CANCER'))
    #adding exact matches of disease codes
    disease_exact <- diagnosis %>% inner_join(disease_codes %>% filter(subs==0), by=c("icd", "icd_version"), multiple="all")
    #searching the sub-codes
    disease_subs <- plyr::adply(disease_codes %>% filter(subs==1), 1, function(dc) {
        return(diagnosis %>% 
            filter(icd_version == dc$icd_version) %>% 
            filter(grepl(paste0("^", dc$icd), icd)) %>% 
            mutate(cohort=dc$cohort)
            )
        }, .parallel=F)
    return(disease_exact %>% bind_rows(disease_subs) %>% select(id, age, date, cohort) %>% 
        arrange(id, age) %>% distinct(id, cohort, .keep_all=TRUE))
}

#' normalize the raw data measurements (quantile normalization by age/sex)
#' @param labs all lab test measurements as returned from get_labs_data function
#' @param visits the age of the patients at each visit to the center, as returned from get_visit_data
#' @param chronic contains onset of chronic disease as returned from get_chronic_diseases function
#' @param age_breaks the age resolution used for computing quantile normalization
normalize_labs <- function(labs, visits, chronic, age_breaks=seq(35, 85, by=5)) {
    labs_ukbb_to_clalit <- data.table::fread(here::here('data/ukbb_lab_field_to_clalit_lab.csv'), header=TRUE)
    l <- labs %>% left_join(labs_ukbb_to_clalit, by="field")

    #need to fix units
    labs_fix_units <- l %>% filter(action != "")
    labs_fix_units$value <- sapply(paste(labs_fix_units$value, labs_fix_units$action), function(x) eval(parse(text=x)))
    l <- l %>% filter(action == "") %>% bind_rows(labs_fix_units)

    message("finished fixing units")
    age_centers <- zoo::rollmean(age_breaks, 2)
    labs_chronic <- l %>% left_join(visits, by = c("id", "timepoint", "sex", "age")) %>% 
        mutate(age_group = cut(age, age_breaks, right=FALSE)) %>%
        left_join(chronic %>% rename(dx_age=age), by="id") %>% 
        mutate(healthy= ifelse(is.na(dx_age) | dx_age > age, 1, 0))

    message("computing lab quantiles")
    normalizer <- .compute_lab_quantiles(labs_chronic %>% filter(healthy==1, !is.na(sex)), quantile_probs=NULL)
    #going over all lab measurements and normalizing them
    norm_labs <- plyr::ddply(labs_chronic %>% filter(!is.na(sex)), 
        plyr::.(track, lab, age_group, sex), 
        function(x) {
            age_sex_lab_normalizer <- normalizer[[paste(x$track[1], x$lab[1], x$age_group[1], x$sex[1], sep=".")]]
            if (is.null(age_sex_lab_normalizer)) {
                return(x %>% mutate(q=NA))	
            }
            return(x %>% mutate(q=age_sex_lab_normalizer(value)))
        }, .parallel=FALSE)	
    norm_labs$lab <- factor(norm_labs$lab, levels=labs_ukbb_to_clalit$lab)
    return(norm_labs)
}

#' computed ecdf for each lab test, age group and sex after unit conversion of lab tests
#' @param labs_chronic all lab test measurements after unit modification and annotation of chronic state
#' @param quantile_probs preset quantiles, to be used when stats are insufficient or predefined resolution is required
.compute_lab_quantiles <- function(healthy_labs, quantile_probs) {
    normalizer <- plyr::dlply(healthy_labs, plyr::.(track, lab, age_group, sex), 
        function(x) {
            if (nrow(x) < 10) {
                return(NULL)
            }
            if (is.null(quantile_probs)) {
                return(ecdf(x$value))
            } else {
                return(quantile(x$value, quantile_probs))
            }
        }, .parallel=TRUE)
    return(normalizer)
}

    
    
get_parents_survival <- function(data) {
    #mothers 
    mothers_death_fields <- grep("f3526", colnames(ukbb_data), value=T)
    mothers_age_fields <- grep("f1845", colnames(ukbb_data), value=T)
    mothers_alive_fields <- grep("f1835", colnames(ukbb_data), value=T)
    #fathers
    fathers_death_fields <- grep("f1807", colnames(ukbb_data), value=T)
    fathers_age_fields <- grep("f2946", colnames(ukbb_data), value=T)
    fathers_alive_fields <- grep("f1797", colnames(ukbb_data), value=T)

    m_age_at_death <- matrixStats::rowMins(as.matrix(data[,mothers_death_fields]), na.rm=TRUE)
    m_age_censoring <- matrixStats::rowMins(as.matrix(data[,mothers_age_fields]), na.rm=TRUE)
    f_age_at_death <- matrixStats::rowMins(as.matrix(data[,fathers_death_fields]), na.rm=TRUE)
    f_age_censoring <- matrixStats::rowMins(as.matrix(data[,fathers_age_fields]), na.rm=TRUE)
    
    parents <- data.frame(id=data$eid, mother_age_at_death=m_age_at_death, mother_last_alive=m_age_censoring,
        father_age_at_death=f_age_at_death, father_last_alive=f_age_censoring) %>% 
        mutate(mdead=!is.infinite(mother_age_at_death), 
               mfollow_time=ifelse(is.infinite(mother_age_at_death), ifelse(is.infinite(mother_last_alive), NA, mother_last_alive), mother_age_at_death),
               fdead=!is.infinite(father_age_at_death), 
               ffollow_time=ifelse(is.infinite(father_age_at_death), ifelse(is.infinite(father_last_alive), NA, father_last_alive), father_age_at_death))
    return(parents)
}

get_socio_demog <- function(data, visits) {
	socio_fields <-  grep('f738', colnames(ukbb_data), value=T)
	socio_tidy <- data %>% select(all_of(c('eid', socio_fields))) %>% 
		reshape2::melt(id.var="eid") %>% 
		rename(id=eid, feature_name=variable) %>%
		filter(!is.na(value)) %>% 
		mutate(feature_name=as.character(feature_name)) %>% 
 		mutate(timepoint = as.numeric(unlist(purrr::map(strsplit(feature_name, "_"), ~ .x[length(.x)-1]))))

	socio_tidy <- socio_tidy %>% left_join(visits) %>% 	mutate(value=as.numeric(value)) %>% 
		select(id, gender, age, socio=value)
	return(socio_tidy)
}



build_cancer_icd9_icd10_dictionary <- function(data) {
    icd10_cancer_codes <- data %>% select("eid", 
        grep('type_of_cancer_icd10_f40006', colnames(data), value=TRUE)) %>% 
        reshape2::melt(id.var="eid") %>% filter(!is.na(value)) %>% distinct(value)
    icd9_cancer_codes <- data %>% select("eid", 
        grep('type_of_cancer_icd9_f40013', colnames(data), value=TRUE)) %>% 
        reshape2::melt(id.var="eid") %>% filter(!is.na(value)) %>% distinct(value)
    cancer_codes <- icd10_cancer_codes %>% mutate(icd_version = 10) %>% 
        bind_rows(icd9_cancer_codes %>% mutate(icd_version=9)) %>% 
        rename(icd=value)
    return(cancer_codes)
}


#note that cancer entries are not updated with latest HESIN data
get_cancer <- function(data, demog) {
	icd_cancer <- data %>% select("eid", 
		grep('type_of_cancer_icd10_f40006', colnames(data), value=TRUE),
		grep('type_of_cancer_icd9_f40013', colnames(data), value=TRUE),
		grep('age_at_cancer_diagnosis_f40008', colnames(data), value=TRUE))
	colnames(icd_cancer) <- gsub('age_at_cancer_diagnosis_f40008_', 'age_', 
		gsub('type_of_cancer_icd9_f40013_', 'cancer9_', 
		gsub('type_of_cancer_icd10_f40006_', 'cancer10_', colnames(icd_cancer))))
	
	entries <- purrr::map_df(strsplit(grep("age", colnames(icd_cancer)[-1], invert=TRUE, value=TRUE), "_"), ~ data.frame(type=.x[1], tp=.x[2], idx=.x[3]))
	icd_cancer <- plyr::adply(entries, 1, function(x) { 
		d <- icd_cancer[,c("eid", paste0(x$type, "_", x$tp, "_", x$idx), paste0("age_", x$tp, "_", x$idx))]
		colnames(d) <- c("id", "cancer", "age")
		d <- d %>% filter(!is.na(cancer)) %>% select(id, coding=cancer, age)
		if(x$type == "cancer10") {
			d <- d %>% left_join(icd10_codes %>% select(coding, meaning)) %>% mutate(coding_ref=19)
		} else {
			d <- d %>% left_join(icd9_codes %>% select(coding, meaning)) %>% mutate(coding_ref=87)
		}
		return(d)
	}) %>% select(id, coding, age, meaning)


	###########################
	#	self_reported
	###########################
	self_rep_cancer <- data %>% select("eid", 
		grep('cancer_code_selfreported_f20001', colnames(data), value=TRUE),
		grep('interpolated_age_of_participant_when_cancer_first_diagnosed_f20007', colnames(data), value=TRUE))
	colnames(self_rep_cancer) <- gsub('interpolated_age_of_participant_when_cancer_first_diagnosed_f20007_', 'age_', gsub('cancer_code_selfreported_f20001_', 'cancer_', colnames(self_rep_cancer)))

	entries <- purrr::map_df(strsplit(colnames(self_rep_cancer)[-1], "_"), ~ data.frame(type=.x[1], tp=.x[2], idx=.x[3])) %>% distinct(tp, idx)
	self_rep_cancer <- plyr::adply(entries, 1, function(x) { 
		d <-  self_rep_cancer[,c("eid", paste0("cancer_", x$tp, "_", x$idx), paste0("age_", x$tp, "_", x$idx))]
		colnames(d) <- c("id", "cancer", "age")
		return(d %>% filter(!is.na(cancer)))
	}) %>% select(id, coding=cancer, age) %>% mutate(coding=as.numeric(coding)) %>% 
	left_join(cancer_codes %>% select(coding, meaning)) %>% mutate(coding_ref=)

	return(cancer_merged)
}
