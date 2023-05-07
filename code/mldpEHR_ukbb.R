source(here::here("code/init.R"))
source(here::here("code/ukbb_preprocessing.R"))
source(here::here("code/models.R"))
library(labNorm)
library(mldpEHR)

ukbb_data <- load_data()
ukbb_demog <- get_demog_data(ukbb_data) %cache_df% here('output/ukbb_demog.csv') %>% as_tibble()
ukbb_diagnosis <- get_diagnosis_data(ukbb_data, ukbb_demog)

ukbb_visits <- get_visit_data(ukbb_demog) %cache_df% here('output/ukbb_visits.csv') %>% as_tibble()
ukbb_labs <- get_labs_data(ukbb_data, ukbb_visits) %cache_df% here('output/ukbb_labs.csv') %>% as_tibble() %>% 
    mutate(sex=c('male', 'female')[sex]) %>% 
    inner_join(ln_ukbb_labs() %>% mutate(field=as.numeric(ukbb_code)) %>% select(field))

ukbb_labs$q <- ln_normalize_multi_ukbb(ukbb_labs %>% select(id, lab_code=field, age, sex, value))

ukbb_diseases <- get_diseases(ukbb_diagnosis, build_cancer_icd9_icd10_dictionary(ukbb_data)) %cache_df% here('output/ukbb_diseases.csv') %>% as_tibble()
models_dir <- 'data/models/'
predictors <- c('longevity', 'diabetes', 'ckd', 'copd', 'cvd', 'liver') %>% 
	purrr::set_names() %>% 
	purrr::map(function(m) 
	{
		readr::read_rds(paste0(models_dir, m, '.rds')) %>% 
			purrr::imap( ~ c(.x, age=as.numeric(.y), feature_names=list(unique(unlist(purrr::map(.x$model, ~ .x$feature_names))))))
	})

potential_features <- unique(unlist(purrr::map(predictors, function(predictor) {
	purrr::map(predictor, function(p) {
		p$feature_names
	})
})))


#building features to be used by all predictors (longevity, diseases)
ukbb_to_clalit <- tgutil::fread('data/ukbb_lab_field_to_clalit_lab.csv')
features <- purrr::map2_df(predictors[[1]], names(predictors[[1]]), function(model, age_model) {
	message(age_model)
	age_model <- as.numeric(age_model)
	labs_features <- ukbb_labs %>% filter(age<age_model, age>age_model-5, !is.na(q)) %>% 
        left_join(ukbb_to_clalit %>% select(field, track), by="field") %>% 
		mutate(feature=paste0(track, '.quantiles_1_years_minus1095')) %>% 
		filter(feature %in% potential_features) %>% 
		group_by(id, feature) %>% summarize(value=mean(q), .groups="drop")

	disease_features <- ukbb_diseases %>% filter(age <= age_model) %>% 
		mutate(feature=paste0('WZMN.', cohort, '_minus43800_0')) %>% 
		filter(feature %in% potential_features) %>% 
		distinct(id, feature) %>% 
		mutate(value=1)

	ids <- unique(c(labs_features$id, disease_features$id))

	#adding female/male/age info
	features_tidy <- data.frame(id=ids, feature="age", value=age_model) %>% 
		bind_rows(ukbb_demog %>% filter(id %in% ids) %>% mutate(feature="male", value= sex==1) %>% select(id, feature, value)) %>% 
		bind_rows(labs_features) %>% 
		bind_rows(disease_features)
				
	#moving from tidy format
	features <- features_tidy %>% pivot_wider(id_cols='id', names_from='feature') %>% 
		mutate(sex=2-male)

	#setting missing diesease values to 0
	disease_feature_names <- grep('WZMN.disease', colnames(features), value=TRUE)
	features[,disease_feature_names][is.na(features[,disease_feature_names])] <- 0

	#adding missing features
	missing_features <- setdiff(potential_features, colnames(features))
	features[,missing_features] <- NA
	
	#requiring RBC
 	features <- features %>% filter(!is.na(lab.101.quantiles_1_years_minus1095))
 	return(features)
}) %cache_df% here('output/ukbb_mldp_features.csv') %>% as_tibble()


predictor_scores <- purrr::map2_df(predictors, names(predictors), ~ mldp_predict_multi_age(features, .x) %>% mutate(predictor=.y))
pop <- predictor_scores %>% filter(predictor == "longevity") %>% select(id, age, sex, longevity=score, longevity_q=quantile) %>% 
	mutate(sex=factor(c('male', 'female')[sex], levels=c('male', 'female'))) %>% 
	left_join(predictor_scores %>% filter(predictor != "longevity") %>% 
		select(id, age, predictor, score) %>% 
		left_join(ukbb_diseases %>% select(id, disease_age=age, predictor=cohort)) %>% 
		mutate(score = ifelse(!is.na(disease_age) & disease_age < age, NA, score)) %>% 
		pivot_wider(id_cols=c("id", "age"), names_from="predictor", values_from="score")
	)

#survival <- get_survival()
