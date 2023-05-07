#' @param demog demographic data returned from demog_data function
#' @param labs contains normalized lab values on extended population
#' @param diseases contains onset of disease for all patients
compute_models_score <- function(demog, labs, diseases, disease_models=c('diabetes', 'ckd', 'copd', 'ncvd', 'liver')) 
{
    longevity_score <- compute_longevity_score(demog, labs, diseases)
    disease_scores <- compute_diseases_score(demog, labs, diseases, disease_models)

    pop_tidy <- longevity_score %>% mutate(disease="longevity") %>% 
        bind_rows(disease_scores %>% filter(!sick) %>% select(-disease_age, -sick)) %>% 
        mutate(disease=factor(disease, levels=c(disease_models, 'longevity'))) %>% 
        group_by(id, sex, age, disease) %>% 
        summarize_at(vars(score), mean, .groups="drop") %>% 
        ungroup
    pop <- pop_tidy %>% pivot_wider(id_cols=c(id, sex, age), names_from=disease, values_from=score)
    pop <- pop %>% group_by(age, sex) %>% mutate(longevity_q=ecdf(longevity)(longevity)) %>% ungroup
    return(pop)
}

compute_longevity_score <- function(demog, labs, diseases) {
    message("computing longevity score")
    multi_models <- readr::read_rds(here::here('data/models/longevity.rds'))
    return(compute_model_score(multi_models, demog, labs, diseases))
}

compute_diseases_score <- function(demog, labs, diseases, disease_models) {
    purrr::map_df(disease_models, function(disease) {
        message("computing ", disease)
        multi_models <- readr::read_rds(here::here(paste0('data/models/', disease, '.rds')))
        ds <- compute_model_score(multi_models, demog, labs, diseases)
        return(ds %>% mutate(disease=disease))
    }) %>%  left_join(diseases %>% select(id, disease=cohort, disease_age=age), by=c("id", "disease")) %>%
        mutate(sick=ifelse(!is.na(disease_age) & age-disease_age > 0, 1, 0),
            disease=factor(disease, levels=unique(disease)))
}

compute_model_score <- function(multi_models, demog, labs, diseases) {
    res <- purrr::map2_df(multi_models, names(multi_models), function(model, age_model) {
        cat(age_model, "...")
        age_model <- as.numeric(age_model)
        labs_data <- labs %>% filter(age<age_model)
        if (nrow(labs_data) == 0) {
            return(data.frame())
        }
        disease_data <- diseases %>% filter(age <= age_model)
        features <- build_features(model, age_model, demog, labs_data, disease_data)
        #requiring RBC
        features <- features %>% filter(!is.na(lab.101.quantiles_1_years_minus1095))
        if (nrow(features) == 0) {
            return(data.frame())
        }
        score <- purrr::map_df(1:length(model$model), function(cv) {
            predictor_features <- features[,model$model[[cv]]$feature_names]
            score = predict(model$model[[cv]], data.matrix(predictor_features))
            return(features %>% mutate(sex = c('male', 'female')[2-male]) %>% 
                select(id, sex) %>% 
                mutate(age=age_model, score=score, cv=cv))
        })
        return(score)
    })
    return(res)
}

build_features <- function(model, age, demog, labs_data, disease_data)
{
    feature_set <- model$feature.set
    lab_features <- labs_data %>% mutate(dt = (age-!!age)*365) %>% 
        mutate(dt_breaks=cut(dt, feature_set$lab_time_breaks, 
            labels=gsub("-", "minus", head(feature_set$lab_time_breaks, -1)), right=FALSE)) %>% 
        filter(!is.na(dt_breaks)) %>% 
        mutate(feature_name=paste0(track, '.quantiles_1_years_', dt_breaks)) %>% 
        rename(feature_value=q) %>% 
        filter(feature_name %in% feature_set$feature_labels$feature) %>% 
        select(id, feature_name, feature_value) %>% 
        group_by(id, feature_name) %>% 
        summarize(feature_value = mean(feature_value), .groups="drop")
    dx_features <- disease_data %>% mutate(track = paste0('WZMN.', cohort)) %>% 
        inner_join(feature_set$dx_tracks, by="track") %>% 
        mutate(feature_name = paste0(track, "_", gsub("-", "minus", dt1), "_", gsub("-", "minus", dt2))) %>% 
        distinct(id, feature_name) %>% 
        mutate(feature_value = 1) %>% 
        filter(id %in% lab_features$id)

    ids <- unique(c(lab_features$id, dx_features$id))

    #adding female/male/age info
    features_tidy <- data.frame(id=ids, feature_name="age", feature_value=age) %>% 
        bind_rows(
            demog %>% 
            filter(id %in% ids) %>% 
            select(id, sex) %>% 
            mutate(feature_name="female", feature_value=sex == 2) %>% 
            select(id, feature_name, feature_value)
        ) %>% 
            bind_rows(
            demog %>% 
            filter(id %in% ids) %>% 
            select(id, sex) %>% 
            mutate(feature_name="male", feature_value=sex == 1) %>% 
            select(id, feature_name, feature_value)
        ) %>% 
        bind_rows(lab_features) %>% 
        bind_rows(dx_features)

    #moving from tidy format
    features <- features_tidy %>% pivot_wider(id_cols=id, names_from=feature_name, values_from=feature_value)
    
    #adding missing columns
    all_features <- unique(unlist(purrr::map(model$model, ~ .x$feature_names)))
    missing_features <- setdiff(all_features, colnames(features))
    missing_data <- data.frame(matrix(NA, ncol=length(missing_features), nrow=nrow(features)))
    colnames(missing_data) <- missing_features
    full_features <- cbind(features, missing_data)

    #setting missing diesease values to 0
    dx_feature_names <- feature_set$feature_labels %>% filter(type == "diagnosis") %>% pull(feature)
    full_features[,dx_feature_names][is.na(full_features[,dx_feature_names])] <- 0

    return(full_features) 
}


#note that model predicts q longevity score from disease scores 
compute_pop_residuals <- function(pop)
{
	#compute residuals of longevity score using diseases
	residuals_model <- readr::read_rds(paste0(models_dir, 'residual_xgboost_model.rds'))
	pre_res <- pop %>% pivot_longer(c("diabetes", "ckd", "copd", "ncvd", "liver"), names_to = "disease", values_to = "score") %>%                        
            replace_na(replace = list(score = 1)) %>% 
            mutate(disease=factor(disease, levels=unique(disease))) %>% 
            pivot_wider(names_from=disease, values_from=score) %>%
            mutate(g=ifelse(gender == 'male', 1, 2))
    pre_res$prediction <- predict(residuals_model, xgboost::xgb.DMatrix(data=pre_res %>% select(age, gender=g, diabetes, ckd, copd, cardio_hr=ncvd, liver) %>% as.matrix()))
    return(pre_res %>% mutate(residual_q=q-prediction))
}
########
# Markovian Lifelong Disease Predisposition models
compute_lifelong_disease_risk <- function(disease_models=c('diabetes', 'ckd', 'copd', 'ncvd', 'liver')) {
    return(purrr::map_df(disease_models, function(disease) {
        mm <- data.table::fread(here::here(paste0('data/models/', disease, '.markov.lifelong.csv')))
        return(mm %>% mutate(disease=disease))
    }) %>% filter(sbin != "no_cbc") %>% mutate(qmax=1/20*as.numeric(sbin), qmin=qmax-1/20))
}


########
# Markovian longevity models
compute_lifelong_longevity_risk <- function(target_pop) {
    mm <- tgutil::fread(here::here('data/models/longevity.markov.lifelong.csv')) %>% #reading in markov models
        filter(sbin != "death", sbin != "no_cbc") %>% # removing death and no-cbc source bins
        mutate(qmax=1/20*as.numeric(sbin), qmin=qmax-1/20) %>% 
        mutate(q=pmean(qmax, qmin)) %>% 
        rename(sex=gender) %>% 
        group_by(sex, age) %>% 
        summarize(model=list(approxfun(c(0, q, 1), c(min(death), death, max(death)))), .groups="drop")
    target_risk <- plyr::ddply(target_pop, plyr::.(age, sex), function(x) {
        f <- mm %>% filter(age == x$age[1], sex==x$sex[1]) %>% pull(model) %>% .[[1]]
        return(x %>% mutate(risk=f(x %>% pull(longevity_q))))
    }, .parallel=T)
}

########
# umap projection
umap_projection <- function(pop) {
	all_umap_models <- readr::read_rds(paste0(data_dir, 'longevity_models/umap_models.rds'))
	names(all_umap_models) <- as.numeric(names(all_umap_models))
	umap_models <- all_umap_models[as.character(sort(pop %>% distinct(age) %>% pull(age)))]
	projection <- runonce::save_run(purrr::map2_df(umap_models, as.numeric(names(umap_models)), function(u, a) {
		message(a)
		d <- pop %>% filter(age == a) %>% na.omit()
		p <- predict(u, d %>% select(q, diabetes, ckd, copd, cardio_hr=ncvd, liver))
		d$x=p[,1]
		d$y=p[,2]
		return(d)
	}), file=paste0(processed_dir, 'ukbb_clalit_umap_projection.rds'))
	u <- projection %>% reshape2::melt(id.vars=c("id", "gender", "age", "x", "y", "longevity")) %>%
		mutate(variable=factor(variable,levels=c('q', 'diabetes', 'ckd', 'copd', 'ncvd', 'liver')))

	#fixing orientation:
	u$x[u$age %in% c(45, 50, 55, 70)] <- -u$x[u$age %in% c(45, 50, 55, 70)]
	u$y[u$age %in% c(55, 60, 70)] <- -u$y[u$age %in% c(55, 60, 70)]
	g <- ggplot(u %>% filter(age<75), aes(x=x, y=y, colour=value)) + geom_point(size=0.005, alpha=0.3) + 
		scale_color_gradientn(colors=ocean.balance(20)) + facet_grid(variable~age) + theme_bw() +
		theme(strip.background=element_blank())
	png(paste0(fig_dir, 'umap.png'), width=1200, height=1080)
	print(g)
	dev.off()
	return(umap_projection)
}




