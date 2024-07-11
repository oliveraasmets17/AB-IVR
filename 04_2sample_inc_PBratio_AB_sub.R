#perform 2-sample instrumental variable regression
#Two-Sample Two-Stage Least Squares
#instrument: AB; exposure: MB (indexes, genuses, species); outcome: diseases
# AB -> MB regression in MB cohort
# AB -> disease regression in EGC cohort (200k; excluding MB cohort)
# AB usage 2004-2014; disease incidence 2014-2022; exclude recent (6 months) AB users
# sensitivity analysis with AB usage in AB subgroups (macrolides J01FA, fluoroquinolones J01MA, and penicilline J01CR)


####DATA + LIBRARIES####
setwd("C:/Users/Administrator/Documents/01_projektid/09_AB_instrument")
load("02_input/andmed.RData")
load("02_input/andmed_200k.RData")
#AB data in large cohort
ab_200k = readRDS("02_input/MB_AB_data_200k.rds")
#AB data in MB cohort
medic = readRDS("02_input/MB_medications_data.rds")

prevotella_bacteroides = readRDS("02_input/metaphlan_PB_ratio.rds")

library(dplyr)
library(Epi)
library(ggplot2)
library(xlsx)
library(RNOmni)


####PREPROCESSING####
####AB subgroups####
# macrolides (J01FA)
# fluoroquionolones (J01MA)
# penicilline (J01CR)


ab_200k = as.data.frame(ab_200k)
ab_200k = ab_200k %>% select("skood", "MedicationAssembled atcDate", "MedicationAssembled atcDateSource",
                             "MedicationAssembled atc code", "MedicationAssembled atc name",
                             "MedicationAssembled icd10 code",
                             "MedicationAssembled source name", "Baseline_date")
ab_200k$med_date = ab_200k$'MedicationAssembled atcDate'
ab_200k$med_source = ab_200k$'MedicationAssembled atcDateSource'
ab_200k$atc_code = ab_200k$'MedicationAssembled atc code'
ab_200k$atc_name = ab_200k$'MedicationAssembled atc name'
ab_200k$icd10_code = ab_200k$'MedicationAssembled icd10 code'
ab_200k$med_source2 = ab_200k$'MedicationAssembled source name'

ab_200k[c("MedicationAssembled atcDate", "MedicationAssembled atcDateSource",
          "MedicationAssembled atc code", "MedicationAssembled atc name",
          "MedicationAssembled icd10 code",
          "MedicationAssembled source name", "MB_sample_date")] = NULL

#exclude self-reported
ab_200k = ab_200k %>% filter(!med_source2 == "Doonori vastustik")
#exclude E-tervis
ab_200k = ab_200k %>% filter(!med_source2 == "E-Tervis")

#AB data
ab_200k$time_gap = (ab_200k$Baseline_date - ab_200k$med_date)/365.25
ab_200k$recent = (ab_200k$time_gap < 0.5 & ab_200k$time_gap >= 0)*1


#detect recent AB-users (for excluding)
exclude = unique(ab_200k$skood[ab_200k$recent == 1])
ab_200k = ab_200k %>% filter(time_gap > 0.5)
#exclude as information the AB-usage after MB sampling (but don't exclude the individuals)
ab_200k = ab_200k %>% filter(time_gap <= 10)

#detect the usage of J01FA
ab_200k$J01FA <- grepl("\\J01FA", ab_200k$atc_code)*1
#data-frame of J01FA users (one row per usage; one individual can have several rows if used J01FA several times)
J01FA_AB = ab_200k %>% filter(J01FA == 1)

#detect the usage of J01MA
ab_200k$J01MA <- grepl("\\J01MA", ab_200k$atc_code)*1
#data-frame of J01MA users
J01MA_AB = ab_200k %>% filter(J01MA == 1)

#detect the usage of J01CR
ab_200k$J01CR <- grepl("\\J01CR", ab_200k$atc_code)*1
#data-frame of J01CR users
J01CR_AB = ab_200k %>% filter(J01CR == 1)

#count of J01FA per person (find by counting how many rows one individual had in J01FA_AB)
J01FA_usage_200k = J01FA_AB %>% count(skood)
J01FA_usage_200k = rename(J01FA_usage_200k, J01FA_n = n)

#count of J01MA per person
J01MA_usage_200k = J01MA_AB %>% count(skood)
J01MA_usage_200k = rename(J01MA_usage_200k, J01MA_n = n)

#count of J01CR_AB per person
J01CR_usage_200k = J01CR_AB %>% count(skood)
J01CR_usage_200k = rename(J01CR_usage_200k, J01CR_n = n)


##merge AB with phen
andmed_200k = merge(J01FA_usage_200k, andmed_200k, by="skood", all.y=T)
andmed_200k = merge(J01MA_usage_200k, andmed_200k, by="skood", all.y=T)
andmed_200k = merge(J01CR_usage_200k, andmed_200k, by="skood", all.y=T)

#substitute NA with 0 (as non-usage)
andmed_200k$J01FA_n[is.na(andmed_200k$J01FA_n)==T] = 0
andmed_200k$J01MA_n[is.na(andmed_200k$J01MA_n)==T] = 0
andmed_200k$J01CR_n[is.na(andmed_200k$J01CR_n)==T] = 0

##AB subgroups in MB cohort

medic = as.data.frame(medic)
medic = medic %>% select("skood", "MedicationAssembled atcDate", "MedicationAssembled atcDateSource",
                         "MedicationAssembled atc code", "MedicationAssembled atc name",
                         "MedicationAssembled icd10 code", "MedicationAssembled icd10 name",
                         "MedicationAssembled source name", "MB_sample_date")
medic$med_date = medic$'MedicationAssembled atcDate'
medic$med_source = medic$'MedicationAssembled atcDateSource'
medic$atc_code = medic$'MedicationAssembled atc code'
medic$atc_name = medic$'MedicationAssembled atc name'
medic$icd10_code = medic$'MedicationAssembled icd10 code'
medic$icd10_name = medic$'MedicationAssembled icd10 name'
medic$med_source2 = medic$'MedicationAssembled source name'
medic$mb_date = medic$'MB_sample_date'

medic[c("MedicationAssembled atcDate", "MedicationAssembled atcDateSource",
        "MedicationAssembled atc code", "MedicationAssembled atc name",
        "MedicationAssembled icd10 code", "MedicationAssembled icd10 name",
        "MedicationAssembled source name", "MB_sample_date")] = NULL

#exclude self-reported
medic = medic %>% filter(!med_source2 == "Doonori vastustik")
#exclude E-tervis
medic = medic %>% filter(!med_source2 == "E-Tervis")

medic$J01 <- grepl("\\J01", medic$atc_code)*1
medic_AB = medic %>% filter(J01 == 1)

medic_AB$time_gap = (medic_AB$mb_date - medic_AB$med_date)/365.25
medic_AB$recent = (medic_AB$time_gap < 0.5 & medic_AB$time_gap >= 0)*1

#detect recent AB-users (for excluding)
exclude = unique(medic_AB$skood[medic_AB$recent == 1])
#exclude as information the AB-usage after MB sampling (but don't exclude the individuals)
medic_AB = medic_AB %>% filter(time_gap > 0.5)
medic_AB = medic_AB %>% filter(time_gap <= 10)

#AB data
medic_AB$J01FA <- grepl("\\J01FA", medic_AB$atc_code)*1
J01FA_AB = medic_AB %>% filter(J01FA == 1)

medic_AB$J01MA <- grepl("\\J01MA", medic_AB$atc_code)*1
J01MA_AB = medic_AB %>% filter(J01MA == 1)

medic_AB$J01CR <- grepl("\\J01CR", medic_AB$atc_code)*1
J01CR_AB = medic_AB %>% filter(J01CR == 1)

#count of J01FA per person
J01FA_usage_MB = J01FA_AB %>% count(skood)
J01FA_usage_MB = rename(J01FA_usage_MB, J01FA_n = n)

#count of J01MA per person
J01MA_usage_MB = J01MA_AB %>% count(skood)
J01MA_usage_MB = rename(J01MA_usage_MB, J01MA_n = n)

#count of J01CR_AB per person
J01CR_usage_MB = J01CR_AB %>% count(skood)
J01CR_usage_MB = rename(J01CR_usage_MB, J01CR_n = n)

##merge AB with phen
andmed = merge(J01FA_usage_MB, andmed, by="skood", all.y=T)
andmed = merge(J01MA_usage_MB, andmed, by="skood", all.y=T)
andmed = merge(J01CR_usage_MB, andmed, by="skood", all.y=T)

#substitute NA with 0 (as non-usage)
andmed$J01FA_n[is.na(andmed$J01FA_n)==T] = 0
andmed$J01MA_n[is.na(andmed$J01MA_n)==T] = 0
andmed$J01CR_n[is.na(andmed$J01CR_n)==T] = 0


####other preprocessing####
exclude = andmed$skood
EGC = andmed_200k
EGC = EGC %>% filter(!skood %in% exclude)
MB = merge(prevotella_bacteroides, andmed, by="skood")
rm(andmed_200k, andmed, exclude, prevotella_bacteroides, J01CR_AB, J01CR_usage_200k,
   J01FA_AB, J01FA_usage_200k, J01MA_AB, J01MA_usage_200k, ab_200k,
   J01CR_usage_MB, J01MA_usage_MB, J01FA_usage_MB, medic, medic_AB)
# two datasets: MB and EGC

#gender-variable
MB$gender = factor(MB$gender, levels = c(0,1), labels = c("Man", "Woman"))

EGC$gender = EGC$gender-1
EGC$gender = factor(EGC$gender, levels = c(0,1), labels = c("Man", "Woman"))

#age
EGC$age = 2014- EGC$birthyear

#################
#### FILTERS ####
#################

#which diseases to include
diseases_inc = c(names(EGC[match("incident_B18",names(EGC)):match("incident_PCOS",names(EGC))]))
#exclude any from previous list?
exclude_diseases = c("incident_i11", "incident_j20","incident_j35", "incident_j45")

#exposures
MB$PB_ratio_int = RankNorm(MB$PB_ratio)
exposures = c("PB_ratio")


#minimum number of cases (for main analysis - does not apply for sub-group analysis)
min_cases = 50

#max number of AB used (max in MB cohort 41)
AB_max = 5

#age exclusion (MB sample age range 23-89)
min_age = 23
max_age = 50

#AB_subgroup: J01CR_n, J01MA_n, J01FA_n
AB_sub = "J01CR_n"
EGC$AB_n=NULL
EGC$AB_n = EGC[,AB_sub]
MB$AB_n=NULL
MB$AB_n = MB[,AB_sub]

#####################
#### FILTERS END ####
#####################

####APPLY FILTERS TO DATASETS###
# diseases to be analyzed
diseases_inc = diseases_inc[!diseases_inc %in% exclude_diseases]

# max N of AB used
EGC = EGC %>% filter(AB_n < AB_max+1)
MB = MB %>% filter(AB_n < AB_max+1)

#age
EGC = EGC %>% filter(age >= min_age) %>% filter(age <= max_age)
MB = MB %>% filter(Age_at_MBsample >= min_age) %>% filter(Age_at_MBsample <= max_age)


####ANALYSIS####

for (j in 1:length(exposures)){
  #j-th exposure
  exposure = exposures[j]
  
  #make empty data-frame for storing the results 
  #each exposure will get it's own table with results for each disease as rows
  res_inc = matrix(NA, nrow = length(diseases_inc), ncol = 10)
  colnames(res_inc) <- c("N", "N_ctrl", "beta", "CI_low", "CI_high", "p_value", "exposure", "beta_exp", "p_value_exp", "group")
  rownames(res_inc) <- diseases_inc
  
  #create a list of prevalent diseases based on the list of incident diseases
  #each prevalent disease is named without prefix "incident_" in our dataframe
  #example: "incident_K58" for indicent disease and "K58" for prevalent disease
  diseases_prev = sub(".*_", "", diseases_inc)
  
  #count the number of incident cases (excluding prevalent cases from dataset beforehand)
  #make a for-loop separately for each disease
  for (i in 1:length(diseases_inc)){
    #exclude the prevalent cases and create a temporary data-frame
    temporary = EGC %>% filter(!EGC[diseases_prev[i]] == 1)
    #count the number of incident cases (when prevalent excluded)
    res_inc[diseases_inc[i], "N"] = table(temporary[diseases_inc[i]])["1"]
    #count the number of controls (when prevalent excluded)
    res_inc[diseases_inc[i], "N_ctrl"] = table(temporary[diseases_inc[i]])["0"]
  }
  
  #make into dataframe
  res_inc = as.data.frame(res_inc)
  
  #filter based on minimum number of incident cases
  res_inc_filtered = res_inc %>% filter(N > min_cases-1)
  
  #create a list of diseases after filtering
  diseases_inc_filtered = rownames(res_inc_filtered)
  #filtered list of prevalent diseases based on incident diseases
  diseases_prev_filtered = sub(".*_", "", diseases_inc_filtered)
  
  # step1 - y
  # step2 - x
  
  # step2 of 2-sample-2-stage least squares (association between instrument and exposure; same for each disease)
  step2 = lm(scale(MB[, match(exposure, names(MB))])~AB_n, data=MB)
  
  #for loop for step1 and combining the two steps for each disease
  for (i in 1:length(diseases_inc_filtered)){
    
    #exclude the prevalent cases and create a temporary data-frame
    temporary = EGC %>% filter(!EGC[diseases_prev_filtered[i]] == 1)
    
    #pursue step1 (association between instrument and outcome (disease)) in the dataset where prevalent cases are excluded
    step1 = glm(temporary[, match(diseases_inc_filtered[i], names(temporary))]~AB_n, family="binomial", data=temporary)
    
    #calculate beta_tsls (two sample two stage least squares) (exposure = microbiome; outcome = disease)
    #beta_exposure_outcome = beta_outcome / beta_exposure
    beta_tsls = step1$coefficients[2]/step2$coefficients[2]
    ##extract information from models to calculate sd for beta_tsls
    #standard deviation of step1 coefficient (instrument on disease beta)
    step1_sd = summary(step1)$coefficients["AB_n","Std. Error"]
    #standard deviation of step2 coefficient (instrument on MB beta)
    step2_sd = summary(step2)$coefficients["AB_n","Std. Error"]
    # step2 coefficient (instrument on MB)
    step2_beta = summary(step2)$coefficients["AB_n","Estimate"]
    
    #sd of tstsls beta based on Pacini and Windmeijer (2016)
    sd_tsls = sqrt((step1_sd**2 + beta_tsls**2 * step2_sd**2)/(step2_beta**2))
    
    #calculate p-value based on z-score
    z_score = beta_tsls/sd_tsls
    p_value = 2*pnorm(abs(z_score), lower.tail = FALSE)
    
    #95% confidence interval for beta_tsls
    CI_low = beta_tsls + 1.96*sd_tsls
    CI_high = beta_tsls - 1.96*sd_tsls
    
    #store values to the results table
    res_inc[diseases_inc_filtered[i], "beta"] = beta_tsls
    res_inc[diseases_inc_filtered[i], "CI_low"] = CI_low
    res_inc[diseases_inc_filtered[i], "CI_high"] = CI_high
    res_inc[diseases_inc_filtered[i], "p_value"] = p_value
    
  }
  
  #name of exposure
  res_inc$exposure = exposure
  #beta from step2 (association between instrument and exposure) - strength of the instrument
  res_inc$beta_exp = summary(step2)$coefficients["AB_n","Estimate"]
  #p-value from step2 (association between instrument and exposure) - strength of the instrument
  res_inc$p_value_exp = summary(step2)$coefficients["AB_n","Pr(>|t|)"]
  #group (currently using all data, but later might want to do sub-analyses by gender (or other groupings))
  res_inc$group = "All"
  #names of diseases (have to be same as row-names)
  res_inc$disease = diseases_inc
  
  
  ####export the results table####
  setwd("C:/Users/Administrator/Documents/01_projektid/09_AB_instrument/03_output/03_PBratio_ABsub/")
  #file-name contains information about exposure and filters (min number of cases, max number of AB, min age, max age)
  write.xlsx(res_inc, paste(paste("res_inc", AB_sub, exposure, min_cases, AB_max, min_age, max_age, sep="_"), ".xlsx", sep=""))
  
  ####PLOT THE RESULT####
  pdf(file = paste(paste("res_inc", AB_sub, exposure, min_cases, AB_max, min_age, max_age, sep="_"), ".pdf", sep=""), 
      height=15, width=6)
  
  
  p2 = ggplot(res_inc %>% filter(group=="All"), aes(x=disease, y=beta, ymin=CI_low, ymax=CI_high),asp=1) + 
    geom_pointrange(fatten=2) + 
    geom_hline(yintercept = 0, linetype=2) +
    coord_flip() +
    theme_bw() +
    theme(legend.position="none",aspect.ratio = 2.5)+
    ylab("beta") +
    xlab(" ") + 
    ggtitle(paste("res_inc", "; exposure: ", exposure, "; min cases: ", min_cases, 
                  "; \n max AB: ", AB_max, "; age min: ", min_age, "; age max: ",
                  max_age, "; \n N(EGC): ", nrow(EGC), "; N(MB): ", nrow(MB),
                  "; \n instrument beta: ", round(res_inc$beta_exp, 3), 
                  ", and p-value: ", round(res_inc$p_value_exp, 3), sep=""))
  print(p2)
  
  dev.off()
  
}




