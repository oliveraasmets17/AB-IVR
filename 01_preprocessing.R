####PREPROCESSING####
##working directory
setwd("C:/Users/Administrator/Documents/01_projektid/09_AB_instrument")

##libraries
library(dplyr)

##data
phen = readRDS("02_input/Phenotype_data.rds")
medic = readRDS("02_input/MB_medications_data.rds")
mb1 = readRDS("02_input/metaphlan_diversityDf_unfiltered.rds")
mb2 = readRDS("02_input/metaphlan_genus_CLR_filtered1p.rds")
mb3 = readRDS("02_input/metaphlan_species_CLR_filtered1p.rds")
ab_200k = readRDS("02_input/MB_AB_data_200k.rds")
diag_200k = readRDS("02_input/EstBB_diagnosis_200k.rds")
med_200k = readRDS("02_input/q_medicationHistory_AB_200k.rds")


#exclude older data than 10 years
#exclude last 6 months prior to MB sample
#exclude individuals who got medication in last 6 months

######################
####MB SAMPLE DATA####
######################

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

#AB data
medic$J01 <- grepl("\\J01", medic$atc_code)*1

medic_AB = medic %>% filter(J01 == 1)
medic_AB$time_gap = (medic_AB$mb_date - medic_AB$med_date)/365.25
medic_AB$recent = (medic_AB$time_gap < 0.5 & medic_AB$time_gap >= 0)*1


#detect recent AB-users (for excluding)
exclude = unique(medic_AB$skood[medic_AB$recent == 1])
medic_AB = medic_AB %>% filter(time_gap > 0.5)
#exclude as information the AB-usage after MB sampling (but don't exclude the individuals)
medic_AB = medic_AB %>% filter(time_gap <= 10)
#count the number of AB for each individual
AB_usage = medic_AB %>% count(skood)

##merge AB with phen
andmed = merge(AB_usage, phen, by="skood", all.y=T)
#exclude recent AB users
andmed = andmed %>% filter(!skood %in% exclude)

#when n not present, (i.e. no J01 code for that person), assign 0 as AB usage
andmed$n[is.na(andmed$n)==T] = 0
andmed = rename(andmed, AB_n = n)

##merge with MB data
andmed = andmed %>% select(skood, AB_n, BMI, Age_at_MBsample, gender)

andmed = merge(andmed, mb1, by="skood")
andmed = merge(andmed, mb2, by="skood")
andmed = merge(andmed, mb3, by="skood")

####SAVE DATA####
save(andmed, file="02_input/andmed.RData")

##########################
####FULL POPULATION AB####
##########################

#baseline date: 2015-01-01
#exclude older data than 10 years
#exclude last 6 months prior to MB sample
#exclude individuals who got medication in last 6 months


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
AB_usage_200k = ab_200k %>% count(skood)

##merge AB with phen
andmed_200k = merge(AB_usage_200k, diag_200k, by="skood", all.y=T)
#exclude recent AB users
andmed_200k = andmed_200k %>% filter(!skood %in% exclude)

andmed_200k$n[is.na(andmed_200k$n)==T] = 0
andmed_200k = rename(andmed_200k, AB_n = n)

# add data for gender and birthyear
med_200k$birthyear = med_200k$'Person birthYear'
med_200k$gender = med_200k$'Person gender code'

med_200k = med_200k %>% select(skood, birthyear, gender)

andmed_200k = merge(med_200k, andmed_200k, by="skood")


####SAVE DATA####
save(andmed_200k, file="02_input/andmed_200k.RData")


