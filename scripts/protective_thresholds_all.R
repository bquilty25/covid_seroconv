
##############
# main
##############
source("scripts/utils.R")

## run models
#wv2WB <- calc_wave(dat=datw4, wav=2, preVar="WT", postVar="Beta", threshold=.1)
#wv3WD <- calc_wave(datw4, 3, "WT", "Delta")
#wv4WO <- calc_wave(datw4, 4, "WT", "Omicron")
#wv2BB <- calc_wave(datw4, 2, "Beta", "Beta")
#wv3BD <- calc_wave(datw4, 3, "Beta", "Delta")
#wv4BO <- calc_wave(datw4, 4, "Beta", "Omicron")
#wv4DO <- calc_wave(datw4, 4, "Delta", "Omicron")
#wv4OO <- calc_wave(datw4, 4, "Omicron", "Omicron")

wv2WW <- calc_wave(datw4, 2, "WT", "WT",.1)
wv3WW <- calc_wave(datw4, 3, "WT", "WT",.1)
wv4WW <- calc_wave(datw4, 4, "WT", "WT",.1)

results <- rbind(wv2WW,
                 wv3WW,
                 wv4WW) 
saveRDS(results,here("model_output","wave_results.rds"))
results %>%
  htmlTable::htmlTable(rownames = FALSE, header=c("Wave",
                                                "Variant assessed before wave",
                                                "Variant assessed after wave",
                                                "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)"))


