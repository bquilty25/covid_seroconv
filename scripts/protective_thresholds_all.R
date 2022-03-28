
##############
# main
##############
source("scripts/utils.R")

wv2WB <- calc_wave(2, "WT", "Beta")
wv3WD <- calc_wave(3, "WT", "Delta")
wv4WO <- calc_wave(4, "WT", "Omicron")

wv2BB <- calc_wave(2, "Beta", "Beta")
wv3BD <- calc_wave(3, "Beta", "Delta")
wv4BO <- calc_wave(4, "Beta", "Omicron")

wv4DO <- calc_wave(4, "Delta", "Omicron")

wv4OO <- calc_wave(4, "Omicron", "Omicron")

results <- rbind(wv2WB,wv3WD,wv4WO,
      wv2BB,wv3BD,wv4BO,
      wv4DO,
      wv4OO) 
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

