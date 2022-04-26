
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

results <- rbind(wv2WW$res,
                 wv3WW$res,
                 wv4WW$res) 
saveRDS(results,here("model_output","wave_results.rds"))
(result_tab <- results %>%
  htmlTable::htmlTable(rownames = FALSE, header=c("Wave",
                                                "Variant assessed before wave",
                                                "Variant assessed after wave",
                                                "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)")))


plots <- list(wv2WW$wave_change_plot,
              wv3WW$wave_change_plot,
              wv4WW$wave_change_plot,
              guide_area(),
              wv2WW$wave_plot,
              wv3WW$wave_plot,
              wv4WW$wave_plot) 
  
(wrap_plots(plots,nrow=2,widths = c(1,1,1,0.5))+plot_layout(guides = "collect"))

ggsave(filename = here("results","combined_plot_wt.png"),width = 400,height = 250,units = "mm",dpi=600,bg="white")

write(result_tab, here("results","wave_results.html"))
