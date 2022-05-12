
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

#WT

wv2WW <- calc_wave(datw4, 2, "WT", "WT",.01)
wv3WW <- calc_wave(datw4, 3, "WT", "WT",.01)
wv4WW <- calc_wave(datw4, 4, "WT", "WT",.01)

results <- rbind(wv2WW$res,
                 wv3WW$res,
                 wv4WW$res) 

saveRDS(results,here("model_output","wave_results.rds"))

(result_tab <- results %>%
  htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
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

write(result_tab, here("results","wt_wave_results.html"))

#proportion seropositive pre-wave protected based on thresholds
datw4 %>% drop_na(value) %>% 
  left_join(results %>% 
              select(Wave,a_0.5) %>% 
              mutate(pre_wave=Wave-1),
            by=c("wave"="pre_wave")) %>% 
  mutate(a_0.5=as.integer(word(a_0.5))) %>% 
  filter(value>1.09) %>% 
  group_by(wave) %>% 
  summarise(n=n(),
            protected_50=sum(value>a_0.5)/n)

#WT vs. wave specific variants

wv2beta <- calc_wave(datw4, 2, "WT", "Beta",.1)
wv3delta <- calc_wave(datw4, 3, "WT", "Delta",.1)
wv4omi <- calc_wave(datw4, 4, "WT", "Omicron",.1,browsing = T)

results <- rbind(wv2beta$res,
                 wv3delta$res,
                 wv4omi$res) 

saveRDS(results,here("model_output","wave_results_wave_specific.rds"))

(result_tab <- results %>%
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                    "Variant assessed before wave",
                                                    "Variant assessed after wave",
                                                    "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                    "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                    "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                    "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                    "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)")))

plots <- list(wv2beta$wave_change_plot,
              wv3delta$wave_change_plot,
              wv4omi$wave_change_plot,
              guide_area(),
              wv2beta$wave_plot,
              wv3delta$wave_plot,
              wv4omi$wave_plot) 

(wrap_plots(plots,nrow=2,widths = c(1,1,1,0.5))+plot_layout(guides = "collect"))

ggsave(filename = here("results","combined_plot_wave_specific.png"),width = 400,height = 250,units = "mm",dpi=600,bg="white")

write(result_tab, here("results","wave_specific_results.html"))


# Vaccinated
#WT

wv2WW_vacc <- calc_wave(datw4, vacc=T, 2, "WT", "WT",.01)
wv3WW_vacc <- calc_wave(datw4, vacc=T, 3, "WT", "WT",.01,browsing = F)
wv4WW_vacc <- calc_wave(datw4, vacc=T, 4, "WT", "WT",.01,browsing = F)

(result_tab <- wv4WW_vacc$res %>%
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)")))

plots <- list(wv2WW_vacc$wave_change_plot,
              wv3WW_vacc$wave_change_plot,
              wv4WW_vacc$wave_change_plot,
              guide_area(),
              wv2WW_vacc$wave_plot,
              wv3WW_vacc$wave_plot,
              wv4WW_vacc$wave_plot) 

(wrap_plots(plots,nrow=2,widths = c(1,1,1,0.5))+plot_layout(guides = "collect"))

ggsave(filename = here("results","combined_plot_wt.png"),width = 400,height = 250,units = "mm",dpi=600,bg="white")

write(result_tab, here("results","wt_wave_results.html"))

#proportion seropositive pre-wave protected based on thresholds
datw4_vacc %>% drop_na(value) %>% 
  left_join(results %>% 
              select(Wave,a_0.5) %>% 
              mutate(pre_wave=Wave-1),
            by=c("wave"="pre_wave")) %>% 
  mutate(a_0.5=as.integer(word(a_0.5))) %>% 
  filter(value>1.09) %>% 
  group_by(wave) %>% 
  summarise(n=n(),
            protected_50=sum(value>a_0.5)/n)
