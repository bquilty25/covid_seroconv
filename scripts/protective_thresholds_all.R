
##############
# main
##############
source("scripts/utils.R")
# 
# ## run models
# #wv2WB <- calc_wave(dat=datw4, wav=2, preVar="WT", postVar="Beta", threshold=.1)
# #wv3WD <- calc_wave(datw4, 3, "WT", "Delta")
# #wv4WO <- calc_wave(datw4, 4, "WT", "Omicron")
# #wv2BB <- calc_wave(datw4, 2, "Beta", "Beta")
# #wv3BD <- calc_wave(datw4, 3, "Beta", "Delta")
# #wv4BO <- calc_wave(datw4, 4, "Beta", "Omicron")
# #wv4DO <- calc_wave(datw4, 4, "Delta", "Omicron")
# #wv4OO <- calc_wave(datw4, 4, "Omicron", "Omicron")
# 
# #WT
# 
# wv2WW <- calc_wave(datw4, 2, "WT", "WT",.01)
# wv3WW <- calc_wave(datw4, 3, "WT", "WT",.01)
# wv4WW <- calc_wave(datw4, 4, "WT", "WT",.01)
# 
# results <- rbind(wv2WW$res,
#                  wv3WW$res,
#                  wv4WW$res) 
# 
# saveRDS(results,here("model_output","wave_results.rds"))
# 
# (result_tab <- results %>%
#   htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
#                                                 "Variant assessed before wave",
#                                                 "Variant assessed after wave",
#                                                 "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
#                                                 "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
#                                                 "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
#                                                 "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
#                                                 "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)")))
# 
# plots <- list(wv2WW$wave_change_plot,
#               wv3WW$wave_change_plot,
#               wv4WW$wave_change_plot,
#               guide_area(),
#               wv2WW$wave_plot,
#               wv3WW$wave_plot,
#               wv4WW$wave_plot) 
#   
# (wrap_plots(plots,nrow=2,widths = c(1,1,1,0.5))+plot_layout(guides = "collect"))
# 
# ggsave(filename = here("results","combined_plot_wt.png"),width = 400,height = 250,units = "mm",dpi=600,bg="white")
# 
# write(result_tab, here("results","wt_wave_results.html"))
# 
# #proportion seropositive pre-wave protected based on thresholds
# datw4 %>% drop_na(value) %>% 
#   left_join(results %>% 
#               select(Wave,a_0.5) %>% 
#               mutate(pre_wave=Wave-1),
#             by=c("wave"="pre_wave")) %>% 
#   mutate(a_0.5=as.integer(word(a_0.5))) %>% 
#   filter(value>1.09) %>% 
#   group_by(wave) %>% 
#   summarise(n=n(),
#             protected_50=sum(value>a_0.5)/n)
# 
# #WT vs. wave specific variants
# 
# wv2beta <- calc_wave(datw4, 2, "WT", "Beta",.1)
# wv3delta <- calc_wave(datw4, 3, "WT", "Delta",.1)
# wv4omi <- calc_wave(datw4, 4, "WT", "Omicron",.1)
# 
# results2 <- rbind(wv2beta$res,
#                  wv3delta$res,
#                  wv4omi$res) 
# 
# saveRDS(results,here("model_output","wave_results_wave_specific.rds"))
# 
# (result_tab <- results %>%
#     htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
#                                                     "Variant assessed before wave",
#                                                     "Variant assessed after wave",
#                                                     "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
#                                                     "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
#                                                     "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
#                                                     "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
#                                                     "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)")))
# 
# plots <- list(wv2beta$wave_change_plot,
#               wv3delta$wave_change_plot,
#               wv4omi$wave_change_plot,
#               guide_area(),
#               wv2beta$wave_plot,
#               wv3delta$wave_plot,
#               wv4omi$wave_plot) 
# 
# (wrap_plots(plots,nrow=2,widths = c(1,1,1,0.5))+plot_layout(guides = "collect"))
# 
# ggsave(filename = here("results","combined_plot_wave_specific.png"),width = 400,height = 250,units = "mm",dpi=600,bg="white")
# 
# write(result_tab, here("results","wave_specific_results.html"))
# 
# 
# # Vaccinated
# #WT
# 
# wv2WW_vacc <- calc_wave(datw4, vacc=T, 2, "WT", "WT",.01)
# wv3WW_vacc <- calc_wave(datw4, vacc=T, 3, "WT", "WT",.01,browsing = F)
# wv4WW_vacc <- calc_wave(datw4, 4, "WT", "WT",.01,browsing = T)
# 

results3 <-
  crossing(
    wav = c(2, 3, 4),
    preVar = c("WT"),
    vacc_agnostic_thresh = c(FALSE, TRUE),
    sero_pos_pre = c(FALSE, TRUE)
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
  filter(!(wav < 4 & !vacc_agnostic_thresh)) %>%
  #filter(sero_pos_pre) %>% 
  #filter(wav==4) %>% 
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      datw4,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = F
    )
  )

(result_tab <- map(results3,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results3,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("0_0","1_1","2_2","Total")) %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
    arrange(sero_pos_pre,Wave,pre,post,vacc_agnostic_thresh) %>% 
    filter(sero_pos_pre==T) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Vaccine agnostic threshold",
                                                  "Only seropositives pre-wave",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "N",
                                                  "N increased",
                                                  "Doses pre-wave",
                                                  "Doses post-wave",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (2.5% CrI)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (97.5% CrI)",
                                                  "N in subgroup")))


write(result_tab, here("results","wt_wave_results.html"))

#wave specific titres
results4 <-
  crossing(
    wav = c(2, 3, 4),
    vacc_agnostic_thresh = c(TRUE),
    sero_pos_pre = c(TRUE)
  ) %>%
  mutate(preVar = c("Beta","Delta","Omicron"),
         postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
  filter(!(wav < 4 & !vacc_agnostic_thresh)) %>%
  #filter(sero_pos_pre) %>% 
  #filter(wav==4) %>% 
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      datw4,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = F
    )
  )

(result_tab <- map(results4,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results4,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("0_0","1_1","2_2","Total")) %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
    arrange(sero_pos_pre,Wave,pre,post,vacc_agnostic_thresh) %>% 
    filter(sero_pos_pre==T) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Vaccine agnostic threshold",
                                                  "Only seropositives pre-wave",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "N",
                                                  "N increased",
                                                  "Doses pre-wave",
                                                  "Doses post-wave",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (2.5% CrI)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (97.5% CrI)",
                                                  "N in subgroup")))
