
##############
# main
##############
source("scripts/utils.R")
# 
# ## run models

results3 <-
  crossing(
    wav = c(2, 3, 4),
    preVar = c("WT"),
    vacc_agnostic_thresh =  c(FALSE, TRUE),
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
      dat %>% filter(age=="adult"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = T,
      diag=T
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
    wav = 4,#c(2, 3, 4),
    vacc_agnostic_thresh = c(TRUE),
    sero_pos_pre = c(TRUE)
  ) %>%
  mutate(preVar = "Omicron",#c("Beta","Delta","Omicron"),
         postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
  filter(!(wav < 4 & !vacc_agnostic_thresh)) %>%
  #filter(sero_pos_pre) %>% 
  #filter(wav==4) %>% 
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="adult"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = T
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


#### Children ----
#wave specific titres
results_children <-
  crossing(
    age = "child",
    wav = c(2,3,4),
    vacc_agnostic_thresh = c(TRUE),
    sero_pos_pre = c(TRUE),
    preVar = c("WT","Beta","Omicron")
  ) %>%
  mutate(postVar = preVar) %>%
  select(age, wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat,
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = F,
      diag = T
    )
  )

(result_tab <- map(results_children,1) %>%
    bind_rows(.id="group") %>% 
    select(-n) %>% 
    left_join(map(results_children,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("Total"),pre=="WT") %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
    arrange(sero_pos_pre,Wave,pre,post,vacc_agnostic_thresh) %>% 
    filter(sero_pos_pre==T) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Age",
                                                  "Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Vaccine agnostic threshold",
                                                  "Only seropositives pre-wave",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Max - min",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "N increased",
                                                  "Doses pre-wave",
                                                  "Doses post-wave",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (2.5% CrI)",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (97.5% CrI)",
                                                  "N")))

write(result_tab, here("results","wt_wave_results_children.html"))


#### Vacc difference ----

results_vacc_diff <-
  crossing(
    wav = c(4),
    preVar = "WT",#c("WT","Omicron"),
    vacc_agnostic_thresh =  c(TRUE),
    sero_pos_pre = c(FALSE),
    vacc_diff=c(T)
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre,vacc_diff) %>%
  filter(!(wav < 4 & !vacc_agnostic_thresh)) %>%
  #filter(sero_pos_pre) %>% 
  #filter(wav==4) %>% 
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="adult"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      vacc_diff = .x$vacc_diff,
      threshold = .01,
      waning = F,
      browsing = T,
      diag=T
    )
  )

#### Vacc difference ----

results_waning <-
  crossing(
    wav = c(4),
    preVar = "WT",#c("WT","Omicron"),
    vacc_agnostic_thresh =  c(TRUE),
    sero_pos_pre = c(FALSE),
    vacc_diff=F,
    waning=T
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre,vacc_diff,waning) %>%
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="adult"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      vacc_diff = .x$vacc_diff,
      threshold = .01,
      waning = .x$waning,
      browsing = T,
      diag=T
    )
  )

#### children ----

results3 <-
  crossing(
    wav = c(2, 3, 4),
    preVar = c("WT"),
    vacc_agnostic_thresh =  c(TRUE),
    sero_pos_pre = c(FALSE, TRUE)
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
  #filter(sero_pos_pre) %>% 
  #filter(wav==4) %>% 
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="child"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      vacc_agnostic_thresh = .x$vacc_agnostic_thresh,
      sero_pos_pre = .x$sero_pos_pre,
      threshold = .01,
      browsing = F,
      diag=T
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
