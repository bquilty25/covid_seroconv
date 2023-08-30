
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
    vacc_agnostic_thresh = c(FALSE,TRUE),
    sero_pos_pre = c(TRUE)
  ) %>%
  mutate(preVar = "WT",#"Omicron",#c("Beta","Delta","Omicron"),
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
    wav = c(2,3,4,5),
    vacc_agnostic_thresh = c(TRUE),
    sero_pos_pre = c(TRUE),
    preVar = c("WT","Beta","Omicron")
  ) %>%
  mutate(postVar = preVar,
         age ="child") %>%
  # bind_rows(crossing(
  #        wav = c(2, 3, 4),
  #        vacc_agnostic_thresh = c(TRUE),
  #        sero_pos_pre = c(TRUE)
  #      ) %>%
  #          mutate(preVar = "WT",
  #                 postVar = preVar,
  #                 age="child")) %>% 
  select(wav, age, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre) %>%
 # slice(16) %>% 
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
      browsing = T,
      diag=T
    )
  )

(result_tab <- map(results_children,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results_children,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("Total")) %>% 
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

write(result_tab, here("results","wt_wave_results_children.html"))


#### Vacc difference ----

results_vacc_diff <-
  crossing(
    wav = c(2,3,4,5),
    preVar = c("WT","Omicron"),
    sero_pos_pre = c(FALSE),
    vacc_diff=c(F),
    data.frame(age=c("adult","adult","child"), 
             vacc_agnostic_thresh =  c(TRUE,FALSE,TRUE)),
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar, vacc_agnostic_thresh, sero_pos_pre,vacc_diff, age) %>%
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
      vacc_diff = .x$vacc_diff,
      threshold = .01,
      waning = F,
      browsing = F,
      diag=F
    )
  )

(results_vacc_diff_tab <- map(results_vacc_diff,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results_vacc_diff,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
     filter(sero_pos_pre==T,
            pre=="WT",
            !(doses_pre=="0"&age=="child"),
            doses_pre==doses_post|doses_pre=="Overall"
            ) %>% 
    select(age, Wave, pre, post, vacc_agnostic_thresh,  n.x, n_increased, c, a, diff, exp_tm, doses_pre, doses_post, n.y, proportion_protected, count_protected) %>% 
    arrange(age,Wave,pre,post,vacc_agnostic_thresh,doses_pre) %>% 
    mutate(across(c(age, Wave, pre, post, n.x, n_increased, vacc_agnostic_thresh, a, c, diff, exp_tm, proportion_protected, count_protected),unfill_vec)) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Age",
                                                  "Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Threshold including vaccinated",
                                                  "N",
                                                  "N increased",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (median %, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (median %, 95% CrI)",
                                                  "Difference between minimal and maximal (median %, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "Doses pre",
                                                  "Doses post",
                                                  "N in subgroup",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI)",
                                                  "Count of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI)")))

write(results_vacc_diff_tab, here("results","wave_results_vacc_diff_tab.html"))


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

results_children <-
  crossing(
   wav=c(2:5),
   preVar=c("WT","Beta","Delta","Omicron"),
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar) %>%
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="child"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      threshold = .01,
      browsing = F,
      diag=F,
      n_iter=50000
    )
  )

(result_tab <- map(results_children,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results_children,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("0_0","1_1","2_2","Total")) %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
    separate(exp_tm,into = c("x","y"),sep=",") %>% 
    mutate(y=extract_numeric(y),
           y=ifelse(y>2500,">2500",y)) %>% 
    mutate(exp_tm = paste0(x,", ",y,")")) %>% 
    filter(sero_pos_pre==T) %>% 
    select(Wave, pre, post, a, c, diff, exp_tm, n.x, n_increased, doses_pre, doses_post, proportion_protected, count_protected) %>% 
    arrange(Wave,pre,post) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "Difference between minimal and maximal (%, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "N",
                                                  "N increased",
                                                  "Doses pre",
                                                  "Doses post",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI",
                                                  "Count of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI)")))


write(result_tab, here("results","wave_results_children_all.html"))

results_children_for_fig <-
  crossing(
    wav=c(2:5),
    preVar=c("WT")
  ) %>%
  mutate(postVar = preVar) %>%
  select(wav, preVar, postVar) %>%
  rowwise() %>%
  group_split() %>%
  map(
    ~ calc_wave(
      dat %>% filter(age=="child"),
      age=.x$age,
      wav = .x$wav,
      preVar = .x$preVar,
      postVar = .x$postVar,
      threshold = .01,
      browsing = F,
      diag=F,
      n_iter = 50000
    )
  )

(results_tab_children <- map(results_children_for_fig,1) %>%
    bind_rows(.id="group") %>% 
    left_join(map(results_children_for_fig,2) %>% 
                bind_rows(.id="group"),by="group")  %>% 
    select(-group) %>% 
    filter(doses%in%c("0_0","1_1","2_2","Total")) %>% 
    separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% 
    filter(sero_pos_pre==T) %>% 
    select(Wave, pre, post, a, c, diff, exp_tm, n.x, n_increased, doses_pre, doses_post, proportion_protected, count_protected) %>% 
    arrange(Wave,pre,post) %>% 
    htmlTable::htmlTable(rnames = FALSE, header=c("Wave",
                                                  "Variant assessed before wave",
                                                  "Variant assessed after wave",
                                                  "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                  "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                  "Difference between minimal and maximal (%, 95% CrI)",
                                                  "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                  "N",
                                                  "N increased",
                                                  "Doses pre",
                                                  "Doses post",
                                                  "Proportion of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI",
                                                  "Count of seropositives with pre-wave antibody titres higher than threshold (median, 95% CrI)")))

write(results_tab_children, here("results","wave_results_children_wt.html"))
