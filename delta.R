pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra)
remotes::install_github("njtierney/mmcc")
library(mmcc)

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Load and clean data
dat <- read_xlsx("SARS_CoV2_data_3_waves_15dec2021.xlsx",sheet = 2) %>% 
  pivot_longer(CoV2S_1w_m:delta_3w_m) %>% 
  mutate(variant=case_when(str_detect(name,"CoV2")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta"),
         variant=fct_relevel(variant,"WT","Beta","Delta"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name)

dat %>%
  ggplot(aes(x=name,y=value,group=pid_child))+
  geom_point(alpha=0.2)+
  geom_path(alpha=0.2)+
  facet_wrap(~variant,scale="free_x",labeller = labeller(`fct_rev(age)`=Hmisc::capitalize))+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))+
  scale_colour_brewer(type="qual",guide=F,direction=-1)

ggsave("change_new.png",width=400,height=150,units="mm",dpi=600,bg="white")

dat %>% pivot_wider(values_from=value,names_from = wave)
  
