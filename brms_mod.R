library(tidyverse)
library(brms)
library(bayesplot)


brms_dat <- wave2_dat %>% rename("Y"=increase_2_v_1) %>% mutate(titre=log(`1`))
### Fit single country
form_1 <- bf(Y ~  a + (c-a)/(1+exp(-b*(titre-tm))),
             a ~ 1 ,
             c ~ 1 ,
             b ~ 1,
             tm ~ 1 ,
             nl = TRUE)

priors1 <- c(
  prior(normal(0, 13000), nlpar = "a", lb=0),
  prior(normal(0, 13000), nlpar = "c", lb=0),
  prior(lognormal(2, 5), nlpar = "b", lb=0),
  prior(lognormal(2, 5), nlpar = "tm", lb=0))

mod1 <- brm(form_1, data = brms_dat,
            prior = priors1, seed = 1234,
            family = "bernoulli",
            chains = 4, cores=4, sample_prior = "no")
summary(mod1)
plot(conditional_effects(mod1), points = TRUE)
