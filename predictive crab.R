library(lme4)
library(tidyverse)
source("ggplot_paper_theme.R") #for custom plot themes
set.seed(12345)

#LOAD IN DATA------------------------------------------------------------------
#load in experimental data used to make the model
mycrab <- read.csv("greencrabforaging.csv") %>% 
  filter(Clam == "VC", Sed == "no") %>% 
  #standardize continuous predictors and calculate proportions and ratios
  mutate(clam_alive = Density - Clam_Consumed,
         stand_density = c(scale(Density)),
         stand_cw = c(scale(CW)),
         stand_claw = c(scale(Claw_width)),
         prop_eaten = Clam_Consumed/Density,
         claw_ratio = Claw_width/CW,
         stand_claw_ratio = c(scale(claw_ratio)))

#load in the field data on sex differences in real crab populations
m_to_f <- read_csv("crsdata.csv") %>% 
  transmute(site = Site, 
            day = day,
            month = month,
            year = year,
            male_total = `Male total (#)`,
            female_total = `Female total (#)`,
            total = male_total + female_total,
            prop_female = female_total/(female_total + male_total)) %>% 
  unite("unique_id",c(site, day, month, year), remove = FALSE) %>% 
  #remove all the empty rows from the weirdly formatted excel sheet
  filter(!is.na(male_total))
#the warnings here are fine and just refer to the issue with column names
#in the original dataset we were given

#first day only from each trap set
m_to_f_first <- m_to_f %>% 
  group_by(site, month, year) %>% 
  arrange(day) %>% 
  distinct(site, month, year, .keep_all = TRUE)

m_to_f_accum <- m_to_f %>% 
  group_by(site, month, year) %>% 
  summarize(male_accum = sum(male_total),
            female_accum = sum(female_total),
            prop_female_accum = female_accum/(female_accum + male_accum),
            total_accum = male_accum + female_accum,
            unique_id = first(unique_id)) %>% 
  ungroup()

crab_sizes <- read_csv("crab_population.csv") %>% 
  group_by(sex) %>% 
  summarize(n = n(),
            mean_claw = mean(claw),
            mean_cw = mean(cw),
            sd_claw = sd(claw),
            sd_cw = sd(cw),
            se_claw = sd_claw/sqrt(n),
            se_cw = sd_cw/sqrt(n),
            ci_claw = se_claw*1.96,
            ci_cw = se_cw*1.96) %>% 
  ungroup()

#MODELS-----------------------------------------------------------------------
#this is the model for carapace width
mod_final <- glm(prop_eaten ~ stand_density + stand_cw + Sex + Sex:stand_cw, 
           data = mycrab, 
           weights = Density, 
           family = binomial(link = "logit"))
mod_final_summ <- summary(mod_final)

#and this is the model for claw width
mod_final2 <- glm(prop_eaten ~ stand_density + stand_claw + Sex + 
                    Sex:stand_claw, 
                 data = mycrab, 
                 weights = Density, 
                 family = binomial(link = "logit"))
mod_final2_summ <- summary(mod_final2)
#CARAPACE WIDTH BOOTSTRAP-------------------------------------------------------
#and we can take a look at the model estimates for the proportion of clams that
#females and males will eat at the mean available density, after accounting for
#differences in claw width
female_rate <- plogis(mod_final_summ$coefficients[1])
male_rate <- plogis(mod_final_summ$coefficients[1] + 
                      mod_final_summ$coefficients[4])

#going to rerun the model 1000 times with resampled data so we can 
#extract the prediction for each iteration of how much females are eating 
#compared to males
boot_model <- function(dataframe){
  #this is probably a janky way of doing this but it works!
  #sample from the total number of rows in the df
  xboot <- sample(seq(from = 1, to = nrow(dataframe), by = 1), 
                  replace = TRUE) %>% 
    as_tibble()
  #then extract the rows associated from those row numbers in the df
  boot_df <- dataframe %>% 
    rownames_to_column() %>% 
    mutate(rowname = as.numeric(rowname)) %>% 
    right_join(xboot, by = c("rowname" = "value"))
  #run a model on the resampled data
  mod <- glm(prop_eaten ~ stand_density + stand_cw + Sex + Sex:stand_cw, 
                data = boot_df, 
                weights = Density, 
                family = binomial(link = "logit"))
  #extract the info from the output of the model
  modsummary <- summary(mod)
  #backtransform the model predicted "proportion of available clams eaten" for
  #each sex
  female_rate <- plogis(modsummary$coefficients[1])
  male_rate <- plogis(modsummary$coefficients[1] + modsummary$coefficients[4])
  #extract the percent difference of females relative to males
  percent_diff <- (female_rate - male_rate)/male_rate
  c(male_rate, female_rate, percent_diff)
}

#now run the function 10000 times to bootstrap how much less we predict females
#are eating compared to males - 1000 seemed too small to get a smooth 
#distribution
boot_diff <- replicate(1000, boot_model(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3)
mean_male <- mean(boot_diff$male_rate)
mean_female <- mean(boot_diff$female_rate)

#take a smaller sample for plotting
sample_boot <- replicate(1000, boot_model(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3) %>% 
  #now compare each row to the mean male rate
  mutate(male_percent = (male_rate - mean_male)/mean_male,
         female_percent = (female_rate - mean_male)/mean_male,
         male_percent2 = male_rate,
         female_percent2 = (female_rate - male_rate)/male_rate)

ggplot(boot_diff, aes(percent_diff)) +
  geom_histogram()

upper_bound <- quantile(boot_diff$percent_diff, 0.975)
lower_bound <- quantile(boot_diff$percent_diff, 0.025)

ci_df <- NULL
for(i in 1:1000){
  male_rate <- sample_boot$male_rate[i]
  female_rate <- sample_boot$female_rate[i]
  df <- tibble(prop_female = seq(from = 0, to = 1, by = 0.01),
               prop_eaten_male = (male_rate - mean_male)/mean_male,
               prop_eaten_female = (female_rate - mean_male)/mean_male,
               prop_eaten_male2 = (male_rate-male_rate)/male_rate,
               prop_eaten_female2 = (female_rate - male_rate)/male_rate,
               #this first delta metric compares both the iteration-specific 
               #male and female rates to the mean male-only population
               #(i.e., it incorporates uncertainty around both the male and 
               #female estimates)
               delta = (1 - prop_female) * prop_eaten_male + 
                 prop_female * prop_eaten_female,
               #this second delta metric 
               delta2 = (1 - prop_female) * prop_eaten_male2 + 
                 prop_female * prop_eaten_female2,
               iteration = i)
  ci_df <- bind_rows(ci_df, df)
}

summ_cis <- ci_df %>% 
  group_by(prop_female) %>% 
  summarize(upper = quantile(delta, 0.975),
            lower = quantile(delta, 0.025),
            median = quantile(delta, 0.5),
            upper2 = quantile(delta2, 0.975),
            lower2 = quantile(delta2, 0.025),
            median2 = quantile(delta2, 0.5))

ggplot(summ_cis, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = median)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population")

ggplot(summ_cis, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.2) +
  geom_line(aes(y = median2)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population")


#now plot the real populations on the hypothetical figure
real_props <- tibble(prop_female = m_to_f$prop_female,
                     unique_id = m_to_f$unique_id,
                     prop_eaten_male = (mean_male - mean_male)/mean_male,
                     prop_eaten_female = (mean_female - mean_male)/mean_male,
                     delta = (1 - prop_female) * prop_eaten_male + 
                        prop_female * prop_eaten_female)

ggplot(summ_cis, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = median)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population") +
  geom_point(data = real_props, aes(y = delta))

real_props_first <- tibble(prop_female = m_to_f_first$prop_female,
                     unique_id = m_to_f_first$unique_id,
                     site = m_to_f_first$site,
                     prop_eaten_male = (mean_male - mean_male)/mean_male,
                     prop_eaten_female = (mean_female - mean_male)/mean_male,
                     delta = (1 - prop_female) * prop_eaten_male + 
                       prop_female * prop_eaten_female)

ggplot(summ_cis, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = median)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population") +
  geom_point(data = real_props_first, aes(y = delta)) +
  ggrepel::geom_label_repel(data = real_props_first,
                            aes(y = delta, label = site, 
                                colour = site),
                            max.overlaps = 100) +
  theme(legend.position = "none")

ggsave("test.png", width = 8, height = 6)




#CLAW WIDTH BOOTSTRAP----------------------------------------------------------
#we can do this one of two ways:
#either assuming a hypothetical population has males and females of the same
#sizes or assuming the male and female size distributions are reflective of
#the mean sizes actually observed

female_rate_claw <- plogis(mod_final2_summ$coefficients[1])
male_rate_claw <- plogis(mod_final2_summ$coefficients[1] + 
                      mod_final2_summ$coefficients[4])

#going to rerun the model 1000 times with resampled data so we can 
#extract the prediction for each iteration of how much females are eating 
#compared to males
boot_model_claw <- function(dataframe){
  #this is probably a janky way of doing this but it works!
  #sample from the total number of rows in the df
  xboot <- sample(seq(from = 1, to = nrow(dataframe), by = 1), 
                  replace = TRUE) %>% 
    as_tibble()
  #then extract the rows associated from those row numbers in the df
  boot_df <- dataframe %>% 
    rownames_to_column() %>% 
    mutate(rowname = as.numeric(rowname)) %>% 
    right_join(xboot, by = c("rowname" = "value"))
  #run a model on the resampled data
  mod <- glm(prop_eaten ~ stand_density + stand_claw + Sex + Sex:stand_claw, 
             data = boot_df, 
             weights = Density, 
             family = binomial(link = "logit"))
  #extract the info from the output of the model
  modsummary <- summary(mod)
  #backtransform the model predicted "proportion of available clams eaten" for
  #each sex
  female_rate <- plogis(modsummary$coefficients[1])
  male_rate <- plogis(modsummary$coefficients[1] + modsummary$coefficients[4])
  #extract the percent difference of females relative to males
  percent_diff <- (female_rate - male_rate)/male_rate
  c(male_rate, female_rate, percent_diff)
}

#now run the function 10000 times to bootstrap how much less we predict females
#are eating compared to males 
boot_diff_claw <- replicate(1000, boot_model_claw(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3)
mean_male_claw <- mean(boot_diff_claw$male_rate)
mean_female_claw <- mean(boot_diff_claw$female_rate)

#for plotting
sample_boot_claw <- replicate(1000, boot_model_claw(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3) %>% 
  #now compare each row to the mean male rate
  mutate(male_percent = (male_rate - mean_male_claw)/mean_male_claw,
         female_percent = (female_rate - mean_male_claw)/mean_male_claw,
         male_percent2 = male_rate,
         female_percent2 = (female_rate - male_rate)/male_rate)

ggplot(boot_diff_claw, aes(percent_diff)) +
  geom_histogram()

upper_bound_claw <- quantile(boot_diff_claw$percent_diff, 0.975)
lower_bound_claw <- quantile(boot_diff_claw$percent_diff, 0.025)

ci_df_claw <- NULL
for(i in 1:1000){
  male_rate <- sample_boot_claw$male_rate[i]
  female_rate <- sample_boot_claw$female_rate[i]
  df <- tibble(prop_female = seq(from = 0, to = 1, by = 0.01),
               prop_eaten_male = (male_rate - mean_male_claw)/mean_male_claw,
               prop_eaten_female = (female_rate - mean_male_claw)/mean_male_claw,
               prop_eaten_male2 = (male_rate-male_rate)/male_rate,
               prop_eaten_female2 = (female_rate - male_rate)/male_rate,
               #this first delta metric compares both the iteration-specific 
               #male and female rates to the mean male-only population
               #(i.e., it incorporates uncertainty around both the male and 
               #female estimates)
               delta = (1 - prop_female) * prop_eaten_male + 
                 prop_female * prop_eaten_female,
               #this second delta metric 
               delta2 = (1 - prop_female) * prop_eaten_male2 + 
                 prop_female * prop_eaten_female2,
               iteration = i)
  ci_df_claw <- bind_rows(ci_df_claw, df)
}

summ_cis_claw <- ci_df_claw %>% 
  group_by(prop_female) %>% 
  summarize(upper = quantile(delta, 0.975),
            lower = quantile(delta, 0.025),
            median = quantile(delta, 0.5),
            upper2 = quantile(delta2, 0.975),
            lower2 = quantile(delta2, 0.025),
            median2 = quantile(delta2, 0.5))

ggplot(summ_cis_claw, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = median)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population")

ggplot(summ_cis_claw, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.2) +
  geom_line(aes(y = median2)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Proportional difference in number of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population")


#CLAW WIDTH WITH SIZE CORRECTION----------------------------------------------
#first take the mean claw size for males and females from a real population and 
#convert it to the standardized form used in the model
stand_mean <- attributes(scale(mycrab$Claw_width))$`scaled:center`
stand_sd <- attributes(scale(mycrab$Claw_width))$`scaled:scale`
female_claw <- (as.data.frame(crab_sizes)[1,3] - stand_mean)/stand_sd
male_claw <- (as.data.frame(crab_sizes)[2,3] - stand_mean)/stand_sd

#calculate the estimated effect for a female and a male crab at their 
#respective mean sizes using the model coefficients, and then use the plogis
#function to backtransform these estimates out of logit-link space
mean_female_claw_cor  <- plogis(mod_final2_summ$coefficients[1] + 
                                 mod_final2_summ$coefficients[3] * female_claw)
mean_male_claw_cor  <- plogis(mod_final2_summ$coefficients[1] + 
                           mod_final2_summ$coefficients[4] + 
                             mod_final2_summ$coefficients[3] * male_claw + 
                             mod_final2_summ$coefficients[5] * male_claw)

#going to rerun the model 1000 times with resampled data so we can 
#extract the prediction for each iteration of how much females are eating 
#compared to males
boot_model_claw_cor <- function(dataframe){
  #this is probably a janky way of doing this but it works!
  #sample from the total number of rows in the df
  xboot <- sample(seq(from = 1, to = nrow(dataframe), by = 1), 
                  replace = TRUE) %>% 
    as_tibble()
  #then extract the rows associated from those row numbers in the df
  boot_df <- dataframe %>% 
    rownames_to_column() %>% 
    mutate(rowname = as.numeric(rowname)) %>% 
    right_join(xboot, by = c("rowname" = "value"))
  #run a model on the resampled data
  mod <- glm(prop_eaten ~ stand_density + stand_claw + Sex + Sex:stand_claw, 
             data = boot_df, 
             weights = Density, 
             family = binomial(link = "logit"))
  #extract the info from the output of the model
  modsummary <- summary(mod)
  #backtransform the model predicted "proportion of available clams eaten" for
  #each sex
  female_rate <- plogis(modsummary$coefficients[1] + 
                          modsummary$coefficients[3] * female_claw)
  male_rate <-plogis(modsummary$coefficients[1] + 
                       modsummary$coefficients[4] + 
                       modsummary$coefficients[3] * male_claw + 
                       modsummary$coefficients[5] * male_claw)
  #extract the percent difference of females relative to males
  percent_diff <- (female_rate - male_rate)/male_rate
  c(male_rate, female_rate, percent_diff)
}

#now run the function 2000 times to bootstrap how much less we predict females
#are eating compared to males - 1000 seemed too small to get a smooth 
#distribution
sample_boot_claw_cor <- replicate(2000, boot_model_claw_cor(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3) %>% 
  #now compare each row to the male rate from the original model
  mutate(male_percent = (male_rate - mean_male_claw_cor)/mean_male_claw_cor,
         female_percent = (female_rate - mean_male_claw_cor)/mean_male_claw_cor,
         #and each row to the male in the bootstrap iteration
         male_percent2 = male_rate,
         female_percent2 = (female_rate - male_rate)/male_rate)

ggplot(sample_boot_claw_cor, aes(percent_diff)) +
  geom_histogram()

upper_bound_claw_cor <- quantile(sample_boot_claw_cor$percent_diff, 0.975)
lower_bound_claw_cor <- quantile(sample_boot_claw_cor$percent_diff, 0.025)

ci_df_claw_cor <- NULL
for(i in 1:2000){
  male_rate <- sample_boot_claw_cor$male_rate[i]
  female_rate <- sample_boot_claw_cor$female_rate[i]
  df <- tibble(prop_female = seq(from = 0, to = 1, by = 0.01),
               prop_eaten_male = (male_rate - mean_male_claw_cor)/mean_male_claw_cor,
               prop_eaten_female = (female_rate - mean_male_claw_cor)/mean_male_claw_cor,
               prop_eaten_male2 = (male_rate-male_rate)/male_rate,
               prop_eaten_female2 = (female_rate - male_rate)/male_rate,
               #this first delta metric compares both the iteration-specific 
               #male and female rates to the mean male-only population
               #(i.e., it incorporates uncertainty around both the male and 
               #female estimates)
               delta = (1 - prop_female) * prop_eaten_male + 
                 prop_female * prop_eaten_female,
               #this second delta metric only incorporates uncertainty around 
               #the female estimate
               delta2 = (1 - prop_female) * prop_eaten_male2 + 
                 prop_female * prop_eaten_female2,
               iteration = i)
  ci_df_claw_cor <- bind_rows(ci_df_claw_cor, df)
}

summ_cis_claw_cor <- ci_df_claw_cor %>% 
  group_by(prop_female) %>% 
  summarize(upper = quantile(delta, 0.975),
            lower = quantile(delta, 0.025),
            median = quantile(delta, 0.5),
            mean = mean(delta),
            upper2 = quantile(delta2, 0.975),
            lower2 = quantile(delta2, 0.025),
            median2 = quantile(delta2, 0.5),
            mean2 = mean(delta2))

og_model_outputs <- tibble(prop_female = seq(from = 0, to = 1, by = 0.001),
                           delta = (1 - prop_female) * (mean_male_claw_cor - mean_male_claw_cor)/mean_male_claw_cor + 
                             prop_female * (mean_female_claw_cor-mean_male_claw_cor)/mean_male_claw_cor)

#now plot the real populations on the hypothetical figure
real_props_claw_cor <- tibble(prop_female = m_to_f_accum$prop_female_accum,
                     unique_id = m_to_f_accum$unique_id,
                     site = m_to_f_accum$site,
                     #set to -Inf so they'll be right along the axis
                     delta = -Inf)


ggplot(summ_cis_claw_cor, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(data = og_model_outputs, aes(y = delta)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Predicted difference in the proportion of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population") +
  scale_shape_identity() +
  geom_point(data = real_props_claw_cor, aes(y = delta), 
             shape = 108, size = 5) +
  theme(legend.position = "none")

#ggsave("Figure6.pdf", height = 5, width = 7)
#ggsave("Figure6.png", height = 5, width = 7)

#SUPPLEMENT: CLAW HEIGHT WITH SIZE CORRECTION EXTREME CIS--------------
#do the exact same thing as the section above except instead of using the mean
#size estimates for males and females, we'll use the upper end of the female 
#CI and lower end of the male CI to see what the smallest possible effect
#would be
stand_mean <- attributes(scale(mycrab$Claw_width))$`scaled:center`
stand_sd <- attributes(scale(mycrab$Claw_width))$`scaled:scale`
female_claw_upper <- (as.data.frame(crab_sizes)[1,3] + 
                  as.data.frame(crab_sizes)[1,9] - stand_mean)/stand_sd
male_claw_lower <- (as.data.frame(crab_sizes)[2,3] - 
                as.data.frame(crab_sizes)[2,9] - stand_mean)/stand_sd

#calculate the estimated effect for a female and a male crab at their 
#respective mean sizes using the model coefficients, and then use the plogis
#function to backtransform these estimates out of logit-link space
female_rate_claw_cor_upper <- plogis(mod_final2_summ$coefficients[1] + 
                                 mod_final2_summ$coefficients[3] * female_claw_upper)
male_rate_claw_cor_lower <- plogis(mod_final2_summ$coefficients[1] + 
                               mod_final2_summ$coefficients[4] + 
                               mod_final2_summ$coefficients[3] * male_claw_lower + 
                               mod_final2_summ$coefficients[5] * male_claw_lower)

#going to rerun the model 1000 times with resampled data so we can 
#extract the prediction for each iteration of how much females are eating 
#compared to males
boot_model_claw_cor_extreme <- function(dataframe){
  #this is probably a janky way of doing this but it works!
  #sample from the total number of rows in the df
  xboot <- sample(seq(from = 1, to = nrow(dataframe), by = 1), 
                  replace = TRUE) %>% 
    as_tibble()
  #then extract the rows associated from those row numbers in the df
  boot_df <- dataframe %>% 
    rownames_to_column() %>% 
    mutate(rowname = as.numeric(rowname)) %>% 
    right_join(xboot, by = c("rowname" = "value"))
  #run a model on the resampled data
  mod <- glm(prop_eaten ~ stand_density + stand_claw + Sex + Sex:stand_claw, 
             data = boot_df, 
             weights = Density, 
             family = binomial(link = "logit"))
  #extract the info from the output of the model
  modsummary <- summary(mod)
  #backtransform the model predicted "proportion of available clams eaten" for
  #each sex
  female_rate <- plogis(modsummary$coefficients[1] + 
                          modsummary$coefficients[3] * female_claw_upper)
  male_rate <-plogis(modsummary$coefficients[1] + 
                       modsummary$coefficients[4] + 
                       modsummary$coefficients[3] * male_claw_lower + 
                       modsummary$coefficients[5] * male_claw_lower)
  #extract the percent difference of females relative to males
  percent_diff <- (female_rate - male_rate)/male_rate
  c(male_rate, female_rate, percent_diff)
}

#now run the function 10000 times to bootstrap how much less we predict females
#are eating compared to males - 1000 seemed too small to get a smooth 
#distribution
boot_diff_claw_cor_extreme <- replicate(1000, 
                                        boot_model_claw_cor_extreme(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3)
mean_male_claw_cor_extreme <- mean(boot_diff_claw_cor_extreme$male_rate)
mean_female_claw_cor_extreme <- mean(boot_diff_claw_cor_extreme$female_rate)

sample_boot_claw_cor_extreme <- replicate(1000, 
                                          boot_model_claw_cor_extreme(mycrab)) %>% 
  t() %>% 
  as_tibble() %>% 
  transmute(male_rate = V1,
            female_rate = V2,
            percent_diff = V3) %>% 
  #now compare each row to the mean male rate
  mutate(male_percent = (male_rate - mean_male_claw_cor_extreme)/
           mean_male_claw_cor_extreme,
         female_percent = (female_rate - mean_male_claw_cor_extreme)/
           mean_male_claw_cor_extreme,
         male_percent2 = male_rate,
         female_percent2 = (female_rate - male_rate)/male_rate)

ggplot(boot_diff_claw_cor_extreme, aes(percent_diff)) +
  geom_histogram()

upper_bound_claw_cor_extreme <- quantile(boot_diff_claw_cor_extreme$percent_diff, 0.975)
lower_bound_claw_cor_extreme <- quantile(boot_diff_claw_cor_extreme$percent_diff, 0.025)

ci_df_claw_cor_extreme <- NULL
for(i in 1:1000){
  male_rate <- sample_boot_claw_cor_extreme$male_rate[i]
  female_rate <- sample_boot_claw_cor_extreme$female_rate[i]
  df <- tibble(prop_female = seq(from = 0, to = 1, by = 0.0001),
               prop_eaten_male = (male_rate - mean_male_claw_cor_extreme)/mean_male_claw_cor_extreme,
               prop_eaten_female = (female_rate - mean_male_claw_cor_extreme)/mean_male_claw_cor_extreme,
               prop_eaten_male2 = (male_rate-male_rate)/male_rate,
               prop_eaten_female2 = (female_rate - male_rate)/male_rate,
               #this first delta metric compares both the iteration-specific 
               #male and female rates to the mean male-only population
               #(i.e., it incorporates uncertainty around both the male and 
               #female estimates)
               delta = (1 - prop_female) * prop_eaten_male + 
                 prop_female * prop_eaten_female,
               #this second delta metric only incorporates uncertainty around 
               #the female estimate
               delta2 = (1 - prop_female) * prop_eaten_male2 + 
                 prop_female * prop_eaten_female2,
               iteration = i)
  ci_df_claw_cor_extreme <- bind_rows(ci_df_claw_cor_extreme, df)
}

summ_cis_claw_cor_extreme <- ci_df_claw_cor_extreme %>% 
  group_by(prop_female) %>% 
  summarize(upper = quantile(delta, 0.975),
            lower = quantile(delta, 0.025),
            median = quantile(delta, 0.5),
            upper2 = quantile(delta2, 0.975),
            lower2 = quantile(delta2, 0.025),
            median2 = quantile(delta2, 0.5))

#now plot the real populations on the hypothetical figure
real_props_claw_cor_extreme <- tibble(prop_female = m_to_f_first$prop_female,
                              unique_id = m_to_f_first$unique_id,
                              site = m_to_f_first$site,
                              #set to -Inf so they'll be right along the axis
                              delta = -Inf)


ggplot(summ_cis_claw_cor_extreme, aes(prop_female)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = median)) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Predicted difference in the proportion of clams\neaten relative to a male-only population",
       x = "Proportion of female crabs in population") +
  scale_shape_identity() +
  geom_point(data = real_props_claw_cor_extreme, aes(y = delta), 
             shape = 108, size = 5) +
  theme(legend.position = "none")
#and they look virtually identical
#ggsave("claw_corrected_predictions_extreme.png", height = 5, width = 7)


#SEX RATIOS AND REPEATED SAMPLING----------------------------------------------
mod_repeat <- glmer(prop_female ~ site + (1|site:month), 
                    weights = total,
                    family = binomial,
                    data = m_to_f)
summary(mod_repeat)

repeat_predict <- ggeffects::ggpredict(mod_repeat, "site") %>% 
  rename(site = x)

ggplot(m_to_f, aes(x = site, y = prop_female)) +
  #plot the raw data
  geom_jitter(width = 0.2, alpha = 0.3) +
  #now plot the  model output
  geom_point(data = repeat_predict,
             aes(y = predicted), size = 3) + 
  geom_errorbar(data = repeat_predict,
              aes(y = predicted, ymin = conf.low, ymax = conf.high), 
              width = 0.2) +
  theme_paper_large() +
  labs(x = "Site",
       y = "Proportion of females in sample")
#ggsave("population_proportions.pdf",width = 9, height = 6)

emmean_ratio <- emmeans::emmeans(mod_repeat, specs = pairwise ~ "site", 
                                 type = "response")
emmean_ratio
plot(emmean_ratio)

m_to_f %>% 
  unite(group, c(site, month, year)) %>% 
  ggplot(aes(day, prop_female)) +
    geom_point() +
    facet_wrap(~group) +
    stat_smooth(method = 'lm')

#sex ratios accumulated
mod_accum <- glmer(prop_female_accum ~ site + (1|site:month), 
                    weights = total_accum,
                    family = binomial,
                    data = m_to_f_accum)
summary(mod_accum)

accum_predict <- ggeffects::ggpredict(mod_accum, "site") %>% 
  rename(site = x)

ggplot(m_to_f_accum, aes(x = site, y = prop_female_accum)) +
  #plot the raw data
  geom_jitter(width = 0.2, alpha = 0.3) +
  #now plot the  model output
  geom_point(data = accum_predict,
             aes(y = predicted), size = 3) + 
  geom_errorbar(data = repeat_predict,
                aes(y = predicted, ymin = conf.low, ymax = conf.high), 
                width = 0.2) +
  theme_paper_large() +
  labs(x = "Site",
       y = "Proportion of females in sample")
#ggsave("population_proportions_accumulated_catch.pdf",width = 9, height = 6)
emmean_ratio_accum <- emmeans::emmeans(mod_accum, specs = pairwise ~ "site", 
                                       type = "response")
emmean_ratio_accum
plot(emmean_ratio_accum)
