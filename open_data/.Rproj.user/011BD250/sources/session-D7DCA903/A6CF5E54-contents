############### PRELIMINARIES ###############
rm(list=ls())

# installing and loading packages
packages <- c("readxl", "ggsignif", "tidyverse", "ggplot2", "nlme", "lme4",
              "sjPlot", "psych", "lsmeans", "ggpubr","methods", "lmerTest", "gtsummary", "flextable")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

# load dataset - make sure working directory is parent folder
load("open_dataset.RData")

# redefine variables for later use
conditioning_concat <- conditioning_concat_export
lucky_concat <- lucky_concat_export
gen_concat <- gen_concat_export
mem_concat <- mem_concat_export
survey_dat <- survey_dat_export %>% mutate(education = as.numeric(education))

# define luckiness data (conditioning)
cs_img_info <- conditioning_concat %>% 
  select(c("drink_2grp", "subId", "stim", "cs_img","expt")) %>% distinct()
cs_img_info <- cs_img_info %>%
  mutate(cs_img_name = ifelse(startsWith(str_sub(cs_img_info$cs_img,-11), "V") ,str_sub(cs_img_info$cs_img,-11),str_sub(cs_img_info$cs_img,-10)))

lucky_concat <- lucky_concat %>%
  rename(cs_img_name = stim)

lucky_compare <- merge(lucky_concat, cs_img_info, by = c("subId", "cs_img_name","drink_2grp", "expt")) %>%
  mutate(expt = factor(expt, levels = c("reward","loss")),
         response = as.numeric(response))

# define RT data (conditioning)
us_rt_dat <- conditioning_concat %>%
  dplyr::filter(us_acc == 1) %>%
  mutate(cs_rt_log = log(cs_rt),
         us_rt_log = log(us_rt)) %>%
  mutate(stim = ifelse(grepl("surprise", trial_cat), "null",
                       ifelse(grepl("alc", trial_cat), "alc", 
                              ifelse(grepl("neut", trial_cat), "neut", 
                                     ifelse(grepl("null", trial_cat), "null", NA)))),
         timing = ifelse(iter %in% c(1:5), 1,
                         ifelse(iter %in% c(6:10), 2,
                                ifelse(iter %in% c(11:15), 3,
                                       ifelse(iter %in% c(16:20), 4,
                                              ifelse(iter %in% c(21:25), 5,
                                                     ifelse(iter %in% c(26:30), 6,
                                                            ifelse(iter %in% c(31:35), 7,
                                                                   ifelse(iter %in% c(36:40), 8, NA)))))))),
         earlyLate = ifelse(trial_num < 59, "early", "late"),
         cong = ifelse(grepl("incong|surprise", trial_cat), "Unexpected", "Expected")) %>%
  filter(trial_cat != "null")

# define choice data (generalization)
gen_final <- gen_concat %>% 
  mutate(stim = ifelse(grepl("alc", trial_cat), "alc",
                       ifelse(grepl("neut", trial_cat), "neut",
                              ifelse(grepl("null", trial_cat), "null",
                                     ifelse(grepl("novel", trial_cat), "novel", NA))))) %>%
  merge(cs_img_info %>% select(-c("drink_2grp")),
        by = c("subId", "stim", "expt"), all.x = TRUE) %>%
  mutate(play_bin = ifelse(gs_choice == -1, 0, 
                           ifelse(gs_choice == 1, 1, NA)),
         deg_plot = ifelse(gs_deg == 999, 0, gs_deg),
         deg_scale = sjmisc::center(as.numeric(deg_plot)),
         stim = factor(stim, levels = c('null', 'novel', 'neut', 'alc')),
         drink_2grp = factor(drink_2grp, levels = c('Light', 'Risky')))

# define RT data (generalization)
probe_rt_dat <- gen_concat %>%
  filter((trial_cat != "gen_novel") & (probe_acc == 1)) %>%
  mutate(gs_rt_log = log(gs_rt),
         probe_rt_log = log(probe_rt)) %>%
  mutate(deg_plot = ifelse(gs_deg == 999, 0, gs_deg),
         trial_cat = factor(trial_cat, levels = c("gen_neut","gen_null", "gen_alc"))) %>%
  filter_all(any_vars(is.na(.)))

# define memory data
mem_agg <- mem_concat %>%
  filter(!(trial_cat == "prac")) %>%
  mutate(acc = ifelse(img_type=="old" & choice %in% c("Confident Old", "Unsure Old"), "hit",
                      ifelse(img_type=="new" & choice %in% c("Confident Old", "Unsure Old"), "fa",
                             ifelse(img_type=="old" & choice %in% c("Confident New", "Unsure New"), "miss",
                                    ifelse(img_type=="new" & choice %in% c("Confident New", "Unsure New"), "cr", NA))))) %>%
  group_by(drink_2grp, subId, trial_cat, img_type, acc, expt) %>%
  dplyr::summarise(num_resp = n()) %>%
  merge(mem_concat %>%
          group_by(subId, audit_tot, trial_cat, img_type) %>%
          dplyr::summarise(tot_img = n()), by = c("subId", "trial_cat", "img_type")) %>%
  mutate(rate = num_resp/tot_img) %>%
  select(-c(img_type,num_resp)) %>% 
  filter(!(is.na(acc))) %>% 
  spread(acc, rate) %>% 
  group_by(subId, trial_cat) %>% 
  dplyr::mutate(hit = sum(hit, na.rm = T), 
                fa = sum(fa, na.rm = T), 
                miss = sum(miss, na.rm = T), 
                cr = sum(cr, na.rm = T)) %>%
  mutate(dummy.fa = ifelse(fa==0, .0001, ifelse(fa==1, .999, fa)),
         dummy.hit = ifelse(hit==0, .0001, ifelse(hit==1, .999, hit)),
         dprime = qnorm(dummy.hit) - qnorm(dummy.fa),
         criterion = -1*rowMeans(cbind(qnorm(dummy.hit/40), qnorm(dummy.fa/40)))) %>%
  mutate(expt = ifelse(expt=="reward", "Reward", "Loss"),
         expt = factor(expt, levels = c("Reward", "Loss")))

# function to mean-center audit within groups
prep_meancent_audit <- function (df, grpvar, grpname, auditvar) {
  dat_subset <- df[which(df[,grpvar]==grpname),]
  dat_subset$audit_meancent <- unlist(dat_subset[,auditvar]) - mean(unlist(dat_subset[,auditvar]), na.rm = T)
  return(dat_subset)
}

############### DEMOGRAPHICS ANALYSIS ###############
table(survey_dat$drink_2grp, survey_dat$expt)

# study-wide table
demog_table <- survey_dat %>%
  select(!c(subId, binger)) %>%
  tbl_strata(
    strata = expt, .tbl_fun = ~.x %>%
      tbl_summary(by = drink_2grp,
                  statistic = list(all_continuous() ~ "{mean} ({sd})"),
                  type = list(education ~ "continuous"),
                  label = list(gender ~ "Sex", age ~ "Age", ethnicity ~ "Ethnicity", race ~ "Race", education ~ "Education Level", audit_tot ~ "AUDIT total score", pss_tot ~ "PSS total score")) %>%
      add_n() %>%
      add_p(test = list(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test")) %>%
      separate_p_footnotes() %>%
      bold_labels())

demog_table

# 1. study 1 tests
survey_dat_reward <- survey_dat %>% filter(expt == "reward")

chisq.test(table(survey_dat_reward$drink_2grp, survey_dat_reward$gender))
chisq.test(table(survey_dat_reward$drink_2grp, survey_dat_reward$ethnicity))
chisq.test(table(survey_dat_reward$drink_2grp, survey_dat_reward$race))

t.test(age ~ drink_2grp, data=survey_dat_reward, var.equal=TRUE)
t.test(education ~ drink_2grp, data=survey_dat_reward, var.equal=TRUE)
t.test(pss_tot ~ drink_2grp, data=survey_dat_reward, var.equal=TRUE)
t.test(audit_tot ~ drink_2grp, data=survey_dat_reward, var.equal=TRUE)
t.test(drinks_per_week ~ drink_2grp, data=survey_dat_reward, var.equal=TRUE)

#2. study 2 tests
survey_dat_loss <- survey_dat %>% filter(expt == "loss")

chisq.test(table(survey_dat_loss$drink_2grp, survey_dat_loss$gender))
chisq.test(table(survey_dat_loss$drink_2grp, survey_dat_loss$ethnicity))
chisq.test(table(survey_dat_loss$drink_2grp, survey_dat_loss$race))

t.test(age ~ drink_2grp, data=survey_dat_loss, var.equal=TRUE)
t.test(education ~ drink_2grp, data=survey_dat_loss, var.equal=TRUE)
t.test(pss_tot ~ drink_2grp, data=survey_dat_loss, var.equal=TRUE)
t.test(audit_tot ~ drink_2grp, data=survey_dat_loss, var.equal=TRUE)
t.test(drinks_per_week ~ drink_2grp, data=survey_dat_loss, var.equal=TRUE)

############### CONDITIONING ANALYSIS ###############

# 1-1. study 1 luckiness rating
lucky_compare_reward <- lucky_compare %>% filter(expt == "reward")
img_luck_mod_reward <- lme(response ~ timing*stim*drink_2grp, random=~1|subId, 
                           dat = lucky_compare_reward %>%
                             mutate(timing = factor(timing, levels = c("post_cond", "baseline")),
                                    stim = factor(stim, levels = c("null", "alc", "neut"))), na.action = na.omit)
car::Anova(img_luck_mod_reward, type=2, test.statistic = "Wald")

# 1-2. posthoc
emmeans(img_luck_mod_reward, pairwise ~ timing|stim)$contrasts
emmeans(img_luck_mod_reward, pairwise ~ timing|stim|drink_2grp)$contrasts

# 1-3. drinking severity effect
cond_reward_light <- prep_meancent_audit(lucky_compare_reward, "drink_2grp", "Light", "audit_tot")
cond_reward_light_mod <- lme(response ~ timing*stim*audit_meancent, random=~1|subId, dat = cond_reward_light)
car::Anova(cond_reward_light_mod, test = 2, test.statistic = "Wald")

cond_reward_risky <- prep_meancent_audit(lucky_compare_reward, "drink_2grp", "Risky", "audit_tot")
cond_reward_risky_mod <- lme(response ~ timing*stim*audit_meancent, random=~1|subId, dat = cond_reward_risky)
car::Anova(cond_reward_risky_mod, test = 2, test.statistic = "Wald")

# 1-4. attentional orienting
us_rt_dat_reward <- us_rt_dat %>% filter(expt == "reward")
implic_cond_iter_reward <- lme(us_rt_log ~ iter*cong*drink_2grp, random=~1|subId, 
                               dat = us_rt_dat_reward %>%
                                 group_by(drink_2grp, subId, iter, cong))
car::Anova(implic_cond_iter_reward, type=2, statistical.test="Wald")

# 1-5. drinking severity effect
cond_reward_light <- prep_meancent_audit(us_rt_dat_reward, "drink_2grp", "Light", "audit_tot")
cond_reward_light_mod <- lme(us_rt_log ~ iter*cong*audit_meancent, random=~1|subId, dat = cond_reward_light)
car::Anova(cond_reward_light_mod, test = 2, test.statistic = "Wald")

cond_reward_risky <- prep_meancent_audit(us_rt_dat_reward, "drink_2grp", "Risky", "audit_tot")
cond_reward_risky_mod <- lme(us_rt_log ~ iter*cong*audit_meancent, random=~1|subId, dat = cond_reward_risky)
car::Anova(cond_reward_risky_mod, test = 2, test.statistic = "Wald")

# 2-1. study 2 luckiness rating
lucky_compare_loss <- lucky_compare %>% filter(expt == "loss")
img_luck_mod_loss <- lme(response ~ timing*stim*drink_2grp, random=~1|subId, 
                         dat = lucky_compare_loss, na.action = na.omit)
car::Anova(img_luck_mod_loss, type=2, test="Wald")

# 2-2. posthoc
emmeans(img_luck_mod_loss, pairwise ~ timing|stim)
emmeans(img_luck_mod_loss, pairwise ~ timing|stim|drink_2grp)

# 2-3. drinking severity effect
cond_loss_light <- prep_meancent_audit(lucky_compare_loss, "drink_2grp", "Light", "audit_tot")
cond_loss_light_mod <- lme(response ~ timing*stim*audit_meancent, random=~1|subId, dat = cond_loss_light)
car::Anova(cond_loss_light_mod, test = 2, test.statistic = "Wald")

cond_loss_risky <- prep_meancent_audit(lucky_compare_loss, "drink_2grp", "Risky", "audit_tot")
cond_loss_risky_mod <- lme(response ~ timing*stim*audit_meancent, random=~1|subId, dat = cond_loss_risky)
car::Anova(cond_loss_risky_mod, test = 2, test.statistic = "Wald")

# 2-4. attentional orienting
us_rt_dat_loss <- us_rt_dat %>% filter(expt == "loss")
implic_cond_iter_loss <- lme(us_rt_log ~ iter*cong*drink_2grp, random=~1|subId, 
                             dat = us_rt_dat_loss %>% 
                               group_by(drink_2grp, subId, iter, cong))
car::Anova(implic_cond_iter_loss, test=2, statistical.test="Wald")

# 2-5. drinking severity effect
cond_loss_light <- prep_meancent_audit(us_rt_dat_loss, "drink_2grp", "Light", "audit_tot")
cond_loss_light_mod <- lme(us_rt_log ~ iter*cong*audit_meancent, random=~1|subId, dat = cond_loss_light)
car::Anova(cond_loss_light_mod, test = 2, test.statistic = "Wald")

cond_loss_risky <- prep_meancent_audit(us_rt_dat_reward, "drink_2grp", "Risky", "audit_tot")
cond_loss_risky_mod <- lme(us_rt_log ~ iter*cong*audit_meancent, random=~1|subId, dat = cond_loss_risky)
car::Anova(cond_loss_risky_mod, test = 2, test.statistic = "Wald")

############### GENERALIZATION ANALYSIS ###############

# 1-1. study 1 choice generalization
gen_final_reward <- gen_final %>% filter(expt == "reward")
binmod_2grp_reward <- glm(play_bin ~ trial_cat*deg_plot*drink_2grp + (trial_num+1), dat = gen_final_reward %>% filter(!(trial_cat == "gen_novel")) %>%
                            mutate(play_bin = ifelse(gs_choice == -1, 0, 
                                                     ifelse(gs_choice == 1, 1, NA)),
                                   deg_plot = ifelse(gs_deg == 999, 0, gs_deg), 
                                   trial_cat = factor(trial_cat, levels = c("gen_neut", "gen_null", "gen_alc")),
                                   drink_2grp = factor(drink_2grp, levels = c("Light", "Risky"))) %>%
                            mutate(trial_cat = ifelse(trial_cat == "gen_neut", "CS+obj",
                                                      ifelse(trial_cat == "gen_alc", "CS+alc",
                                                             ifelse(trial_cat == "gen_null", "CS-", NA))),
                                   trial_cat = factor(trial_cat, levels = c("CS+obj", "CS+alc", "CS-"))),
                          family = "binomial")

summary(binmod_2grp_reward)
car::Anova(binmod_2grp_reward, test.statistic = "Wald")

# 1-2. posthoc
emmeans(binmod_2grp_reward, pairwise ~ trial_cat)$contrasts
emmeans(binmod_2grp_reward, pairwise ~ drink_2grp)$contrasts
emmeans(binmod_2grp_reward, pairwise ~ drink_2grp|trial_cat)$contrasts

emtrends(binmod_2grp_reward, pairwise ~ trial_cat, var = "deg_plot")
emtrends(binmod_2grp_reward, pairwise ~ drink_2grp|trial_cat, var = "deg_plot")

# 1-3. drinking severity effects
gen_concat.light <- as.data.frame(prep_meancent_audit(gen_final_reward, "drink_2grp", "Light", "audit_tot"))

explicit_gen.light <- glm(play_bin ~ deg_plot*trial_cat*audit_meancent + (trial_num+1), 
                          dat = gen_concat.light %>% filter(!(trial_cat == "gen_novel")) %>%
                            mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null"))),
                          family = "binomial")

car::Anova(explicit_gen.light, test = 2, test.statistic = "Wald")
summary(explicit_gen.light)

gen_concat.risky <- prep_meancent_audit(gen_final_reward, "drink_2grp", "Risky", "audit_tot")

explicit_gen.risky <- glm(play_bin ~ deg_plot*trial_cat*audit_meancent + (trial_num+1), 
                          dat = gen_concat.risky %>% filter(!(trial_cat == "gen_novel")) %>%
                            mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null"))),
                          family = "binomial")

car::Anova(explicit_gen.risky, test = 2, test.statistic = "Wald")
summary(explicit_gen.risky)

# 1-4 exploratory: median split
explicit_gen.light_medsplit <- glm(play_bin ~ deg_plot*trial_cat*audit_medsplit + (trial_num+1), 
                                   dat = gen_final_reward %>% filter(drink_2grp == "Light") %>% filter(!(trial_cat == "gen_novel")) %>%
                                     mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null")), audit_medsplit = ifelse(audit_tot<median(audit_tot), "Low", "High")),
                                   family = "binomial")

emtrends(explicit_gen.light_medsplit, pairwise ~ audit_medsplit|trial_cat, var = "deg_plot")

explicit_gen.risky_medsplit <- glm(play_bin ~ deg_plot*trial_cat*audit_medsplit + (trial_num+1), 
                                   dat = gen_final_reward %>% filter(drink_2grp == "Risky") %>% filter(!(trial_cat == "gen_novel")) %>%
                                     mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null")), audit_medsplit = ifelse(audit_tot<median(audit_tot), "Low", "High"),
                                            audit_medsplit = factor(audit_medsplit, levels = c("Low", "High"))),
                                   family = "binomial")

emtrends(explicit_gen.risky_medsplit, pairwise ~ audit_medsplit|trial_cat, var = "deg_plot")

# 1-5. attentional orienting
rt_probe_gradient <- lme(probe_rt_log ~ deg_plot*drink_2grp*probe_type, random = ~1|subId, dat = probe_rt_dat %>%
                           filter(expt == "reward"))
car::Anova(rt_probe_gradient, type = 2, test.statistic = "Wald")

# 2-1. study 2 choice generalization
gen_final_loss <- gen_final %>% filter(expt == "loss")
binmod_2grp_loss <- glm(play_bin ~ trial_cat*deg_plot*drink_2grp + (trial_num+1), dat = gen_final_loss %>% filter(!(trial_cat == "gen_novel")) %>%
                          mutate(play_bin = ifelse(gs_choice == -1, 0, 
                                                   ifelse(gs_choice == 1, 1, NA)),
                                 deg_plot = ifelse(gs_deg == 999, 0, gs_deg), 
                                 trial_cat = factor(trial_cat, levels = c("gen_null", "gen_alc", "gen_neut")),
                                 drink_2grp = factor(drink_2grp, levels = c("Light", "Risky"))),
                        family = "binomial")
car::Anova(binmod_2grp_loss, test = 2, test.statistic = "Wald")
summary(binmod_2grp_loss)

# 2-2. posthoc
emmeans(binmod_2grp_loss, pairwise ~ trial_cat)$contrasts
emmeans(binmod_2grp_loss, pairwise ~ drink_2grp|trial_cat)$contrasts

emtrends(binmod_2grp_loss, pairwise ~ trial_cat, var = "deg_plot")$contrasts
emtrends(binmod_2grp_loss, pairwise ~ drink_2grp|trial_cat, var = "deg_plot")$contrasts

# 2-3. drinking severity effects
gen_concat.light <- as.data.frame(prep_meancent_audit(gen_final_loss, "drink_2grp", "Light", "audit_tot"))
explicit_gen.light <- glm(play_bin ~ deg_plot*trial_cat*audit_meancent + (trial_num+1), 
                          dat = gen_concat.light %>% filter(!(trial_cat == "gen_novel")) %>%
                            mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null"))),
                          family = "binomial")
car::Anova(explicit_gen.light, test = 2, test.statistic = "Wald")

gen_concat.risky <- prep_meancent_audit(gen_final_loss, "drink_2grp", "Risky", "audit_tot")
explicit_gen.risky <- glm(play_bin ~ deg_plot*trial_cat*audit_meancent + (trial_num+1), 
                          dat = gen_concat.risky %>% filter(!(trial_cat == "gen_novel")) %>%
                            mutate(trial_cat = factor(trial_cat, levels = c( "gen_null", "gen_alc", "gen_neut"))),
                          family = "binomial")
car::Anova(explicit_gen.risky, test = 2, test.statistic = "Wald")
summary(explicit_gen.risky)

# 2-4. exploratory: median split
explicit_gen.light_medsplit <- glm(play_bin ~ deg_plot*trial_cat*audit_medsplit + (trial_num+1), 
                                   dat = gen_final_loss %>% filter(drink_2grp == "Light") %>% filter(!(trial_cat == "gen_novel")) %>%
                                     mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null")), audit_medsplit = ifelse(audit_tot<median(audit_tot), "Low", "High")),
                                   family = "binomial")

emtrends(explicit_gen.light_medsplit, pairwise ~ audit_medsplit|trial_cat, var = "deg_plot")

explicit_gen.risky_medsplit <- glm(play_bin ~ deg_plot*trial_cat*audit_medsplit + (trial_num+1), 
                                   dat = gen_final_loss %>% filter(drink_2grp == "Risky") %>% filter(!(trial_cat == "gen_novel")) %>%
                                     mutate(trial_cat = factor(trial_cat, levels = c( "gen_neut", "gen_alc", "gen_null")), audit_medsplit = ifelse(audit_tot<median(audit_tot), "Low", "High")),
                                   family = "binomial")

emtrends(explicit_gen.risky_medsplit, pairwise ~ trial_cat|audit_medsplit, var = "deg_plot")

# 2-5. attentional orienting
rt_probe_gradient <- lme(probe_rt_log ~ deg_plot*drink_2grp*probe_type, random = ~1|subId, dat = probe_rt_dat %>%
                           filter(expt == "loss"))
car::Anova(rt_probe_gradient, type = 2, test.statistic = "Wald")
summary(rt_probe_gradient)


############### MEMORY ANALYSIS ###############

# 1-1. study 1 dprime
mem_agg_reward <- mem_agg %>% filter(expt == "Reward")
dprime_mod_reward <- lme(dprime ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_reward)
car::Anova(dprime_mod_reward, test.statistic = "Wald")

emmeans(dprime_mod_reward, pairwise ~ trial_cat)$contrasts

# 1-2. drinking severity effect
mem_reward.light <- prep_meancent_audit(mem_agg_reward, "drink_2grp", "Light", "audit_tot")
dprime_reward.light <- lme(dprime ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_reward.light)
car::Anova(dprime_reward.light, test = 2, test.statistic = "Wald")

mem_reward.risky <- prep_meancent_audit(mem_agg_reward, "drink_2grp", "Risky", "audit_tot")
dprime_reward.risky <- lme(dprime ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_reward.risky)
car::Anova(dprime_reward.risky, test = 2, test.statistic = "Wald")

# 1-3. study 1 hit
hit_mod_reward <- lme(hit ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_reward)
car::Anova(hit_mod_reward, type=2, test.statistic="Wald")

emmeans(hit_mod_reward, pairwise ~ drink_2grp|trial_cat)$contrasts
emmeans(hit_mod_reward, pairwise ~ trial_cat|drink_2grp)$contrasts

# 1-4. drinking severity effect
mem_concat.light <- as.data.frame(prep_meancent_audit(mem_agg_reward, "drink_2grp", "Light", "audit_tot"))
hit_mod_reward <- lme(hit ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_concat.light)
car::Anova(hit_mod_reward, type=2, test.statistic="Wald")

mem_concat.risky <- as.data.frame(prep_meancent_audit(mem_agg_reward, "drink_2grp", "Risky", "audit_tot"))
hit_mod_reward <- lme(hit ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_concat.risky)
summary(hit_mod_reward)
car::Anova(hit_mod_reward, type=2, test.statistic="Wald")

# 1-5. study 1 fa
fa_mod_reward <- lme(fa ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_reward)
car::Anova(fa_mod_reward)

emmeans(fa_mod_reward, pairwise ~ drink_2grp|trial_cat)$contrasts
emmeans(fa_mod_reward, pairwise ~ trial_cat|drink_2grp)$contrasts

# 1-6. drinking severity effect
mem_concat.light <- as.data.frame(prep_meancent_audit(mem_agg_reward, "drink_2grp", "Light", "audit_tot"))
fa_mod_reward <- lme(fa ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_concat.light)
car::Anova(fa_mod_reward, type=2, test.statistic="Wald")
summary(fa_mod_reward)

mem_concat.risky <- as.data.frame(prep_meancent_audit(mem_agg_reward, "drink_2grp", "Risky", "audit_tot"))
fa_mod_reward <- lme(fa ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_concat.risky)
summary(fa_mod_reward)
car::Anova(fa_mod_reward, type=2, test.statistic="Wald")

# 2-1. study 2 dprime
mem_agg_loss <- mem_agg %>% filter(expt == "Loss")
dprime_mod_loss <- lme(dprime ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_loss)
car::Anova(dprime_mod_loss)

emmeans(dprime_mod_loss, pairwise ~ trial_cat)$contrasts

# 2-2. drinking severity effect
mem_loss.light <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Light", "audit_tot")
dprime_loss.light <- lme(dprime ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.light)
car::Anova(dprime_loss.light, test = 2, test.statistic = "Wald")

mem_loss.risky <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Risky", "audit_tot")
dprime_loss.risky <- lme(dprime ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.risky)
car::Anova(dprime_loss.risky, test = 2, test.statistic = "Wald")

# 2-3. study 2 hit
hit_mod_loss <- lme(hit ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_loss)
car::Anova(hit_mod_loss, type = 2, test.statistic = "Wald")

# 2-4. drinking severity effect
mem_loss.light <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Light", "audit_tot")
hit_loss.light <- lme(hit ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.light)
car::Anova(hit_loss.light, test = 2, test.statistic = "Wald")

mem_loss.risky <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Risky", "audit_tot")
hit_loss.risky <- lme(hit ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.risky)
car::Anova(hit_loss.risky, test = 2, test.statistic = "Wald")

# 2-5. study 2 fa
fa_mod_loss <- lme(fa ~ trial_cat*drink_2grp, random=~1|subId, dat = mem_agg_loss)
car::Anova(fa_mod_loss, type = 2, test.statistic = "Wald")

# 2-6. drinking severity effect
mem_loss.light <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Light", "audit_tot")
fa_loss.light <- lme(fa ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.light)
car::Anova(fa_loss.light, test = 2, test.statistic = "Wald")

mem_loss.risky <- prep_meancent_audit(mem_agg_loss, "drink_2grp", "Risky", "audit_tot")
fa_loss.risky <- lme(fa ~ trial_cat*audit_meancent, random=~1|subId, dat = mem_loss.risky)
car::Anova(fa_loss.risky, test = 2, test.statistic = "Wald")
