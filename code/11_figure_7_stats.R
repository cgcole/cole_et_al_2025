

c_diff_symptom_t_subset <- c_diff_symptom_t %>% 
  filter(mouse.number != 175) %>% 
  group_by(mouse.number, group) %>% 
  summarise(weight_loss_percentage = 100 - min(weight_percentage)) %>% 
  ungroup() %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK")))


c_diff_cfu_t_subset <- c_diff_cfu_t %>% 
  filter(mouse.number != 175) %>% 
  group_by(mouse.number, experiment.sample.name, day, days.post.infection,
           sampleid, group, treatment, Bacteria) %>% 
  summarise(CFU.g = mean(CFU.g),
            n = n()) %>% 
  ungroup() %>% 
  mutate(CFU.g = ifelse(CFU.g <= 1000, 1000, as.numeric(CFU.g))) %>% 
  group_by(day, group) %>% 
  mutate(mean_CFU.g = mean(CFU.g)) %>% 
  ungroup() %>% 
  mutate(group = str_replace(group, "SPF", "SPF + NVMC")) %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + NVMC", "PBS", "BpSCSK")))


kleb_cfu_t_subset <- kleb_cfu_t %>% 
  mutate(group = str_replace(group, "SPF", "SPF + Amp")) %>% 
  group_by(mouse.number, experiment.sample.name, day, days.post.infection,
           sampleid, group, treatment, Bacteria) %>% 
  summarise(CFU.g = mean(CFU.g),
            n = n()) %>% 
  ungroup() %>% 
  mutate(CFU.g = ifelse(CFU.g <= 1000, 1000, as.numeric(CFU.g))) %>% 
  group_by(day, group) %>% 
  mutate(mean_CFU.g = mean(CFU.g)) %>% 
  ungroup() %>% 
  filter(group != "BpKH6") %>% 
  mutate(group = factor(group, levels = c("SPF + Amp", "PBS", "BpSCSK"))) 



##Kleb CFU stats

#day1
kleb_cfu_day1 <- kleb_cfu_t_subset %>% 
  filter(days.post.infection == 1)

lm_model <- lm(CFU.g ~ group, kleb_cfu_day1)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Was not normal

levene_test(kleb_cfu_day1, CFU.g ~ group, center = median)

#not normal and unequal variance
kruskal_test(kleb_cfu_day1, CFU.g ~ group)
kleb_day1_dunn <- dunn_test(kleb_cfu_day1, CFU.g ~ group, p.adjust.method = "none")

kleb_day1_p_val <- c(0.001670212, 0.116024973, 0.116024973)


#day3
kleb_cfu_day3 <- kleb_cfu_t_subset %>% 
  filter(days.post.infection == 3)


lm_model <- lm(CFU.g ~ group, kleb_cfu_day3)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") 

bartlett.test(CFU.g ~ group, kleb_cfu_day3)


#Normal but unequal variance
oneway.test(CFU.g ~ group, data = kleb_cfu_day3, var.equal = FALSE)
kleb_day3_welch <- pairwise.t.test(kleb_cfu_day3$CFU.g, kleb_cfu_day3$group, p.adjust.method = "none", pool.sd = FALSE)

kleb_day3_p_val <- c(0.1507982, 0.1613437, 0.02747215)


#day7
kleb_cfu_day7 <- kleb_cfu_t_subset %>% 
  filter(days.post.infection == 7)

lm_model <- lm(CFU.g ~ group, kleb_cfu_day7)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") 

levene_test(kleb_cfu_day7, CFU.g ~ group, center = median)

#Not normal but equal variance
kruskal_test(kleb_cfu_day7, CFU.g ~ group)
kleb_day7_dunn <- dunn_test(kleb_cfu_day7, CFU.g ~ group, p.adjust.method = "none")

kleb_day7_p_val <- c(0.429522788, 0.048265728, 0.005691759)


p.adjust(c(kleb_day1_p_val, kleb_day3_p_val, kleb_day7_p_val), method = "BH")




##C. difficile CFU stats

#day7
c_diff_cfu_day7 <- c_diff_cfu_t_subset %>% 
  filter(days.post.infection == 7)

lm_model <- lm(CFU.g ~ group, c_diff_cfu_day7)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") 

bartlett.test(CFU.g ~ group, c_diff_cfu_day7)

#Normal and Equal variance
aov_model <- aov(c_diff_cfu_day7$CFU.g ~ c_diff_cfu_day7$group)
summary(aov_model)
TukeyHSD(aov_model)



##C. difficile Weight loss
lm_model <- lm(weight_loss_percentage ~ group, c_diff_symptom_t_subset)
residuals <- residuals(lm_model)
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col = "red") #Was not normal

levene_test(c_diff_symptom_t_subset, weight_loss_percentage ~ group, center = median)

#Not normal but equal variance
kruskal_test(c_diff_symptom_t_subset, weight_loss_percentage ~ group)
dunn_test(c_diff_symptom_t_subset, weight_loss_percentage ~ group, p.adjust.method = "BH")







