

#Subsetting the data
t_subset <- seq_table_filtered %>% 
  filter(rel_abundance > 0) %>% 
  filter(experiment %in% c("CCME11", "CCME16")) %>% 
  filter(day == 28) %>% 
  filter(group %in% c("SPF", "PBS", "BpKH6", "BpSCSK")) %>% 
  mutate(group = factor(group, levels = c("SPF", "PBS", "BpKH6", "BpSCSK"))) %>% 
  mutate(day = factor(day, levels = unique(day))) %>% 
  group_by(group, day) %>% 
  mutate(mouse.total = n_distinct(mouse.number)) %>% 
  ungroup()


##Table 1
t1 <- t_subset %>%
  select(Phylum, Class, Order, Family, asv, group) %>% 
  group_by( Phylum, Class, Order, Family, group) %>% 
  mutate(count = n_distinct(asv)) %>% 
  ungroup() %>% 
  select(-asv) %>% 
  distinct() %>% 
  spread(group, count) %>% 
  mutate_all(list(~ replace(., is.na(.), "-"))) %>% 
  mutate(Phylum = if_else(duplicated(Phylum), NA, Phylum),
         Class = if_else(duplicated(Class), NA, Class),
         Order = if_else(duplicated(Order), NA, Order))

write_xlsx(t1, "./tables/table_1.xlsx")

rm(t1, t_subset)

