library(synapseClient)
synapseLogin()
library(dplyr)
library(reshape2)

library(stringr)
library(tidyr)

FactorCovariates <- c('C4_Cell_Line_ID',
                      'Diffname_short', 'run', 'lane', 'PassageAtThaw', 'PassageAtHarvest',
                      'Cell_Line_Type', 'Cell_Line_of_Origin', 'Tissue_of_Origin', 
                      'Reprogramming_Gene_Combination', 'Culture_Conditions', 'Donor_Life_Stage',
                      'Race', 'Ethnicity' , 'Gender', 'Disease', 'Originating_Lab', 'Donor_ID',
                      'Cell_Type_of_Origin_Level2', 'Reprogramming_Vector_Type', 'Reprogramming_Vector_Type_Level2',
                      'XIST_Level')

# For diffname_short values
diffstates <- c("DE", "EB", "ECTO", "MESO-5", "MESO-15", "MESO-30")

# Taking this from the 
res <- synTableQuery(sprintf("select * from syn3156503"))@values %>%
  filter(Diffname_short %in% diffstates) %>%
  select(one_of(FactorCovariates))

res2 <- melt(res, id.vars = "C4_Cell_Line_ID") %>% 
  select(-C4_Cell_Line_ID) %>%
  unique() %>%
  filter(!(value %in% c("N/A", "")), !is.na(value))

res3 <- res2 %>% 
  left_join(res2, by=c("variable")) %>%
  filter(as.character(value.x) < as.character(value.y)) %>%
  rename(class=variable, variable1=value.x, variable2=value.y) %>%
  mutate(variable1short=str_replace_all(variable1, '[^[:alnum:]]', ''),
         variable2short=str_replace_all(variable2, '[^[:alnum:]]', ''))# %>%
  # select(class, variable1, variable1short, variable2, variable2short)

res3diff <- res3 %>% filter(class == "Diffname_short") %>% mutate(datarestriction="All")
res3rest <- res3 %>% filter(class != "Diffname_short")

origlen <- nrow(res3rest)
restdiffstates <- c("All", diffstates)
res3rest <- res3rest[rep(1:origlen, each=length(restdiffstates)), ] %>%
  mutate(datarestriction=rep(restdiffstates, times=origlen))

res4 <- rbind(res3diff, res3rest) %>%
  unite(tmp, variable1short, variable2short, sep="_vs_", remove=FALSE) %>%
  unite(shortname, datarestriction, tmp, sep="__", remove=FALSE) %>%
  select(class, datarestriction, variable1, variable1short, variable2, variable2short, shortname, -tmp)
