# Create all pairwise comparisons in a standardized format:
#
# datarestriction__variable1_vs_variable2
#
# The final standardized comparisons only contain alphanumeric values.
#
# Data restriction takes on values from differentiation states, 
# or all for global analyses.

library(synapseClient)
synapseLogin()
library(dplyr)
library(reshape2)

library(stringr)
library(tidyr)

# These are the covariates that will be generated pairwise,
# except C4_Cell_Line_ID (used as an ID and then removed)
FactorCovariates <- c('C4_Cell_Line_ID',
                      'Diffname_short', 'Cell_Line_Type', 'Cell_Line_of_Origin', 'Tissue_of_Origin', 
                      'Reprogramming_Gene_Combination', 'Culture_Conditions', 'Donor_Life_Stage',
                      'Race', 'Ethnicity' , 'Gender', 'Disease', 'Originating_Lab',
                      'Cell_Type_of_Origin', 'Cell_Type_of_Origin_Level2', 'Reprogramming_Vector_Type', 
                      'Reprogramming_Vector_Type_Level2', 'XIST_Level')

# For diffname_short values, in order
diffstates <- c("DE", "EB", "ECTO", "MESO-5", "MESO-15", "MESO-30", "SC")

# Taking this from the mRNA table because it has the most samples and has the 
# added Level2 covariates that are not present in cell line or sample proc
# metadata tables
res <- synTableQuery(sprintf("select * from syn3156503"))@values %>%
  filter(Diffname_short %in% diffstates) %>%
  select(one_of(FactorCovariates))

# Only for ordering the diffstates so MESO-5 comes first.
# Need to change them all back afterwards!
res <- res %>%
  mutate(Diffname_short=ifelse(Diffname_short == "MESO-5", "MESO-05", Diffname_short))

# Melt and get unique, get rid of anything N/A
res2 <- melt(res, id.vars = "C4_Cell_Line_ID") %>% 
  select(-C4_Cell_Line_ID) %>%
  unique() %>%
  filter(!(value %in% c("N/A", "")), !is.na(value))

# Merge with self on the covariate type
# Get rid of duplicates and self-comparisons
res3 <- res2 %>% 
  left_join(res2, by=c("variable")) %>%
  filter(as.character(value.x) < as.character(value.y))

# Change back to MESO-5
res3 <- res3 %>%
  mutate(value.x=ifelse(value.x == "MESO-05", "MESO-5", value.x),
         value.y=ifelse(value.y == "MESO-05", "MESO-5", value.y))

# Rename the colulmns, and create short versions of everything
# Short versions contain ONLY alphanumerics
res3 <- res3 %>%
  dplyr::rename(class=variable, variable1=value.x, variable2=value.y) %>%
  mutate(variable1Short=str_replace_all(variable1, '[^[:alnum:]]', ''),
         variable2Short=str_replace_all(variable2, '[^[:alnum:]]', ''))

# For the differentiation state comparisons, the data restriction is All
res3diff <- res3 %>% filter(class == "Diffname_short") %>% mutate(dataRestriction="All")

# For all other comparisons, the data restriction is All + all diffstates
res3rest <- res3 %>% filter(class != "Diffname_short")

origlen <- nrow(res3rest)
restdiffstates <- c("All", diffstates)
res3rest <- res3rest[rep(1:origlen, each=length(restdiffstates)), ] %>%
  mutate(dataRestriction=rep(restdiffstates, times=origlen))

# Combine the diffstate comparisons and the rest
res4 <- rbind(res3diff, res3rest) 

# Make the data restriction also column-friendly (only alphanumerics)
# Then create the final comparison column: datarestriction__variable1_vs_variable2
# Finally, reorder the columns
res4 <- res4 %>%
  mutate(dataRestrictionShort=str_replace_all(dataRestriction, '[^[:alnum:]]', '')) %>%
  unite(tmp, variable1Short, variable2Short, sep="_vs_", remove=FALSE) %>%
  unite(comparison, dataRestrictionShort, tmp, sep="__", remove=FALSE) %>%
  select(class, dataRestriction, dataRestrictionShort, variable1, variable1Short, 
         variable2, variable2Short, comparison, -tmp)

# Upload to a Synapse table
tc <- list(class=TableColumn(name="class", columnType="STRING", maximumSize=200),
           datarestriction=TableColumn(name="dataRestriction", columnType="STRING", maximumSize=200),
           datarestrictionshort=TableColumn(name="dataRestrictionShort", columnType="STRING", maximumSize=200),
           variable1=TableColumn(name="variable1", columnType="STRING", maximumSize=200),
           variable1short=TableColumn(name="variable1Short", columnType="STRING", maximumSize=200),
           variable2=TableColumn(name="variable2", columnType="STRING", maximumSize=200),
           variable2short=TableColumn(name="variable2Short", columnType="STRING", maximumSize=200),
           comparison=TableColumn(name="comparison", columnType="STRING", maximumSize=500))


schema <- TableSchema(name="Differential comparisons", parent="syn1773109", columns=tc)
tbl <- Table(tableSchema=schema, values=res4)
tbl <- synStore(tbl)
