# 0. HOUSEKEEPING ####
# clearing environment and memory
rm(list = ls())
graphics.off()

# loading packages
library(readxl)
library(tidyverse)
library(tableone)
library(countrycode)
library(AMR)
library(taxize)
library(purrr)
library(data.table)
library(future.apply)
library(parallel)
library(future)
library(rstan)
library(bayesplot)
library(mvtnorm)
library(posterior)
library(tibble)
library(ggplot2)
library(forcats)
library(cmdstanr)
library(writexl)
library(arsenal)
library(furrr)
library(sf)
library(scales)
library(rnaturalearth)
library(wbstats)
plan(multisession)

# setting working directory
setwd("C:/Users/jhitch/OneDrive - Nexus365/Vivli data challenges/Vivli 2025/Koen Pouwels's files - VIVLI_2/Student submission/Joint regression modelling")

# 1. READING IN ATLAS DATA ####
df <- read.csv("2025_03_11 atlas_antibiotics.csv",
               colClasses = rep("character", 134),
               na.strings = "")

# 2. DATASET STRUCTURE ####
dim(df)
head(df)
str(df)


### rename columns to lower case and shorten
df <- df %>%
  rename(
    id         = Isolate.Id,
    study      = Study,
    species    = Species,
    family     = Family,
    country    = Country,
    source     = Source,
    gender     = Gender,
    age        = Age.Group,
    speciality = Speciality,
    in_out     = In...Out.Patient,
    year       = Year,
    phenotype  = Phenotype
  )

# 3. PIE CHART SHOWING RAW BLOOD ISOLATE COUNTS BY COUNTRY ####
# filtering for blood isolates (ignoring blood vessel as a source)
df_blood <- df %>%
  filter(source == "Blood")

# counts of blood isolates by country
df_blood_country <- df_blood %>%
  count(country, name = "n") %>%
  # calculating percentage contribution of each country to total blood isolates
  mutate(
    perc  = n / sum(n),
    label = paste0(country, " (", scales::percent(perc, accuracy = 0.1), ")")
  )


# pie chart for blood isolates by country
ggplot(df_blood_country, aes(x = "", y = perc, fill = country)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Counts of blood isolates by country")

# 4. HEATMAP SHOWING RAW BLOOD ISOLATE COUNTS BY COUNTRY ####
# adding ISO3C code to data
df_blood_country <- df_blood_country %>%
  mutate(iso3c = countrycode(country, origin = "country.name", destination = "iso3c"))

# creating a country heat map for blood isolate counts
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  left_join(df_blood_country, by = c("iso_a3_eh" = "iso3c"))

heatmap <- ggplot(world) +
  geom_sf(aes(fill = n), color = "white", size = 0.2) +   # thin borders
  scale_fill_gradient(
    low = "#619cff",   # medium blue
    high = "#ff6f61",  # medium red
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "Global distribution of ATLAS blood isolates",
    fill = "Value"
  )

ggsave("blood_isolates_heatmap.png", plot = heatmap)

# 5. PIE CHART SHOWING RAW BLOOD ISOLATE COUNTS BY WORLD BANK INCOME GROUP ####
# getting World Bank income group data
countries <- wb_countries()
income <- countries[, c("iso3c", "country", "income_level")]

# adding ISO3C to blood isolate data
df_blood <- df_blood %>%
  mutate(iso3c = countrycode(country, origin = "country.name", destination = "iso3c"))

# joining income group to blood isolate data
df_blood <- df_blood %>%
  left_join(income, by = "iso3c")

df_blood %>% filter(income_level == "Not classified" | is.na(income_level)) %>% distinct(iso3c)
# Taiwan and Venezuela are the two countries with missing income group
# From a Google search, Taiwan is high-income
# Venezuela is genuinely currently unclassified by the World Bank due to a lack of economic data

# filling in the gaps for Taiwan
df_blood <- df_blood %>%
  mutate(income_level = ifelse(iso3c == "TWN", "High income", income_level))

# counts of blood isolates by income group
df_blood_income <- df_blood %>%
  group_by(income_level) %>%
  summarise(n = n()) %>%
  ungroup()

# removing countries that are not classified (this is just Venezuela)
df_blood_income <- df_blood %>%
  count(income_level, name = "n") %>%
  filter(income_level != "Not classified") %>%
  # adding each income group's percentage contribution to total blood isolates
  mutate(
    perc  = n / sum(n),
    label = paste0(income_level, " (", scales::percent(perc, accuracy = 0.1), ")")
  )


# creating a pie chart for raw blood isolate counts by World Bank income group
ggplot(df_blood_income, aes(x = "", y = perc, fill = income_level)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Counts of isolates by World Bank income group")


# checking for duplicates
sum(duplicated(df$id)) # no duplicates

# checking for duplicates based on all columns
# sum(duplicated(df$id))  # unhash if needed
# no duplicate isolates found

## recode as na
df <- df %>%
  mutate(
    gender     = na_if(gender, "-"),
    age        = na_if(age, "Unknown"),
    speciality = na_if(speciality, "None Given"),
    speciality = na_if(speciality, "Other"),
    source     = na_if(source, "None Given"),
    in_out     = na_if(in_out, "None Given"),
    in_out     = na_if(in_out, "Other")
  )

#recode as factor
df <- df %>%
  mutate(across(c(study, family, country, gender, speciality, source, in_out), as.factor))


# variables ending "_I" are SIR e.g. "Amikacin_I"
# variables without "_I" are MIC e.g. "Amikacin"
# not subsetting out MIC columns yet because we may create new SIR columns by 
# interpreting the MICs ourselves using the AMR package

# assume one isolate per patient (not sure we can validate this)
# we only have a few patient characteristics e.g. age group, gender, ward

# 6. INSPECTING AND REFORMATTING ISOLATE/PATIENT METADATA VARIABLES ####
# 2.1 Study
# 


# isolate counts and percentages by study
n_study <- df %>%
  count(study, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_study)


# checking duration of each study
n_study_country <- df %>%
  group_by(study) %>%
  summarise(start = min(year), end = max(year), country = unique(country), .groups = "drop")
print(n_study_country)
# SPIDAAR is just Uganda, Malawi, Uganda and Kenya 2021-2023


# 2.2 Species

# checking number of unique species
species <- unique(df$species)
length(species)

# using the AMR package to get a range of potentially useful microorganism properties
# species names and gram stain are likely to be the two most important ones
df <- df %>%
  mutate(mo_amr = as.mo(species),
         short_name_amr = mo_shortname(species),
         full_name_amr = mo_fullname(species),
         species_amr = mo_species(species),
         genus_amr = mo_genus(species),
         family_amr = mo_family(species),
         order_amr = mo_order(species),
         class_amr = mo_class(species),
         phylum_amr = mo_phylum(species),
         kingdom_amr = mo_kingdom(species),
         domain_amr = mo_domain(species),
         type_amr = mo_type(species),
         pathogenicity_amr = mo_pathogenicity(species),
         gramstain_amr = mo_gramstain(species),
         gram_negative_amr = mo_is_gram_negative(species),
         gram_positive_amr = mo_is_gram_positive(species),
         oxygen_tolerance_amr = mo_oxygen_tolerance(species),
         anaerobic_amr = mo_is_anaerobic(species))

# calculating counts and percentages by bacteria family
n_family_amr <- df %>%
  count(family_amr, name = "count") %>%
  # calculating isolate counts and percentages by family
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
head(n_family_amr)

sum(is.na(df$family)) 


# 2.4 country
length(unique(df$country)) 

# calculating counts and percentages by country
n_country <- df %>%
  count(country, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_country)


# 2.5 gender
# isolate counts by gender
n_gender <- df %>%
  count(gender, name = "count") %>%
  ## calculate percentages
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_gender)
# large number of NAs for gender. I will deal with these later

# 2.6 age group
unique(df$age)

# converting age group to ordered factor
df$age_factor <- factor(df$age,
                        levels = c("0 to 2 Years", "3 to 12 Years", "13 to 18 Years"
                                   ,"19 to 64 Years", "65 to 84 Years", "85 and Over"),
                        ordered = TRUE)

# calculating isolate counts and percentages by age group
n_age <- df %>%
  count(age, name = "count") %>%
  mutate(percent = round(100*n()/nrow(df), 1)) %>%
  arrange(desc(count))
print(n_age)



# 2.7 speciality
# calculating isolate counts and percentages by speciality
n_speciality <- df %>%
  count(speciality, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_speciality)



# 2.8 source

# isolate counts and percentages by source
n_source <- df %>%
  count(source, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_source)


# 2.9 in/outpatient

# calculating isolate counts and percentages by inpatient/outpatient
n_in_out <- df %>%
  count(in_out, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_in_out)

# 2.10 year

# calculating isolate counts and percentages by country
n_year <- df %>%
  count(year, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_year)


# 2.11 phenotype (probably won't use this)

# calculating isolate counts and percentages by phenotype
n_phenotype <- df %>%
  count(phenotype, name = "count") %>%
  mutate(percent = round(100 * count / nrow(df), 1)) %>%
  arrange(desc(count))
print(n_phenotype)
# high proportion of missing values

# calculating proportions of isolates with missing values for the different
# isolate/patient metadata variables
n_missing_metadata <- df %>%
  select(1:13, family_amr, gramstain_amr) %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "na_percent") %>%
  arrange(desc(na_percent))
print(n_missing_metadata)

# 7. CREATING VECTORS FOR MIC AND SIR COLUMN NAMES ####
# Naming system: 
# "Amikacin" = MIC
# "Amikacin_I" = SIR

# creating vectors for all MIC and SIR column names
mic_col_names <- grep("^[A-Z].*[^_I]$", colnames(df), value = TRUE)
mic_col_names <- mic_col_names[mic_col_names %in% colnames(df)[14:110]]
sir_col_names <- grep("_I$", colnames(df), value = TRUE)

# 8. INITIAL FILTERS FOR ANALYTICAL SAMPLE ####
# The interpretation of MIC values takes a long time so it is not feasible to
# do this for all isolates. Instead, we apply some initial filters here based on
# our planned analytical sample and interpret MIC values for this subset of
# isolates

# We focus on possible bloodstream infections caused by any gram-negative pathogen
# in adult patients in the period 2018-2023
df_1 <- df %>%
  filter(gram_negative_amr == TRUE, # gram negatives
         source == "Blood", # blood isolates
         year %in% c(2018:2023), # study period
         age %in% c("19 to 64 Years", "65 to 84 Years", "85 and Over", NA)) # adults

# duplicating to preserve raw filtered data table
df_analysis <- df_1

# 9. INTERPRETING MIC VALUES INTO SIR USING AMR PACKAGE ####

# interpreting MIC columns in terms of SIR for relevant antibiotics
# declaring subset of antibiotics
keep <- c("Meropenem", "Piperacillin.tazobactam", "Ampicillin",
          "Tigecycline", "Gentamicin", "Amikacin", "Ceftazidime", "Levofloxacin",
          "Colistin", "Amoxycillin.clavulanate")
mic_col_names <- mic_col_names[mic_col_names %in% keep]

# investigating intrinsic resistance
# now adding column for intrinsic resistance (TRUE or FALSE) based on EUCAST 
# for each antibiotic
for (ab in mic_col_names) {
  # Create column name for intrinsic resistance
  col_name <- paste0(ab, "_intrinsic_resistant")
  
  # Apply the function and create the column
  df_analysis[[col_name]] <- mo_is_intrinsic_resistant(
    x = df_analysis$mo_amr,
    ab = ab
  )
}
summary(df_analysis)

# counting cases where SIR is missing but intrinsic resistance = TRUE meaning
# that the SIR value should really be R
missing_SIR_counts <- sapply(mic_col_names, function(ab) {
  sir_col <- paste0(ab, "_I")
  intrinsic_col <- paste0(ab, "_intrinsic_resistant")
  
  sum(is.na(df_analysis[[sir_col]]) & df_analysis[[intrinsic_col]] == TRUE)
})

# View result
missing_SIR_counts

# interpreting MIC values in terms of SIR using EUCAST 2025 guidelines
# firstly looking at pip-tazo as an example
# A. without intrinsic resistance
df_analysis$Piperacillin.tazobactam_SIR_EUCAST <-
  as.sir(x = as.mic(df_analysis$Piperacillin.tazobactam),
         mo = df_analysis$mo_amr,
         ab = as.ab("Piperacillin.tazobactam"),
         guideline = "EUCAST 2025",
         uti = FALSE,
         add_intrinsic_resistance = FALSE)

# isolate counts by new SIR value for pip/tazo
df_analysis %>%
  count(Piperacillin.tazobactam_SIR_EUCAST)

# isolate counts by original SIR value for pip/tazo
df_analysis %>%
  count(Piperacillin.tazobactam_I)

# B. now with intrinsic resistance
df_analysis$Piperacillin.tazobactam_SIR_EUCAST_IR <-
  as.sir(x = as.mic(df_analysis$Piperacillin.tazobactam),
         mo = df_analysis$mo_amr,
         ab = as.ab("Piperacillin.tazobactam"),
         guideline = "EUCAST 2025",
         uti = FALSE,
         add_intrinsic_resistance = TRUE)

# isolate counts by new SIR value for pip/tazo with intrinsic resistance
df_analysis %>%
  count(Piperacillin.tazobactam_SIR_EUCAST_IR)
# accounting for intrinsic resistance doesn't change anything for pip/tazo but
# this might not be the case for Linezolid

# without intrinsic resistance
df_analysis$Linezolid_SIR_EUCAST <-
  as.sir(x = as.mic(df_analysis$Linezolid),
         mo = df_analysis$mo_amr,
         ab = as.ab("Linezolid"),
         guideline = "EUCAST 2025",
         uti = FALSE,
         add_intrinsic_resistance = FALSE)

# isolate counts by SIR for Linezolid
df_analysis %>%
  count(Linezolid_SIR_EUCAST)

# with intrinsic resistance
df_analysis$Linezolid_SIR_EUCAST_IR <-
  as.sir(x = as.mic(df_analysis$Linezolid),
         mo = df_analysis$mo_amr,
         ab = as.ab("Linezolid"),
         guideline = "EUCAST 2025",
         uti = FALSE,
         add_intrinsic_resistance = TRUE)

# isolate counts by SIR for Linezolid with intrinsic resistance
df_analysis %>%
  count(Linezolid_SIR_EUCAST_IR)

# for linezolid, we go from having mostly NAs with a few S and R to almost all
# having R (because of intrinsic resistance)
# it seems that gram-negatives are intrinsically resistant to linezolid


# now interpreting MIC values for all antibiotics we want to include in the study
# A. Not taking into account intrinsic resistance
for (ab in mic_col_names) {
  mic_col <- ab
  sir_col <- paste0(ab, "_SIR_EUCAST")
  
  # Apply as.sir to get SIR values
  df_analysis[[sir_col]] <- as.sir(
    x = as.mic(df_analysis[[mic_col]]),
    mo = df_analysis$mo_amr,
    ab = as.ab(ab),
    guideline = "EUCAST 2025",
    uti = FALSE,
    add_intrinsic_resistance = FALSE
  )
}

# B. Taking into account intrinsic resistance
for (ab in mic_col_names) {
  mic_col <- ab
  sir_col <- paste0(ab, "_SIR_EUCAST_IR")
  
  # Apply as.sir to get SIR values
  df_analysis[[sir_col]] <- as.sir(
    x = as.mic(df_analysis[[mic_col]]),
    mo = df_analysis$mo_amr,
    ab = as.ab(ab),
    guideline = "EUCAST 2025",
    uti = FALSE,
    add_intrinsic_resistance = TRUE
  )
}

# 10. CHECKING MISSINGNESS FOR AMR VARIABLES ####
# We now have two additional SIR columns for each antibiotic. The first takes
# the MIC values in the raw ATLAS dataset and interprets them into SIR values
# based on latest EUCAST guidelines and ignoring intrinsic resistance. The
# second does the same but takes into account intrinsic resistance i.e. if
# the microorganism is known to be intrinsically resistant to a given 
# antibiotic then it will be classified as resistant regardless of the MIC value
# in the data

# analysing missingness in SIR values (this may inform choice of antibiotics)
# without intrinsic resistance
na_proportions <- sapply(
  df_analysis[grep("_SIR_EUCAST$", names(df_analysis))],
  function(x) 100*mean(is.na(x))
)
na_proportions <- as.data.frame(na_proportions)
na_proportions <- na_proportions %>%
  filter(na_proportions < 3)

# with intrinsic resistance
na_proportions_IR <- sapply(
  df_analysis[grep("_SIR_EUCAST_IR$", names(df_analysis))],
  function(x) 100*mean(is.na(x))
)
na_proportions_IR <- as.data.frame(na_proportions_IR)
na_proportions_IR <- na_proportions_IR %>%
  filter(na_proportions_IR < 3)

nrow(na_proportions)
nrow(na_proportions_IR)
# linezolid is intrinsically resistant to most gram-negatives
# We will exclude Linezolid from the analysis

# 11. CREATING BINARY RESISTANT VS NON-RESISTANT VARIABLES ####
# creating binary variables for resistant vs non-resistant
# without intrinsic resistance
df_analysis <- df_analysis %>%
  mutate(across(
    ends_with("_SIR_EUCAST"),
    ~ if_else(
      as.character(.) %in% c("R", "I"),
      1L,
      if_else(is.na(.), NA_integer_, 0L)
    ),
    .names = "{.col}_bin"
  )) %>%
  rename_with(
    ~ sub("_SIR_EUCAST_bin$", "_bin", .x),
    ends_with("_SIR_EUCAST_bin")
  )

# with intrinsic resistance
df_analysis <- df_analysis %>%
  mutate(across(
    ends_with("_SIR_EUCAST_IR"),
    ~ if_else(
      as.character(.) %in% c("R", "I"),
      1L,
      if_else(is.na(.), NA_integer_, 0L)
    ),
    .names = "{.col}_bin_IR"
  )) %>%
  rename_with(
    ~ sub("_SIR_EUCAST_IR_bin_IR$", "_bin_IR", .x),
    ends_with("_SIR_EUCAST_IR_bin_IR")
  )


# creating binary ICU vs non-ICU variable
df_analysis <- df_analysis %>%
  mutate(icu = case_when(
    is.na(speciality) ~ NA_integer_,
    speciality %in% c("Medicine ICU", "Pediatric ICU", "Surgery ICU") ~ 1L,
    TRUE ~ 0L
  ))
summary(df_analysis$icu)

# checking new variables
# antibiotic with assumed low resistance
#summary(df_analysis$Meropenem_bin) # unhash if needed
#summary(df_analysis$Meropenem_bin_IR) # unhash if needed

# antibiotic with assumed high resistance 
#summary(df_analysis$Ampicillin_bin) # unhash if needed
#summary(df_analysis$Ampicillin_bin_IR) # unhash if needed



# 12. SAVING INTERMEDIATE DATA TABLE ####
# saving this intermediate data table as R and csv files for purpose of descriptive
# analysis. This data table contains SIR, MIC and binary resistant vs non-resistant
# classifications for gram-negative blood isolates and a set of relevant antibiotics
# for bloodstream infections collected in adults 2018-2022. We can filter for 
# specific countries and complete cases later

save(df_analysis, file = "vivli_clean_data_v4.RData")
write.csv(df_analysis, file = "vivli_clean_data_v4.csv", row.names = FALSE)

# 13. REFORMATTING THE DATA SPECIFICALLY FOR MODELLING ####
# creating new version of the intermediate data table to preserve raw intermediate 
# data file
df <- df_analysis
raw <- df_analysis
rm(df_analysis)

# adding ISO3C based on country name
df <- df %>%
  mutate(iso3c = countrycode(country, origin = "country.name", destination = "iso3c"))

# isolate counts by country
country_plot <- df %>%
  group_by(iso3c) %>%
  summarise(n = n()) %>%
  ungroup()

# checking proportions of isolates with NAs for resistance by country ## unhash if needed
# df %>%
#   group_by(country) %>%
#   summarise(across(ends_with("_bin"), ~ 100*mean(is.na(.)), .names = "prop_na_{.col}")) %>%
#   summary()

# ignoring linezolid, the maximum rate of missingness for a given country is 6%
# i.e. almost all gram-negative blood isolates are tested against all of our
# antibiotics of interest

# for simplicity, we want to collapse the age groups, combining 64-84 and >85 
# into an "65 and Over" category
summary(df$age)
# converting age group to character
df$age <- as.character(df$age)
df$age[df$age %in% c("65 to 84 Years", "85 and Over")] <- "65 and Over"
df %>%
  count(age)

# filtering for complete cases based on model variables
df <- df[complete.cases(df[, c("id", 
                               "country", 
                               "year",
                               "gender", 
                               "age", 
                               "icu",
                               "Meropenem_bin", 
                               "Piperacillin.tazobactam_bin",
                               "Ampicillin_bin", 
                               "Tigecycline_bin",
                               "Gentamicin_bin", 
                               "Amikacin_bin", 
                               "Ceftazidime_bin",
                               "Levofloxacin_bin", 
                               "Colistin_bin", 
                               "Amoxycillin.clavulanate_bin")]), ]
# 38,524 complete case isolates


# adding region and income group variables
df <- df %>%
  mutate(region = countrycode(country, origin = "country.name", destination = "un.regionsub.name"))
# this didn't work for Taiwan

countries <- wb_countries()
income <- countries[, c("iso3c", "country", "income_level")]

# joining income group
df <- merge(df, income, by = "iso3c", all.x = TRUE)

# counts of isolates by region
df %>%
  count(region, sort = TRUE)
# Southern Europe is the region with the most isolates

# creating a separate data table for isolates from Southern Europe
df_south_europe <- df %>%
  filter(region == "Southern Europe")

# calculating isolate counts by country within Southern Europe
n_south_europe_country <- df_south_europe %>%
  group_by(country.x) %>%
  summarise(n = n()) %>%
  ungroup()

# adding percentages
n_south_europe_country <- n_south_europe_country %>%
  mutate(perc = n/sum(n),
         label = paste0(country.x, " (", percent(perc, accuracy = 0.1), ")"))

# pie chart for counts of isolates by country within Southern Europe
ggplot(n_south_europe_country, aes(x = "", y = perc, fill = country.x)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Counts of Enterobacteriaceae blood isolates by country")

# looking at counts of isolates by species
df_south_europe %>%
  count(short_name_amr, sort = TRUE)

# creating a data table just model variables
df_save <- df_south_europe %>%
  dplyr::select("full_name_amr",
                "family_amr",
                "country.x",
                "region",
                "year",
                "gender", 
                "age", 
                "icu",
                "Meropenem_bin", 
                "Piperacillin.tazobactam_bin",
                "Ampicillin_bin", 
                "Tigecycline_bin",
                "Gentamicin_bin", 
                "Amikacin_bin", 
                "Ceftazidime_bin",
                "Levofloxacin_bin", 
                "Colistin_bin", 
                "Amoxycillin.clavulanate_bin") %>%
  rename(country = country.x)

# renaming some variables
df_save <- df_save %>%
  rename(icu_bin = icu,
         pathogen = full_name_amr,
         family = family_amr)

n_country <- df_save %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  ungroup()

# 14. SAVING FINAL DATA TABLE FOR DESCRIPTIVE ANALYSIS AND MODELLING ####
summary(df_save)
saveRDS(df_save, "vivli_data_final.RDS")

