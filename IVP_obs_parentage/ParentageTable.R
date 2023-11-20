# ---------------
# Title: Parentage Table
# Date: 17 Nov 2023
# Author: mgranellruiz
# Goal: My goal is to end up with a table that has 1 individual per row and two columns of potential mathers and fathers.
# with the date of birth and possible date of conception
# ---------------

# library ---------------------
library(lubridate)
library(dplyr)
library(ggplot2)
library(lme4)
library(ggstatsplot)
library(fitdistrplus)
library(gamlss)
library(ggside)
library(ggpubr)
library(stringr)
library(gridExtra)
library(ggtext)
library(tidyr)
library(patchwork)
source('/Users/mariagranell/Repositories/data/functions.R')

# path ------------------------
setwd("/IVP_obs_parentage")

# data ------------------------
lifehistory <- read.csv("/Users/mariagranell/Repositories/data/life_history/tbl_Creation/TBL/factchecked_LH_260523.csv")

# First I´ve noticed somethimes I don´t have a AnimalCode, since for this I just need an identifyer
# when they don´t have a code I´ll put AnimalName without spaces.
lh <- lifehistory %>% mutate(AnimalCode = ifelse(is.na(AnimalCode), str_replace(AnimalName, " ",""), AnimalCode))
lh <- lh %>% mutate(AnimalName = ifelse(AnimalCode == "Nil", "Nile", AnimalName)) # Nil is Nile in full name

# INDIVIDUALS ----------------
# Firstly I am only going to select the individuals for which we have a BirthGroup, those are the ones for which we might
# know their mum/dad
# This is the dataframe that we are going to use as a template. It will contain4 collumns
# AnimalName: full name
# AnimalCode: 3 or 4 letter identifier, 3 for males 4 for females
# Mother: if the mother is known for the individual
# FatherGenetics: done by Franca during her master
# PossibleFathers: adult males that could have been present duringthe time of conception in the group
# DOBestimate: Estimated DOB

indv <- lifehistory %>% filter(Tenure_type == "BirthGroup", !is.na(FirstDate)) %>%
  mutate(Year_born = as.numeric(year(DOB_estimate))) %>%
  rename(FatherGenetics = Father) %>%
  dplyr::select(AnimalName,AnimalCode, Mother, FatherGenetics, DOB_estimate, FirstDate, DOB, Group_mb, Tenure_type, Year_born)

nrow(indv %>% distinct(AnimalName))
# 645 indiv, non of them with the name duplicated

# out of curiosity, which year was the best in terms of births? This is very dirty do not trust
indv %>% mutate(year = as.numeric(year(DOB_estimate))) %>% ggplot(aes(x = year)) + geom_histogram(binwidth = 1)

# MOM -----------------------------------
# how many of them have already a mom?
nrow(indv%>% filter(!is.na(Mother)))
# 481 indv have a mom

# lets have a look to the individuals that don´t have a mother
# they cannot be resolved, however you can see from the histogram most of the are
# before 2016, thus the beggining of the project or from CR or IF. Thus is normal we don´t know
nomum <- indv%>% filter(is.na(Mother)) %>% mutate(year = as.numeric(year(FirstDate)))
nomum  %>% ggplot(aes(x = year)) + geom_histogram(binwidth = 1)
#View(nomum %>% filter(year > 2016))

# DAD GENETICS -----------------------------------
# how many of them have already a dad?
nrow(indv%>% filter(!is.na(FatherGenetics)))
# 75 indv have a dad

# These dads should have been dicovered when franca did her thesis
# this is correct only dad for babies born between 2010 and 2015
indv%>% filter(!is.na(FatherGenetics)) %>% mutate(year = as.numeric(year(FirstDate)))%>%
  ggplot(aes(x = year)) + geom_histogram(binwidth = 1)

# POSSIBLE DADYS -----------------------------------

# Now the fun part starts. I am going to make a list of the years and,
# and a list of all the males older than 3 year old present in that group during the mating season for that year.
# Broad mating season (MS) definition: from March to August.

# list of year based on the individuals
year_list <- indv %>% mutate(Year_born = as.numeric(year(FirstDate))) %>% distinct(Year_born)

# MS period, between March and August
begginingMS <- "-03-01"
endingMS <- "-08-30"

# end table
ProbFathers <- data.frame()

for (i in seq_len(nrow(year_list))){
    year <- year_list[i,1]

    aa <- lh %>%
      mutate(age_MS = add_age(DOB_estimate, date = paste0(year, begginingMS), unit = "Years")) %>%
      filter(
        Sex == "M", # only males
        age_MS > 3, # older than 3 years during the MS
        StartDate_mb < paste0(year, endingMS), # they start in the group before the MS ends
        EndDate_mb > paste0(year, begginingMS) # they leave the gp after the beggining of the MS
      ) %>%
      group_by(Group_mb) %>%
      summarize(Year = year, PotentialFathers = paste(AnimalCode, collapse = ";"))

    # Combine the results
    ProbFathers <- rbind(ProbFathers, aa)
}
rm(year, begginingMS, endingMS, year_list, aa, i)

# Now we have to combine the ProbFathers df with the Individuals
head(ProbFathers)
head(indv)

# Merge dataframes
tbl_parentage<- indv %>% left_join(ProbFathers, by = c("Year_born" = "Year", "Group_mb")) %>%
  dplyr::select(AnimalName, AnimalCode, Mother, FatherGenetics, PotentialFathers, DOB, DOB_estimate, Group_mb)

#write.csv(tbl_parentage, "/Users/mariagranell/Repositories/hyRAD/IVP_obs_parentage/tbl_obs_parentage.csv", row.names = F)

## CURIOSITY CHECK

# FATHERS IN POTENTIAL FATHERS?

ff <- tbl_parentage %>% filter(!is.na(FatherGenetics)) %>%
  left_join(lh %>% dplyr::select(AnimalName, AnimalCode) %>% distinct(.)
              %>% rename(FatherCode = AnimalCode),
                 by = c("FatherGenetics" = "AnimalName"))

# how many Father from Genetics we have that have no code
nrow(ff %>% filter(!is.na(FatherGenetics) & is.na(FatherCode)))
# 17 out of 75 individuals have a father with no code name.

# who are they? two individuals. Sanbonani (no match) and Nile (probably Nil) # Nil-Nile corrected!
print(ff %>% filter(!is.na(FatherGenetics) & is.na(FatherCode)) %>% distinct(FatherGenetics))

# How many fathers appear in the potential fathers
ff %>% mutate(poteualgen = ifelse(str_detect(PotentialFathers, FatherCode), "yes", "no")) %>%
  group_by(poteualgen) %>% summarize( n = n())
View(ff %>% mutate(poteualgen = ifelse(str_detect(PotentialFathers, FatherCode), "yes", "no")))
# only two cases when the father was not named in the potential fathes. For Vin and Zink
# Seems like for Vin, Voldemort was visiting LT before we spotted him. which make sense. LT is not very habituated
# For Zink, Yst apparently never went to KB but he must have!

# WHICH DATES MOST INDV WERE BORN
# Mos individuals were born between Oct-Dec. If we assume 7 monts pregnancy, most conceptions happen March-May
tbl_parentage %>% filter(!is.na(DOB)) %>%
  mutate(birth = ymd(paste("2021", month(DOB), day(DOB), sep = "-")),
                         conception =ymd(paste("2021", month(DOB) -7, day(DOB), sep = "-"))) %>%
  pivot_longer(cols = c("birth", "conception"), names_to = "Event", values_to = "Dates") %>%
   ggplot(aes(x = Dates, fill = Event)) + geom_histogram(binwidth = 1) +
   coord_flip() +
   scale_x_date(date_breaks = "1 month", date_labels = "%B")   # Set breaks and labels for each month


# HOW MANY INDIVIDUALS DO WE WANT TO SEQUENCE
# 155 Babies have been borned since 2021-01-01
nrow(lh %>% filter(DOB_estimate > 2021-01-01, Tenure_type == "BirthGroup"))
# 438 individuals have been part of the IVP since 2021-01-01
nrow(lh %>% filter(EndDate_mb > 2021-01-01) %>% distinct(AnimalCode))




