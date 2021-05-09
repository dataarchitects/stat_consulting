######################### loading the required libraries ########################################################
library(dplyr)
library(lubridate)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(psych)
library(mclust)
library(HDclassif)
library(tidyr)
#################################################################################################################
#setting seed to ensure reproducability
set.seed(123)
#reading the dataset of vaccination data
data = read.csv("/home/shourya/stat-consulting-report/corona-data.csv")
#selecting columns required for the study
vaccine_data =  data[,c("date","location","continent","reproduction_rate","icu_patients","hosp_patients","total_vaccinations","people_fully_vaccinated","people_vaccinated_per_hundred","population","population_density","aged_65_older","gdp_per_capita","life_expectancy","extreme_poverty","hospital_beds_per_thousand")]
#ensuring the date column is correctly formatted
vaccine_data$date <- ymd(vaccine_data$date)
#extracting year and month from the date column
vaccine_data$year <- year(vaccine_data$date)
vaccine_data$month <- month(vaccine_data$date)
#filling null values of vaccination with 0
vaccine_data$total_vaccinations[is.na(vaccine_data$total_vaccinations)] <- 0
#filling null values of people fully vaccinated with 0
vaccine_data$people_fully_vaccinated[is.na(vaccine_data$people_fully_vaccinated)] <- 0
#filling the null values of available bed with 0
vaccine_data$hospital_beds_per_thousand[is.na(vaccine_data$hospital_beds_per_thousand)] <- 0
###############################################################################################################

########################## creating dataset for average monthly vaccinations #################################

cumm_data <- vaccine_data %>% group_by(location) %>% summarize(Sum_Vaccinations = sum(total_vaccinations, na.rm=TRUE), Sum_tot_Vacc = sum(people_fully_vaccinated, na.rm = TRUE))
cumm_data$Sum_Vaccinations <- cumm_data$Sum_Vaccinations/15
cumm_data$Sum_tot_Vacc <- cumm_data$Sum_tot_Vacc/15
join_cols <- vaccine_data[,c("location","life_expectancy","gdp_per_capita","population_density","population","aged_65_older","people_vaccinated_per_hundred")]
join_cols <- join_cols[!duplicated(join_cols$location),]
cumm_data <- merge(x = cumm_data, y = join_cols, by = "location", all.x = TRUE)
cumm_data[is.na(cumm_data)] <- 0
#cumm_data <- na.omit(cumm_data)
rownames(cumm_data) <- cumm_data$location
keeps <- "location"
cumm_data <- subset(cumm_data, select = -c(location))
#remove the observation world and from continents
row.names.remove <- c("World","Africa","Asia","Europe","North America","Oceania","South America", "European Union")
cumm_data <- cumm_data[!(row.names(cumm_data) %in% row.names.remove), ]

###############################################################################################################

################################################# question 1 ###################################################
#selecting data from target countries in the study
target <- c("Belgium","India","Brazil","United States","South Africa")
fil_cn <- filter(vaccine_data, location %in% target)
#filling nulls by latest available value for percentage of population vaccinated
fil_cn <- fil_cn %>% group_by(location) %>% fill(people_vaccinated_per_hundred, .direction = "down")
#creating variable for cumulitive figure of individual vaccines administered
plot_data <- fil_cn %>% group_by(location) %>% mutate(cum_vac = cumsum(total_vaccinations))
#creating the plot of cummulitive figure of individual vaccines administered per country
ggplot(data = plot_data, aes(x = date, y = cum_vac, color = location)) + geom_line() + ggtitle("Progress of Vaccinations") + xlab("Date") + ylab("Cumulitive Vaccinations")
#creating the plot for percentage of population vaccinated
ggplot(data = fil_cn, aes(x = date, y = people_vaccinated_per_hundred, color = location)) + geom_line() + ggtitle("Progress of Vaccinations") + xlab("Date") + ylab("% population vaccinated")
#cummulitive sum of vaccinations
summed_data <- plot_data %>% group_by(location) %>% summarise(sum_tot_vac = sum(total_vaccinations))
#latest information about percentage of population vaccinated
percentage_data <- plot_data %>% arrange(location,people_vaccinated_per_hundred) %>% group_by(location) %>% summarise_all(last)
#regression analysis of total vaccinations
vac_effect <- lm(Sum_Vaccinations ~ life_expectancy + gdp_per_capita + population_density + aged_65_older, data = cumm_data)
summary(vac_effect)

######################################## question 2 ############################################################
#estimating the number of clusters appropiate for the dataset
wssplot <- function(data, nc=15, seed=123){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of groups",
       ylab="Sum of squares within a group")}
wssplot(cumm_data, nc = 20)

#creating k means clustering object and visualization of the solution
res.km <- kmeans(scale(cumm_data), 3, nstart = 20)
fviz_cluster(res.km, data = cumm_data)
#creating a dataframe of cluster number
clus_join <- cumm_data
clus_join$location <- rownames(cumm_data)
clus_df <- data.frame(res.km$cluster)
clus_df$location <- rownames(clus_df)
clus_df <- merge(x = clus_df, y = clus_join, by = "location", all.x = TRUE)
#filtering dataset for the first clustering and summarizing required information
clus_1 <- clus_df %>% filter(res.km.cluster==1) %>% summarise(mean_ma = mean(Sum_Vaccinations), mean_tv = mean(Sum_tot_Vacc), mean_gdp = mean(gdp_per_capita.y), mean_pd = mean(population_density.y), mean_le = mean(life_expectancy.y))
#filtering dataset for the second clustering and summarizing required information
clus_2 <- clus_df %>% filter(res.km.cluster==2) %>% summarise(mean_ma = mean(Sum_Vaccinations), mean_tv = mean(Sum_tot_Vacc), mean_gdp = mean(gdp_per_capita.y), mean_pd = mean(population_density.y), mean_le = mean(life_expectancy.y))
#filtering dataset for the third clustering and summarizing required information
clus_3 <- clus_df %>% filter(res.km.cluster==3) %>% summarise(mean_ma = mean(Sum_Vaccinations), mean_tv = mean(Sum_tot_Vacc), mean_gdp = mean(gdp_per_capita.y), mean_pd = mean(population_density.y), mean_le = mean(life_expectancy.y))

########################################## Question 3 #########################################################
vaccine_data$icu_patients[is.na(vaccine_data$icu_patients)] <- 0
vaccine_data$hosp_patients[is.na(vaccine_data$hosp_patients)] <- 0
################ global effect of vaccination on ICU patients and hospitalizations ##########################

summed_global <- vaccine_data %>% group_by(date) %>% summarise(sum_hosp = sum(hosp_patients), sum_icu = sum(icu_patients), sum_vac = cumsum(total_vaccinations))

mod_glo_icu <- lm(sum_icu ~ sum_vac, data = summed_global)
summary(mod_glo_icu)
mod_glo_hosp <- lm(sum_hosp ~ sum_vac, data = summed_global)
summary(mod_glo_hosp)

######################## effect of vaccination in belgium ###################################################

belg_data <- vaccine_data %>% filter(location == "Belgium")
summed_belg <- belg_data %>% group_by(date) %>% summarise(sum_hosp = sum(hosp_patients), sum_icu = sum(icu_patients), sum_vac = cumsum(total_vaccinations))

bel_glo_icu <- lm(sum_icu ~ sum_vac, data = summed_belg)
summary(bel_glo_icu)
bel_glo_hosp <- lm(sum_hosp ~ sum_vac, data = summed_belg)
summary(bel_glo_hosp)

#studying effect of percentage population vaccinated on hospitalization and icu cases
belg_data$hosp_patients[is.na(belg_data$hosp_patients)] <- 0
belg_data$people_vaccinated_per_hundred[is.na(belg_data$people_vaccinated_per_hundred)] <- 0
#belg_data$people_vaccinated_per_hundred[is.na(belg_data$people_vaccinated_per_hundred)] <- 0
belg_data <- belg_data %>% fill(people_vaccinated_per_hundred, .direction = "down")

bel_hosp_per <- lm(hosp_patients ~ people_vaccinated_per_hundred, data = belg_data)
summary(bel_hosp_per)

bel_icu_per <- lm(icu_patients ~ people_vaccinated_per_hundred, data = belg_data)
summary(bel_icu_per)


##################################### effect of vaccination in israel ###################################

is_data <- vaccine_data %>% filter(location == "Israel")
summed_is <- is_data %>% group_by(date) %>% summarise(sum_hosp = sum(hosp_patients), sum_icu = sum(icu_patients), sum_vac = cumsum(total_vaccinations))
is_data <- is_data %>% fill(people_vaccinated_per_hundred, .direction = "down")


is_glo_icu <- lm(sum_icu ~ sum_vac, data = summed_is)
summary(is_glo_icu)
is_glo_hosp <- lm(sum_hosp ~ sum_vac, data = summed_is)
summary(is_glo_hosp)

is_per_icu <- lm(icu_patients ~ people_vaccinated_per_hundred, data = is_data)
summary(is_per_icu)

is_per_hosp <- lm(hosp_patients ~ people_vaccinated_per_hundred, data = is_data)
summary(is_per_hosp)
