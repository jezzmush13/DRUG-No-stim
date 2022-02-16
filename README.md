# DRUG-No-stim
Most recent code for analysis

#GCaMP Calcium Analysis DRUG-No stim


#load the library 
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)

##################################################################################
# SECTION 1 OF THE CODE, CLEANING THE DATA

#load the data 
filenames <- list.files(pattern="*.csv") #read all the .csv files 


#create a dataframe; Remove missing values
df<-purrr::map_df(filenames, read_csv, .id = 'filename') 
df<- na.omit(df) #remove missing values 

#change the names to lowercase
df$roiName<- tolower(df$roiName)
df$Spot<- tolower(df$Spot)
df$Condition<- tolower(df$Condition)
df$drug<- tolower(df$drug)

#Change nominal and date format
df$filename = factor(df$filename)
df$roiName = factor(df$roiName)
df$Animal = factor(df$Animal)
df$drug = factor(df$drug)
df$Spot = factor(df$Spot)
df$filename = factor(df$filename)
df$date = factor(df$date)
view(df)

#Modifying variables and cleaning data frame
df_short <- df[ ,c('date', 'drug', 'Condition', 'Spot', 'Animal', 'roiName', 'amplitude', 'halfWidth','peakAUC')]
#Replace after drugs for baseline
df_short$drug <- str_replace_all(df_short$drug, 'after drugs', 'baseline')
#adding unique Spot
df_short$Unique_Spot <-paste(df_short$Animal, df_short$Spot, sep= "_")
#adding unique ROI name
df_short$Unique_ROIname <-paste(df_short$Animal, df_short$Spot, df_short$roiName, sep= "_")
#adding Duration of event
df_short$Duration <-df_short$halfWidth*2
#Replacing the labiling of "link" to process in pericyte structure
df_short$roiName <- str_replace_all(df_short$roiName, 'l1', 'plink')
#create a new variable for process and soma
df_short$ROI<- ifelse(grepl(pattern = "p", df_short$roiName, ignore.case = T),"Process","Soma")
#verifying categorical variables of the data frame
unique(df_short$drug)
unique(df_short$ROI)
unique(df_short$Animal)
df_short <- df_short[ ,c('date', 'Animal', 'Unique_Spot', 'Unique_ROIname', 'ROI', 'Condition', 'drug', 'amplitude', 'Duration','peakAUC')]
view(df_short)

#Remove negative amplitudes before calculating the mean
df_short<- df_short%>%
  filter(amplitude>0)

# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(df_short, c("date", 'Animal', "Unique_Spot", "Unique_ROIname", "ROI", "Condition", "drug"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(Duration,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude))


view(ROI.means)

#Converting to nominal factor the modified data frame
ROI.means$date = factor(ROI.means$date)
ROI.means$Animal = factor(ROI.means$Animal)
ROI.means$Unique_Spot= factor(ROI.means$Unique_Spot)
ROI.means$Unique_ROIname = factor(ROI.means$Unique_ROIname)
ROI.means$ROI = factor(ROI.means$ROI)
ROI.means$Condition = factor(ROI.means$Condition)
ROI.means$drug = factor(ROI.means$drug)

summary(ROI.means) 

########################################################################################################
#SECTION 2 OF THE CODE STATISTICS


#Analizyng capillary pericyte Calcium response to different drugs in NO STIM conditions

#FILTER DATA SET OF NO STIM. 
data.no.stim <- ROI.means%>%
  filter(Condition == "no stim")
view(data.no.stim)
summary(data.no.stim)

data.no.stim$date = factor(data.no.stim$date)
data.no.stim$Animal = factor(data.no.stim$Animal)
data.no.stim$Unique_Spot = factor(data.no.stim$Unique_Spot)
data.no.stim$Unique_ROIname = factor(data.no.stim$Unique_ROIname)
data.no.stim$ROI = factor(data.no.stim$ROI)
data.no.stim$drug  = factor(data.no.stim$drug)


#AMPLITUDE

library(plyr)
ddply(data.no.stim, ~ drug, function(data) summary(data$amp_mean))
ddply(data.no.stim, ~ drug, summarise, amp.mean=mean(amp_mean), amp.sd=sd(amp_mean))

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$amp_mean)
hist(data.no.stim[data.no.stim$drug == "peg",]$amp_mean)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$amp_mean)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$amp_mean)
boxplot(amp_mean ~ drug, data=data.no.stim, xlab="Drug", ylab="amplitude") # boxplots to see residuals

#Testing for normality in the residuals
library(plyr)
m = aov(amp_mean ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(amp_mean ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(amp_mean ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#See if data fits log normality 
# Kolmogorov-Smirnov test for log-normality is an option to use, but I just used it if justification is needed.


#MODIFICATION OF THE DATA WITH LOGNORMAL DISTRIBUTION
#creating log of amp_mean
data.no.stim$logamp_mean = log(data.no.stim$amp_mean) # log transform
View(data.no.stim) # verify

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$logamp_mean)
hist(data.no.stim[data.no.stim$drug == "peg",]$logamp_mean)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$logamp_mean)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$logamp_mean)
boxplot(logamp_mean ~ drug, data=data.no.stim, xlab="Drug", ylab="amplitude") # boxplots

#Testing for normality in the residuals with Log transformation
library(plyr)
m = aov(logamp_mean ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(logamp_mean ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(logamp_mean ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#DOING LINEAR MIXED MODELS 
# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls to compare each of the levels evenly to each other
#in pair wise comparisons.
contrasts(data.no.stim$date) <- "contr.sum"
contrasts(data.no.stim$Unique_Spot) <- "contr.sum"
contrasts(data.no.stim$Unique_ROIname) <- "contr.sum"
contrasts(data.no.stim$ROI)<- "contr.sum"
contrasts(data.no.stim$Condition) <- "contr.sum"
contrasts(data.no.stim$drug) <- "contr.sum"

# LMM with date, Animal, Spot and ROI as random effect for DRUG factor
m = lmer(logamp_mean ~ (drug) + (1|Animal) + (1|Unique_Spot:date) + (1|Unique_ROIname), data=data.no.stim)
Anova(m, type=2, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))








#FREQUENCY 


library(plyr)
ddply(data.no.stim, ~ drug, function(data) summary(data$nEvents))
ddply(data.no.stim, ~ drug, summarise, nEvents=mean(nEvents), nEvents.sd=sd(nEvents, na.rm=FALSE))
#The standard deviation of a length-one or zero-length vector is NA.

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$nEvents)
hist(data.no.stim[data.no.stim$drug == "peg",]$nEvents)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$nEvents)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$nEvents)
boxplot(nEvents ~ drug, data=data.no.stim, xlab="Drug", ylab="amplitude") # boxplots to see residuals

#Testing for normality in the residuals
library(plyr)
m = aov(nEvents ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(nEvents ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(nEvents ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#See if data fits log normality 
# Kolmogorov-Smirnov test for log-normality is an option to use, but I just used it if justification is needed.

#MODIFICATION OF THE DATA WITH LOGNORMAL DISTRIBUTION
#creating log of amp_mean
data.no.stim$lognEvents = log(data.no.stim$nEvents) # log transform
View(data.no.stim) # verify
summary(data.no.stim)

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$lognEvents)
hist(data.no.stim[data.no.stim$drug == "peg",]$lognEvents)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$lognEvents)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$lognEvents)
boxplot(lognEvents ~ drug, data=data.no.stim, xlab="Drug", ylab="Frequency") # boxplots to see residuals

#Testing for normality in the residuals with Log transformation
library(plyr)
m = aov(lognEvents ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(lognEvents ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(lognEvents ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#DOING LINEAR MIXED MODELS 
# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls to compare each of the levels evenly to each other
#in pair wise comparisons.
contrasts(data.no.stim$date) <- "contr.sum"
contrasts(data.no.stim$Unique_Spot) <- "contr.sum"
contrasts(data.no.stim$Unique_ROIname) <- "contr.sum"
contrasts(data.no.stim$ROI)<- "contr.sum"
contrasts(data.no.stim$Condition) <- "contr.sum"
contrasts(data.no.stim$drug) <- "contr.sum"

# LMM with date, Animal, Spot and ROI as random effect for DRUG factor
m = lmer(lognEvents ~ (drug)   + (1|Animal) + (1|Unique_Spot:date) + (1|Unique_ROIname), data=data.no.stim)
Anova(m, type=2, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))






#Peak Area Under the Curve 

library(plyr)
ddply(data.no.stim, ~ drug, function(data) summary(data$AUC_mean))
ddply(data.no.stim, ~ drug, summarise, AUC.mean=mean(AUC_mean), AUC.sd=sd(AUC_mean))

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$AUC_mean)
hist(data.no.stim[data.no.stim$drug == "peg",]$AUC_mean)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$AUC_mean)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$AUC_mean)
boxplot(AUC_mean ~ drug, data=data.no.stim, xlab="Drug", ylab="Auc") # boxplots to see residuals

#Testing for normality in the residuals
library(plyr)
m = aov(AUC_mean ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(AUC_mean ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(AUC_mean ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#See if data fits log normality 
# Kolmogorov-Smirnov test for log-normality is an option to use, but I just used it if justification is needed.


#MODIFICATION OF THE DATA WITH LOGNORMAL DISTRIBUTION
#creating log of amp_mean
data.no.stim$logAUC_mean = log(data.no.stim$AUC_mean) # log transform
View(data.no.stim) # verify

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$logAUC_mean)
hist(data.no.stim[data.no.stim$drug == "peg",]$logAUC_mean)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$logAUC_mean)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$logAUC_mean)
boxplot(logAUC_mean ~ drug, data=data.no.stim, xlab="Drug", ylab="logAUC_mean") # boxplots

#Testing for normality in the residuals with Log transformation
library(plyr)
m = aov(logAUC_mean ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(logAUC_mean ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(logAUC_mean ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#DOING LINEAR MIXED MODELS 
# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls to compare each of the levels evenly to each other
#in pair wise comparisons.
contrasts(data.no.stim$date) <- "contr.sum"
contrasts(data.no.stim$Unique_Spot) <- "contr.sum"
contrasts(data.no.stim$Unique_ROIname) <- "contr.sum"
contrasts(data.no.stim$ROI)<- "contr.sum"
contrasts(data.no.stim$Condition) <- "contr.sum"
contrasts(data.no.stim$drug) <- "contr.sum"

# LMM with date, Animal, Spot and ROI as random effect for DRUG factor
m = lmer(logAUC_mean ~ (drug) + (1|Animal) + (1|Unique_Spot:date) + (1|Unique_ROIname), data=data.no.stim)
Anova(m, type=2, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))








#PEAK WIDTH 


library(plyr)
ddply(data.no.stim, ~ drug, function(data) summary(data$dur_mean))
ddply(data.no.stim, ~ drug, summarise, dur.mean=mean(dur_mean), dur.sd=sd(dur_mean))

# histograms for two factors
hist(data.no.stim[data.no.stim$drug == "baseline",]$dur_mean)
hist(data.no.stim[data.no.stim$drug == "peg",]$dur_mean)
hist(data.no.stim[data.no.stim$drug == "nimodipine",]$dur_mean)
hist(data.no.stim[data.no.stim$drug == "pyr3",]$dur_mean)
boxplot(AUC_mean ~ drug, data=data.no.stim, xlab="Drug", ylab="Peak Duration") # boxplots to see residuals

#Testing for normality in the residuals
library(plyr)
m = aov(dur_mean ~ drug, data=data.no.stim) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(dur_mean ~ drug, data=data.no.stim, center=mean) # Levene's test
leveneTest(dur_mean ~ drug, data=data.no.stim, center=median) # Brown-Forsythe test

#DATA IS NORMALLY DISTRIBUTED. NO NEED OF LOG TRANSFORMATION.

#DOING LINEAR MIXED MODELS 
# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls to compare each of the levels evenly to each other
#in pair wise comparisons.
contrasts(data.no.stim$date) <- "contr.sum"
contrasts(data.no.stim$Unique_Spot) <- "contr.sum"
contrasts(data.no.stim$Unique_ROIname) <- "contr.sum"
contrasts(data.no.stim$ROI)<- "contr.sum"
contrasts(data.no.stim$Condition) <- "contr.sum"
contrasts(data.no.stim$drug) <- "contr.sum"

# LMM with date, Animal, Spot and ROI as random effect for DRUG factor
m = lmer(dur_mean ~ (drug) + (1|Animal) + (1|Unique_Spot:date) + (1|Unique_ROIname), data=data.no.stim)
Anova(m, type=2, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))


##########################################################################################


#Analizyng capillary pericyte calcium signaling in different structures.




