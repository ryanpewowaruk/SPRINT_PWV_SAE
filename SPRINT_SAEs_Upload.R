## Set up - Read in and clean dataa

library(haven)
library(ggplot2)
library(Hmisc)
library(MASS)
library(splines)
library(contrast)
library(rms)


setwd("C:/Users/pewow/Documents/SPRINT/SAEs")

# Read in stiffness data
PWV <- read.csv("SPRINT_Supiano_PWV_Mechanisms.csv")

# read in SAE data
SAEs <- read_sas("safety_events_v2.sas7bdat")
SAEs$maskid <- SAEs$MASKID

# Read in CVD events data for final follow-up time for each participant
CVD <- read.csv("outcomes.csv")

# Read in frailty data
frailty <- read.csv("sprint_frailty_classic.csv")


# Read in demographics and baseline descriptors
DEMO_raw <-  read_sas("baseline.sas7bdat")

# Make dataframe with necesseary demographic and co-variate data
DEMO <- data.frame(
  maskid = DEMO_raw$maskid,
  age = DEMO_raw$age,
  female = DEMO_raw$female,
  race = DEMO_raw$race4,
  bmi = DEMO_raw$BMI,
  smoke = factor(DEMO_raw$smoke_3cat),
  cvd = DEMO_raw$sub_cvd,
  ckd = DEMO_raw$sub_ckd,
  fbg = DEMO_raw$GLUR,
  n_BP_med = DEMO_raw$n_agents,
  fr_risk = DEMO_raw$fr_risk10yrs,
  SBP = DEMO_raw$sbp,
  DBP = DEMO_raw$dbp,
  scr = DEMO_raw$screat)

DEMO$MAP <- as.numeric( with(DEMO, (2*DBP + SBP)/3) ) # Calculate MAP

DEMO$race2 <- DEMO$race == "BLACK" # Make 2-category race variable Black and non-Black

PWV$year<- as.numeric( sub('.', '', PWV$Visit) ) # Remove 'T' from visit string

year0 <- PWV$year == 0
# Sort PWV data into baseline (year 0) and follow up (years 1-3)
baseline <- subset(PWV, year0)

# Merge datasets into one
baseline <- merge(baseline, DEMO, by="maskid") # merge baseline 

data <- merge(baseline, CVD,  by="maskid")

data <- merge(data, frailty, by="maskid")

# Calculate follow-up time for each participant
data$t_fu <- with(data, pmax(t_cvddeath, t_nonmiacs)) # no participants had CVD death and ACS so this is convenient way to calculate max FU time

data <- subset(data, data$t_fu>0) # exclude participants w/o follow-up time

data$randAssign <- factor(data$randAssign) # convert treatment assignment from numberic into factor

# initalize variables for calculating SAEs
data$n_SAE <- 0 # Serious adverse events
data$n_ERV <- 0 # Emergency room visits
data$n_LAB <- 0 # Abnormal labs
data$n_OH <- 0 # Orthostatic hypotension

data$SAE_hypo <- 0 # hypotension SAE
data$SAE_synco <- 0 # syncope SAE
data$SAE_brady <- 0 # bradycardia SAE
data$SAE_electro <- 0 # electrolyte SAE
data$SAE_fall <- 0 # injurious fall SAE
data$SAE_aki <- 0 # acute kidney injury SAE

data$OH_dizzy <- 0 # reported dizziness during OH


# Loop to calculate number of SAEs for each participant
for (i in 1:length(data$maskid) )  {
  
  id = data$maskid[i]
  
  temp_SAE <- subset(SAEs, SAEs$maskid == id ) 
  temp_SAE <- subset(temp_SAE, temp_SAE$EVENTDAYS <= data$t_fu[i] )

  print( c(id, i,  max(temp_SAE$OH_SYMPTOMS, na.rm=TRUE) ) )
  
  n_SAE <- sum( as.numeric( temp_SAE$TYPE == "SAE" ) )
  data$n_SAE[i] <- n_SAE
  data$SAE[i] <- as.numeric(n_SAE > 0) 
  
  data$SAE_hypo[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_HYPOTEN  ), na.rm=TRUE )} else {0}
  data$SAE_synco[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_SYNCOPE  ), na.rm=TRUE )} else {0}
  data$SAE_brady[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_BRADY  ), na.rm=TRUE )} else {0}
  data$SAE_electro[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_ELECTRO  ), na.rm=TRUE )} else {0}
  data$SAE_fall[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_INJFALL  ), na.rm=TRUE )} else {0}
  data$SAE_aki[i] <- if(n_SAE>0) {max(as.numeric( temp_SAE$MC_AKI  ), na.rm=TRUE )} else {0}
  
  
  n_ERV <- sum( as.numeric( temp_SAE$TYPE == "ERV" ) )
  data$n_ERV[i] <- n_ERV
  data$ERV[i] <- as.numeric(n_ERV > 0) 
  
  n_OH <- sum( as.numeric( temp_SAE$TYPE == "OH" ) )
  data$n_OH[i] <- n_OH
  data$OH[i] <- as.numeric(n_OH > 0) 
  
  data$OH_dizzy[i] <- if(n_OH>0) {max(as.numeric( temp_SAE$OH_SYMPTOMS  ), na.rm=TRUE ) } else {0}
  
  n_LAB <- sum( as.numeric( temp_SAE$TYPE == "LAB" ) )
  data$n_LAB[i] <- n_LAB
  data$LAB[i] <- as.numeric(n_LAB > 0) 
  
}


# Clean Dataset
MISSING_covar <- is.na(data$age) | # remove participants with missing co-variates
  is.na(data$female) |
  is.na(data$race) |
  is.na(data$MAP) |
  is.na(data$n_BP_med) |
  is.na(data$smoke) |
  is.na(data$scr) 
sum(MISSING_covar)

data <- subset(data, subset = !MISSING_covar)


# Define datadist so RMS functions work
dd <- datadist(baseline)
options(datadist='dd')

## Make Table 1

Table1 <- data.frame(
  Tx = data$randAssign,
  age = data$age,
  female = data$female,
  race = data$race,
  bmi = data$bmi,
  smoke = data$smoke,
  cvd = data$cvd,
  ckd = data$ckd,
  n_BP_med = data$n_BP_med,
  SBP = data$SBP,
  DBP = data$DBP,
  scr = data$scr,
  fi_ctns = data$frailty_ctns,
  fi_cat = data$frailty_catg,
  PWV = data$PWV,
  PWV_Total = data$PWV_Total,
  PWV_struct = data$PWV_struct,
  PWV_LD = data$PWV_LD,
  t_fu = data$t_fu
  )

# report out table 1 for each treatment arm
describe(subset(Table1, Table1$Tx==0) )
describe(subset(Table1, Table1$Tx==1) )

## Histograms of Event Counts

hlm <- coord_cartesian(ylim = c(0, 652), xlim=c(0,11)) # set axis coordinates for histogram

# histogram for SAEs
ggplot(data, aes(x=n_SAE)) + geom_histogram( binwidth=1) +  
  labs(title = "SAE", x="Number of Events") + hlm + 
  theme_bw()

# histogram for ER Visits
ggplot(data, aes(x=n_ERV)) + geom_histogram( binwidth=1) +  
  labs(title = "ER Visits", x="Number of Events") +hlm + 
  theme_bw()

# histogram for abnormal labs
ggplot(data, aes(x=n_LAB)) + geom_histogram( binwidth=1) +  
  labs(title = "Abnormal Labs", x="Number of Events") + hlm + 
  theme_bw()

# histogram for orthostatic hypotension
ggplot(data, aes(x=n_OH)) + geom_histogram( binwidth=1) +  
  labs(title = "Orthostatic Hypotension", x="Number of Events") + hlm + 
  theme_bw()


## Convert time into Visits becayse OH and labs only assessed at regularly scheduled visits

for (i in 1:length(data$maskid) )  {
  
  t <- data$t_fu[i]
  
  # OH (1-month, 6-months, 1-year, annually after)
  
  if (t<365/2 - 60) {data$visit_OH[i] <- 1}
  else if (t<365 - 60 ) {data$visit_OH[i] <- 2}
  else if (t<365*2 - 60 ) {data$visit_OH[i] <- 3}
  else if (t<365*3 - 60 ) {data$visit_OH[i] <- 4}
  else {data$visit_OH[i] <- 5}
  
  # Labs (1-month, 6-months, 1-year, annually after)
  
  if (t<365/4 - 30) {data$visit_LAB[i] <- 1}
  else if (t<365/2 - 30 ) {data$visit_LAB[i] <- 2}
  else if (t<365 - 60 ) {data$visit_LAB[i] <- 3}
  else if (t<365*1.5 - 60 ) {data$visit_LAB[i] <- 4}
  else if (t<365*2 - 60 ) {data$visit_LAB[i] <- 5}
  else if (t<365*2.5 - 60 ) {data$visit_LAB[i] <- 6}
  else if (t<365*3 - 60 ) {data$visit_LAB[i] <- 7}
  else if (t<365*3.5 - 60 ) {data$visit_LAB[i] <- 8}
  else {data$visit_LAB[i] <- 9}
  
}

describe(data$visit_OH)
describe(data$visit_LAB)

## Make Table 2 - Number of Adverse Events

Table2 <- data.frame(
  Tx = data$randAssign,
  SAE = data$SAE,
  SAE_hypo = data$SAE_hypo,
  SAE_synco = data$SAE_synco,
  SAE_brady = data$SAE_brady,
  SAE_electro = data$SAE_electro,
  SAE_fall = data$SAE_fall,
  SAE_aki = data$SAE_aki,
  LAB = data$LAB,
  OH = data$OH,
  OH_dizzy = data$OH_dizzy
)

# SAEs that don't fall into a pre-defined category
Table2$SAE_undefined <- with(Table2, SAE & (SAE_hypo + SAE_synco + SAE_brady + SAE_electro + SAE_fall + SAE_aki) == 0 )

# Report out SAEs for each treatment group
describe(subset(Table2, Table2$Tx==0) )
describe(subset(Table2, Table2$Tx==1) )

# Calculate p-value for treatment group difference for SAEs
f_SAE_justTx <- glm.nb(n_SAE ~ randAssign, 
                       offset(log(t_fu)), data=data)
summary(f_SAE_justTx)

# Calculate p-value for treatment group difference for OH
f_OH_justTx <- glm.nb(n_SAE ~ randAssign, 
                       offset(log(visit_OH)), data=data)
summary(f_OH_justTx)

## Spike Histograms

# Make spike histogram to add to plots for Total PWV
pwv_tot <- data.frame( spikecomp(baseline$PWV_Total ) )
pwv_tot$y <- pwv_tot$y.Freq

# Make spike histogram to add to plots for structural PWV
pwv_struct <- data.frame( spikecomp(baseline$PWV_struct ) )
pwv_struct$y <- pwv_struct$y.Freq

# Make spike histogram to add to plots for load-dependent PWV
pwv_ld <- data.frame( spikecomp(baseline$PWV_LD ) )
pwv_ld$y <- pwv_ld$y.Freq

# Make spike histogram to add to plots for mean arterial pressure (MAP)
map <- data.frame( spikecomp(baseline$MAP ) )
map$y <- map$y.Freq


## Define knot locations for restricted cubic splines (knots set at 10th and 90th percentile)
k_age <- quantile(data$age, c(0.1, 0.9))
k_MAP <- quantile(data$MAP, c(0.1, 0.9))
k_n_BP_med <- quantile(data$n_BP_med, c(0.1, 0.9))
k_scr <- quantile(data$scr, c(0.1, 0.9), na.rm=TRUE)
k_PWV_tot <- quantile(data$PWV_Total, c(0.1, 0.9))
k_PWV_struct <- quantile(data$PWV_struct, c(0.1, 0.9))
k_PWV_LD <- quantile(data$PWV_LD, c(0.1, 0.9))
k_FI <- quantile(data$frailty_ctns, c(0.1, 0.9), na.rm=TRUE)
 
## Make Prediction Data Frame for plotting splines of PWV
data_pred <- data.frame(
  randAssign = factor(rep(0, 101) ),
  age = rep(median(data$age), 101),
  female = rep(median(data$female), 101),
  race = rep("WHITE", 101),
  smoke = factor(rep(1, 101) ),
  MAP = rep(median(data$MAP), 101),
  n_BP_med = rep(median(data$n_BP_med), 101),
  scr = rep(median(data$scr, na.rm=TRUE), 101),
  PWV_Total = seq(from = quantile(data$PWV_Total, 0.025 ), to = quantile(data$PWV_Total, 0.975 ), length.out=101) ,
  PWV_struct = seq(from = quantile(data$PWV_struct, 0.025 ), to = quantile(data$PWV_struct, 0.975 ), length.out=101) ,
  PWV_LD = seq(from = quantile(data$PWV_LD, 0.025 ), to = quantile(data$PWV_LD, 0.975 ), length.out=101),
  t_fu = 10*rep(365, 101)
)

# Make separate data frame for plotting MAP spines
data_pred_MAP <- data.frame(
  randAssign = factor(rep(0, 101) ),
  age = rep(median(data$age), 101),
  female = rep(median(data$female), 101),
  race = rep("WHITE", 101),
  smoke = factor(rep(1, 101) ),
  MAP = seq(from = quantile(data$MAP, 0.025 ), to = quantile(data$MAP, 0.975 ), length.out=101),
  n_BP_med = rep(median(data$n_BP_med), 101),
  scr = rep(median(data$scr, na.rm=TRUE), 101),
  PWV_Total = seq(from = quantile(data$PWV_Total, 0.025 ), to = quantile(data$PWV_Total, 0.975 ), length.out=101) ,
  PWV_struct = seq(from = quantile(data$PWV_struct, 0.025 ), to = quantile(data$PWV_struct, 0.975 ), length.out=101) ,
  PWV_LD = seq(from = quantile(data$PWV_LD, 0.025 ), to = quantile(data$PWV_LD, 0.975 ), length.out=101),
  t_fu = 10*rep(365, 101)
  # FI <- quantile(data$age, c(0.1, 0.5, 0.9))
  
)

## SAE Negative Binomial Regressions

ylm <- coord_cartesian(ylim = c(0, 1.2)) # set plot axis coordinates

# regression w/o MAP
f_SAE_noMAP <- glm.nb(n_SAE ~ randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                  ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + 
                  ns(scr, df=2, Boundary.knots = k_scr), 
                offset(log(t_fu)), data=data)

# regression w MAP
f_SAE <- glm.nb(n_SAE ~ randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                      ns(MAP, df=2, Boundary.knots = k_MAP) + 
                  ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + 
                  ns(scr, df=2, Boundary.knots = k_scr), 
                    offset(log(t_fu)), data=data)
anova(f_SAE_noMAP, f_SAE) # calculate p-value for adding MAP

# Make splines for plotting number of events
pred_SAE_map <- cbind(data_pred_MAP, predict(f_SAE, data_pred_MAP, type = "link", se.fit=TRUE))
pred_SAE_map <- within(pred_SAE_map, {
  n_SAE <- exp(fit) # convert log(events) to events
  LL <- exp(fit - 1.96 * se.fit) # calculate lower 95% CI
  UL <- exp(fit + 1.96 * se.fit) # calculate upper 95% CI
})

# Plot splines
ggplot(pred_SAE_map, aes(x=MAP, y=n_SAE)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0), linewidth=2) # last line is to add histogram


# regression w total PWV
f_SAE_tot <- glm.nb(n_SAE ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                      ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                    offset(log(t_fu)), data=data)

anova(f_SAE, f_SAE_tot)


pred_SAE_tot <- cbind(data_pred, predict(f_SAE_tot, data_pred, type = "link", se.fit=TRUE))
pred_SAE_tot <- within(pred_SAE_tot, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_tot, aes(x=PWV_Total, y=n_SAE)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)

contrast(f_SAE_tot, 
         list(randAssign=0, age=median(data$age), female=median(data$female),
              race=median(data$race), smoke=1, MAP=median(data$MAP),
              n_BP_med = median(data$n_BP_med), scr = median(data$scr), 
              PWV_Total = quantile(data$PWV_Total, 0.25) ),
         
         list(randAssign=0, age=median(data$age), female=median(data$female),
              race=median(data$race), smoke=1, MAP=median(data$MAP),
              n_BP_med = median(data$n_BP_med), scr = median(data$scr), 
              PWV_Total = quantile(data$PWV_Total, 0.75) )
         )

# regression w structural PWV
f_SAE_struct <- glm.nb(n_SAE ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)


anova(f_SAE, f_SAE_struct)

pred_SAE_struct <- cbind(data_pred, predict(f_SAE_struct, data_pred, type = "link", se.fit=TRUE))
pred_SAE_struct <- within(pred_SAE_struct, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_struct, aes(x=PWV_struct, y=n_SAE)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)

# regression w load-dependent PWV
f_SAE_LD <- glm.nb(n_SAE ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                     ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                   offset(log(t_fu)), data=data)

anova(f_SAE, f_SAE_LD)

pred_SAE_LD <- cbind(data_pred, predict(f_SAE_LD, data_pred, type = "link", se.fit=TRUE))
pred_SAE_LD <- within(pred_SAE_LD, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_LD, aes(x=PWV_LD, y=n_SAE)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)

## Orthostatic Hypotension Negative Binomial Regressions

# regression w/o MAP
f_OH_noMAP <- glm.nb(n_OH ~ randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                 ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + 
                 ns(scr, df=2, Boundary.knots = k_scr),
               offset(log(visit_OH)), data=data)

# regression w MAP
f_OH <- glm.nb(n_OH ~ randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                 ns(MAP, df=2, Boundary.knots = k_MAP) + 
                 ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + 
                 ns(scr, df=2, Boundary.knots = k_scr),
                offset(log(visit_OH)), data=data)

anova(f_OH_noMAP, f_OH)

pred_OH_map <- cbind(data_pred_MAP, predict(f_OH, data_pred_MAP, type = "link", se.fit=TRUE))
pred_OH_map <- within(pred_OH_map, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_map, aes(x=MAP, y=n_OH)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)


# regression w total PWV
f_OH_tot <- glm.nb( n_OH ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
  ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                    offset(log(visit_OH)), data=data)

anova(f_OH, f_OH_tot)

pred_OH_tot <- cbind(data_pred, predict(f_OH_tot, data_pred, type = "link", se.fit=TRUE))
pred_OH_tot <- within(pred_OH_tot, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_tot, aes(x=PWV_Total, y=n_OH)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)

# regression w structural PWV
f_OH_struct <- glm.nb(n_OH ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                        ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                      offset(log(visit_OH)), data=data)

anova(f_OH, f_OH_struct)

pred_OH_struct <- cbind(data_pred, predict(f_OH_struct, data_pred, type = "link", se.fit=TRUE))
pred_OH_struct <- within(pred_OH_struct, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_struct, aes(x=PWV_struct, y=n_OH)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)

# regression w load-dependent PWV
f_OH_LD <- glm.nb(n_OH ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) + randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                    ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                  offset(log(visit_OH)), data=data)

anova(f_OH, f_OH_LD)

pred_OH_LD <- cbind(data_pred, predict(f_OH_LD, data_pred, type = "link", se.fit=TRUE))
pred_OH_LD <- within(pred_OH_LD, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_LD, aes(x=PWV_LD, y=n_OH)) +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
  geom_line( size = 1) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0), linewidth=2)


## SAEs - Assess for interaction with Treatment Group

data_pred_Tx <- data.frame(
  randAssign = factor(c(rep(0,101), rep(1,101) ) ),
  age = rep(median(data$age), 2*101),
  female = rep(median(data$female), 2*101),
  race = rep("WHITE", 2*101),
  smoke = factor(rep(1, 2*101) ),
  MAP = rep(median(data$MAP), 2*101),
  n_BP_med = rep(median(data$n_BP_med), 2*101),
  scr = rep(median(data$scr, na.rm=TRUE), 2*101),
  PWV_Total = rep(seq(from = quantile(data$PWV_Total, 0.025 ), to = quantile(data$PWV_Total, 0.975 ), length.out=101) , 2),
  PWV_struct = rep(seq(from = quantile(data$PWV_struct, 0.025 ), to = quantile(data$PWV_struct, 0.975 ), length.out=101), 2) ,
  PWV_LD = rep(seq(from = quantile(data$PWV_LD, 0.025 ), to = quantile(data$PWV_LD, 0.975 ), length.out=101), 2) ,
  # FI <- quantile(data$age, c(0.1, 0.5, 0.9))
  
)

data_pred_Tx_MAP <- data_pred_Tx
data_pred_Tx_MAP$MAP <- rep(seq(from = quantile(data$MAP, 0.025 ), to = quantile(data$MAP, 0.975 ), length.out=101) , 2)

# regression w MAP*treatment interaction
f_SAE_MAP_Tx <- glm.nb(n_SAE ~ ns(MAP, df=2, Boundary.knots=k_MAP) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)

anova(f_SAE, f_SAE_MAP_Tx)

pred_SAE_MAP_Tx <- cbind(data_pred_Tx_MAP, predict(f_SAE_MAP_Tx, data_pred_Tx_MAP, type = "link", se.fit=TRUE))
pred_SAE_MAP_Tx <- within(pred_SAE_MAP_Tx, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_MAP_Tx, aes(x=MAP, y=n_SAE, group=randAssign) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)


f_SAE_tot_Tx <- glm.nb(n_SAE ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                      ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                    offset(log(t_fu)), data=data)

anova(f_SAE_tot, f_SAE_tot_Tx)

pred_SAE_tot_Tx <- cbind(data_pred_Tx, predict(f_SAE_tot_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_SAE_tot_Tx <- within(pred_SAE_tot_Tx, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_tot_Tx, aes(x=PWV_Total, y=n_SAE, group=randAssign) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)


f_SAE_struct_Tx <- glm.nb(n_SAE ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)


anova(f_SAE_struct, f_SAE_struct_Tx)

pred_SAE_struct_Tx <- cbind(data_pred_Tx, predict(f_SAE_struct_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_SAE_struct_Tx <- within(pred_SAE_struct_Tx, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_struct_Tx, aes(x=PWV_struct, y=n_SAE, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)



f_SAE_LD_Tx <- glm.nb(n_SAE ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                     ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                   offset(log(t_fu)), data=data)

anova(f_SAE_LD, f_SAE_LD_Tx)

pred_SAE_LD_Tx <- cbind(data_pred_Tx, predict(f_SAE_LD_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_SAE_LD_Tx <- within(pred_SAE_LD_Tx, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_LD_Tx, aes(x=PWV_LD, y=n_SAE, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)




## OH - Assess for interaction with Treatment Group

f_OH_MAP_Tx <- glm.nb(n_OH ~ ns(MAP, df=2, Boundary.knots=k_MAP) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(visit_OH)), data=data)

anova(f_OH, f_OH_MAP_Tx)

pred_OH_MAP_Tx <- cbind(data_pred_Tx_MAP, predict(f_OH_MAP_Tx, data_pred_Tx_MAP, type = "link", se.fit=TRUE))
pred_OH_MAP_Tx <- within(pred_OH_MAP_Tx, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_MAP_Tx, aes(x=MAP, y=n_OH, group=randAssign) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)



f_OH_tot_Tx <- glm.nb(n_OH ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(visit_OH)), data=data)

anova(f_OH_tot, f_OH_tot_Tx)

pred_OH_tot_Tx <- cbind(data_pred_Tx, predict(f_OH_tot_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_OH_tot_Tx <- within(pred_OH_tot_Tx, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_tot_Tx, aes(x=PWV_Total, y=n_OH, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)


f_OH_struct_Tx <- glm.nb(n_OH ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                        ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                      offset(log(visit_OH)), data=data)

anova(f_OH_struct, f_OH_struct_Tx)

pred_OH_struct_Tx <- cbind(data_pred_Tx, predict(f_OH_struct_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_OH_struct_Tx <- within(pred_OH_struct_Tx, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_struct_Tx, aes(x=PWV_struct, y=n_OH, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)


f_OH_LD_Tx <- glm.nb(n_OH ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) * randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                           ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                         offset(log(visit_OH)), data=data)

anova(f_OH_LD, f_OH_LD_Tx)

pred_OH_LD_Tx <- cbind(data_pred_Tx, predict(f_OH_LD_Tx, data_pred_Tx, type = "link", se.fit=TRUE))
pred_OH_LD_Tx <- within(pred_OH_LD_Tx, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_LD_Tx, aes(x=PWV_LD, y=n_OH, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=randAssign), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0, group=1), linewidth=2)


ggplot(pred_OH_LD_Tx, aes(x=PWV_LD, y=n_OH, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25) +
  geom_line( size = 1, aes(color=randAssign)) + coord_cartesian(ylim=c(0,1.25), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) )


## SAEs - Assess for interaction with Frailty

describe(subset(Table1$fi_ctns, Table1$fi_cat=="Fit") )
describe(subset(Table1$fi_ctns, Table1$fi_cat=="Pre-frail") )
describe(subset(Table1$fi_ctns, Table1$fi_cat=="Frail") )


q_FI <- c(0.07, 0.15, 0.25)

data_pred_FI <- data.frame(
  randAssign = factor(rep(0,3*101) ),
  age = rep(median(data$age), 3*101),
  female = rep(median(data$female), 3*101),
  race = rep("WHITE", 3*101),
  smoke = factor(rep(1, 3*101) ),
  MAP = rep(median(data$MAP), 3*101),
  n_BP_med = rep(median(data$n_BP_med), 3*101),
  scr = rep(median(data$scr, na.rm=TRUE), 3*101),
  PWV_Total = rep(seq(from = quantile(data$PWV_Total, 0.025 ), to = quantile(data$PWV_Total, 0.975 ), length.out=101) , 3),
  PWV_struct = rep(seq(from = quantile(data$PWV_struct, 0.025 ), to = quantile(data$PWV_struct, 0.975 ), length.out=101), 3) ,
  PWV_LD = rep(seq(from = quantile(data$PWV_LD, 0.025 ), to = quantile(data$PWV_LD, 0.975 ), length.out=101), 3),
  frailty_ctns <- c( rep(q_FI[1], 101), rep(q_FI[2], 101), rep(q_FI[3], 101) ),
  FI_cat <- c( rep("Fit", 101), rep("Pre-Frail", 101), rep("Frail", 101)) 
  
)


data_pred_FI_MAP <- data_pred_FI
data_pred_FI_MAP$MAP <- rep(seq(from = quantile(data$MAP, 0.025 ), to = quantile(data$MAP, 0.975 ), length.out=101) , 3)

f_SAE_MAP_FI <- glm.nb(n_SAE ~ ns(MAP, df=2, Boundary.knots=k_MAP) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)

f_SAE_MAP_FI_int <- glm.nb(n_SAE ~ ns(MAP, df=2, Boundary.knots=k_MAP) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                             ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                           offset(log(t_fu)), data=data)

anova(f_SAE_MAP_FI, f_SAE_MAP_FI_int)

pred_SAE_MAP_FI <- cbind(data_pred_FI_MAP, predict(f_SAE_MAP_FI_int, data_pred_FI_MAP, type = "link", se.fit=TRUE))
pred_SAE_MAP_FI <- within(pred_SAE_MAP_FI, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_MAP_FI, aes(x=MAP, y=n_SAE, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")



f_SAE_tot_FI <- glm.nb(n_SAE ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)

f_SAE_tot_FI_int <- glm.nb(n_SAE ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                            ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                          offset(log(t_fu)), data=data)

anova(f_SAE_tot_FI, f_SAE_tot_FI_int)

pred_SAE_tot_FI <- cbind(data_pred_FI, predict(f_SAE_tot_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_SAE_tot_FI <- within(pred_SAE_tot_FI, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_tot_FI, aes(x=PWV_Total, y=n_SAE, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")



f_SAE_struct_FI <- glm.nb(n_SAE ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(t_fu)), data=data)

f_SAE_struct_FI_int <- glm.nb(n_SAE ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                             ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                           offset(log(t_fu)), data=data)

anova(f_SAE_struct_FI, f_SAE_struct_FI_int)

pred_SAE_struct_FI <- cbind(data_pred_FI, predict(f_SAE_struct_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_SAE_struct_FI <- within(pred_SAE_struct_FI, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(pred_SAE_struct_FI, aes(x=PWV_struct, y=n_SAE, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")


f_SAE_LD_FI <- glm.nb(n_SAE ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                            ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                          offset(log(t_fu)), data=data)

f_SAE_LD_FI_int <- glm.nb(n_SAE ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                                ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                              offset(log(t_fu)), data=data)

anova(f_SAE_LD_FI, f_SAE_LD_FI_int)

pred_SAE_LD_FI <- cbind(data_pred_FI, predict(f_SAE_LD_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_SAE_LD_FI <- within(pred_SAE_LD_FI, {
  n_SAE <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_SAE_LD_FI, aes(x=PWV_LD, y=n_SAE, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")


## OH - Assess for interaction with Frailty

f_OH_MAP_FI <- glm.nb(n_OH ~ ns(MAP, df=2, Boundary.knots=k_MAP) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(visit_OH)), data=data)

f_OH_MAP_FI_int <- glm.nb(n_OH ~ ns(MAP, df=2, Boundary.knots=k_MAP) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                             ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                           offset(log(visit_OH)), data=data)

anova(f_OH_MAP_FI, f_OH_MAP_FI_int)

pred_OH_MAP_FI <- cbind(data_pred_FI_MAP, predict(f_OH_MAP_FI_int, data_pred_FI_MAP, type = "link", se.fit=TRUE))
pred_OH_MAP_FI <- within(pred_OH_MAP_FI, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_MAP_FI, aes(x=MAP, y=n_OH, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$MAP, c(0.025, 0.975) ) ) +
  labs(x = "MAP (mmHg)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=map, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")




f_OH_tot_FI <- glm.nb(n_OH ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                         ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                       offset(log(visit_OH)), data=data)

f_OH_tot_FI_int <- glm.nb(n_OH ~ ns(PWV_Total, df=2, Boundary.knots=k_PWV_tot) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                             ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                           offset(log(visit_OH)), data=data)

anova(f_OH_tot_FI, f_OH_tot_FI_int)

pred_OH_tot_FI <- cbind(data_pred_FI, predict(f_OH_tot_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_OH_tot_FI <- within(pred_OH_tot_FI, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_tot_FI, aes(x=PWV_Total, y=n_OH, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_Total, c(0.025, 0.975) ) ) +
  labs(x = "Total PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_tot, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")



f_OH_struct_FI <- glm.nb(n_OH ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                            ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                          offset(log(visit_OH)), data=data)

f_OH_struct_FI_int <- glm.nb(n_OH ~ ns(PWV_struct, df=2, Boundary.knots=k_PWV_struct) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                                ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                              offset(log(visit_OH)), data=data)

anova(f_OH_struct_FI, f_OH_struct_FI_int)

pred_OH_struct_FI <- cbind(data_pred_FI, predict(f_OH_struct_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_OH_struct_FI <- within(pred_OH_struct_FI, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(pred_OH_struct_FI, aes(x=PWV_struct, y=n_OH, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_struct, c(0.025, 0.975) ) ) +
  labs(x = "Structural PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_struct, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")


f_OH_LD_FI <- glm.nb(n_OH ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) + ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                        ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                      offset(log(visit_OH)), data=data)

f_OH_LD_FI_int <- glm.nb(n_OH ~ ns(PWV_LD, df=2, Boundary.knots=k_PWV_LD) * ns(frailty_ctns, df=2, Boundary.knots=k_FI) +  randAssign + ns(age, df=2, Boundary.knots = k_age) + female + race + 
                            ns(MAP, df=2, Boundary.knots = k_MAP) + ns(n_BP_med, df=2, Boundary.knots = k_n_BP_med) +  smoke + ns(scr, df=2, Boundary.knots = k_scr), 
                          offset(log(visit_OH)), data=data)

anova(f_OH_LD_FI, f_OH_LD_FI_int)

pred_OH_LD_FI <- cbind(data_pred_FI, predict(f_OH_LD_FI_int, data_pred_FI, type = "link", se.fit=TRUE))
pred_OH_LD_FI <- within(pred_OH_LD_FI, {
  n_OH <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(pred_OH_LD_FI, aes(x=PWV_LD, y=n_OH, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25, show.legend=FALSE) +
  geom_line( size = 1, aes(color=FI_cat), show.legend=FALSE) + coord_cartesian(ylim=c(0,1.2), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  geom_segment(data=pwv_ld, aes(x, y/10 +0, xend=x, yend=0, group="FIT"), linewidth=2) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")

## Get Legend to add to figures
ggplot(pred_OH_LD_FI, aes(x=PWV_LD, y=n_OH, group=FI_cat) ) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=FI_cat), alpha = .25) +
  geom_line( size = 1, aes(color=FI_cat)) + coord_cartesian(ylim=c(0,1), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of OH") +
  theme_bw() +
  theme( text = element_text(size = 18) ) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")

ggplot(pred_SAE_LD_Tx, aes(x=PWV_LD, y=n_SAE, group=randAssign)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill=randAssign), alpha = .25) +
  geom_line( size = 1, aes(color=randAssign)) + coord_cartesian(ylim=c(0,1), xlim=quantile(data$PWV_LD, c(0.025, 0.975) ) ) +
  labs(x = "Load-Dep. PWV (m/s)", y = "Number of SAEs") +
  theme_bw() +
  theme( text = element_text(size = 18) )

## Scratch Paper


data_ERV <- subset(SAEs, SAEs$TYPE=="ERV")
data_LAB <- subset(SAEs, SAEs$TYPE=="LAB")
data_OH <- subset(SAEs, SAEs$TYPE=="OH")
data_SAE <- subset(SAEs, SAEs$TYPE=="SAE")

data$serious <- with(data, SAE_HOSPITAL + SAE_DISABLED + SAE_FATALITY + SAE_LFTHREAT + SAE_IMPMEDEVENT )

data_serious <- subset(data, data$serious>0)

data_related <- subset(data, !(data$RELATED == "" ) )
data_related$RELATED <- as.numeric(data_related$RELATED)
data_related <- subset(data_related, data_related$RELATED == 1 ) 

data_related <- rbind(data_related, data_LAB)

data_related <- rbind(data_related, data_OH)


describe(data_ERV)
describe(data_LAB)
describe(data_OH)
describe(data_SAE)
describe(data_related)
describe(data_serious)


data_SAE$rank <- rank(data_SAE$EVENTDAYS, ties.method = "first")

ggplot(data_SAE, aes(x=EVENTDAYS / 365, y=rank)) + geom_line(size = 1) +
  labs(title = "SAEs", x = "Years", y = "Number of Events") +
  theme_bw() + geom_smooth(method = "lm", color = "pink")

data_OH$rank <- rank(data_OH$EVENTDAYS, ties.method = "first")

ggplot(data_OH, aes(x=EVENTDAYS / 365, y=rank)) + geom_line(size = 1) +
  labs(title = "Orthostatic Hypotension", x = "Years", y = "Number of Events") +
  theme_bw() + geom_smooth(method = "lm", color = "pink")

for (i in 1:length(data_OH$maskid) )  {
  t <- data_OH$EVENTDAYS[i] /365
  
if (t<=1/12) {data_OH$visit[i] <- t*12 }
  else if (t <= 1/2) { data_OH$visit[i] <- (t-1/12 ) *(12/5) +1 } 
  else if ( t<= 1) { data_OH$visit[i] <- (t-1/2 ) *2 +2}
  else {data_OH$visit[i] <- t+2}

}

ggplot(data_OH, aes(x=visit, y=rank)) + geom_line(size = 1) +
  labs(title = "Orthostatic Hypotension", x = "Visits", y = "Number of Events") +
  theme_bw() + geom_smooth(method = "lm", color = "pink")


data_LAB$rank <- rank(data_LAB$EVENTDAYS, ties.method = "first")

ggplot(data_LAB, aes(x=EVENTDAYS / 365, y=rank)) + geom_line(size = 1) +
  labs(title = "Orthostatic Hypotension", x = "Years", y = "Number of Events") +
  theme_bw() + geom_smooth(method = "lm", color = "pink")

## Scratch
ggplot(data=data, aes(x=n_OH, y=n_SAE)) + geom_hex() + geom_smooth()