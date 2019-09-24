# Load the libraries -----------------------------------------------------------

library(here)
library(Rssa)
library(boot)
library(lubridate)
library(dplyr)
library(zoo)


# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment J

expJ_1sec <- read.csv("ExpJ-processed.csv", skip = 6, header=TRUE)
head(expJ_1sec)
expJ_1sec_sub <- expJ_1sec[,c(1,3,7)]
colnames(expJ_1sec_sub) <- c("time_sec","miss_ratio","bedload_rate")
head(expJ_1sec_sub)
expJ_1sec_sub$bedload_rate <- ifelse(is.nan(expJ_1sec_sub$bedload_rate),NA,expJ_1sec_sub$bedload_rate) # Convert NaN (from Excel) to NA
expJ_1sec_sub$miss_ratio <- ifelse(is.nan(expJ_1sec_sub$miss_ratio),NA,expJ_1sec_sub$miss_ratio) # Convert NaN (from Excel) to NA

# Are there any missing values in the time series itself?
max(expJ_1sec_sub$time_sec - lag(expJ_1sec_sub$time_sec),na.rm = TRUE)
min(expJ_1sec_sub$time_sec - lag(expJ_1sec_sub$time_sec),na.rm = TRUE)
# No, the subtracted lags are all equal to 1

# Filter the data:
good <- complete.cases(expJ_1sec_sub)
expJ_1sec_sub <- expJ_1sec_sub[good,]
expJ_1sec_sub <- filter(expJ_1sec_sub, miss_ratio <1)
head(expJ_1sec_sub)
tail(expJ_1sec_sub)

# Create the timeseries:
expJ_1sec_sub <- ts(expJ_1sec_sub$bedload_rate, start = c(1970,01,01))
plot(expJ_1sec_sub)

# Examime the series for cyclical behavoir:
period_expJ <- spec.pgram(expJ_1sec_sub, spans = 5, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expJ
# Period in the data is 28,800 seconds

# SSA --------------------------------------------------------------------------

# Basic SSA
expJ_s <- ssa(expJ_1sec_sub) # Create an SSA object
expJ_s # Look inside the SSA object for relevant details
plot(expJ_s) # These are the eignevalues; the signal is contained in the first 10 eigentriples (at the break-in-slope of the plot)
plot(expJ_s, type = "vectors", idx = 1:10) # These are the first eigenvectors given by the vector indicies
plot(expJ_s, type = "paired", idx = 1:10, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(expJ_s)) # This produces a w-correlation matrix
grouping.auto(expJ_s, grouping.method = c("wcor"))

# Apply IOSSA
expJ_ios <- iossa(expJ_s, nested.groups = list(1, 2:11), kappa = 2, maxiter = 100, tol = 1e-5)
plot(expJ_ios, type = "vectors", idx = 1:10)
plot(expJ_ios, type = "paired", idx = 1:10, plot.contrib = FALSE)

# Stop an determine the eigenvectors that represent the trend:
expJ_recon <- reconstruct(expJ_ios, groups = list(1)) #list the eigenvectors that represent the trend -- "expJ_recon" contains the residuals
expJ_trend <- expJ_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. expJ_recon$F1 + expJ_recon$F2...)
plot(expJ_trend) #plot the trend
expJ_res <- residuals(expJ_recon) # Extract the residuals
plot(expJ_res) # Plot the results and look for cyclical behavoir in the data
period_expJ <- spec.pgram(expJ_res, spans = 5, detrend = FALSE, log = "yes") # Plot a periodogram and look for periods in the data
period_expJ # Looks like 43200 second period (almost)

# Stop and determine the cyclical behavior in the data:
expJ_s2 <- ssa(expJ_res,L=43200) # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(expJ_s2) #These are the eignevalues; the signal is in the first 9 eigentriples
plot(expJ_s2, type = "vectors", idx = 1:10)
plot(expJ_s2, type = "paired", idx =1:10, plot.contrib = FALSE)
plot(wcor(expJ_s2, groups = as.list(1:10)))
grouping.auto(expJ_s2, grouping.method = c("wcor"))

# Apply IOSSA
expJ_ios_s2 <- iossa(expJ_s2, nested.groups = list(1:2,3:4,5:6,7:8,9:10), kappa = 2, maxiter = 100, tol = 1e-5)
plot(expJ_ios_s2, type = "vectors", idx = 1:10)
plot(expJ_ios_s2, type = "paired", idx = 1:10, plot.contrib = FALSE)
grouping.auto(expJ_ios_s2, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasonexpJ <- reconstruct(expJ_ios_s2, groups = list(1:2,3:4,5:6,7:8,9:10))
parestimate(expJ_ios_s2, groups = list(1:2,3:4,5:6,7:8,9:10),method = "esprit")
plot(seasonexpJ)
saveRDS(seasonexpJ, "output/seasonexpJ.rds")


plot(expJ_trend + seasonexpJ$F1 + seasonexpJ$F2 + seasonexpJ$F3 + seasonexpJ$F4 + seasonexpJ$F5)
plot(expJ_1sec_sub)

# Examine the noise ------------------------------------------------------------

# Noise Envelope:
reexpJ <- residuals(expJ_ios_s2) #extract the residuals
env_expJ <- ssa(reexpJ^2, L=3600)
rsd_expJ <- sqrt(reconstruct(env_expJ, groups=list(1))$F1)


g.noise.plot <- function(res,rsd){
        plot(res, type="l",
             main = "Noise Envelope of Bedload",
             xlab = "Time (seconds)",
             ylab = "Bedload Rate (g/sec)")
        lines(rsd, type="l",col="red")
        lines(-rsd, type="l",col="red")
}


g.noise.plot(reexpJ,rsd_expJ)



saveRDS(reexpJ, "output/reexpJ.rds")
saveRDS(rsd_expJ, "output/rsd_expJ.rds")


# Test for white noise
parestimate(expJ, groups = list(Trend = 1, Seasonality = c(2:3,4:5,6:7,8:9,10:11)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(reexpJ,type="Ljung",lag=17298)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(reexpJ[,2],type="Ljung",lag=42)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation



46285.636


boot.mean <- function(x,i){boot.mean <- mean(x[i],na.rm=T)}

gres_boot_expJ <- boot(abs(reexpJ), boot.mean, R = 10000)
print(gres_boot_expJ)
sres_boot_expJ <- boot(abs(reexpJ[,2]), boot.mean, R = 10000)
print(sres_boot_expJ)







# For convience, repeat the reconstruction in one step:
recon_expJ <- reconstruct(ios, L = 43200,
                          groups = list(Trend = 1, Seasonality = c(2:3,4:5)))

recon_expJ_final <- plot(recon_expJ, plot.method = "xyplot", superpose = FALSE,
                         auto.key = list(columns = 3),
                         col = c("blue", "green", "red", "violet"),
                         lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expJ_final)
saveRDS(recon_expJ_final, "output/recon_expJ_final.rds")


