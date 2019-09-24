# Singular Spectral Analysis of Marwan's flume data containing mixtures of sand
# and gravel. The signal is decomposed into trend, seasonal, and noise. The
# seasonal data are then examined to identify the phase, amplitude, and period
# of the signal. The magntidue of the noise is also examined. Each experiment is
# held in a seperate script. The results are compiled in a Sweave document.

# Stephen Bird

# 2019-05-07
# 2019-08-06

# Load the libraries -----------------------------------------------------------

library(here)
library(Rssa)
library(boot)
library(zoo)


# Read Marwan's data, subset and rename 3 columns ------------------------------

# Load Experiment J. Gravel feed is 11 g/ms and there is no sand feed. The
# experiment was run for 25 hrs with a data point every 15 min.

expJ_15min <- read.csv("ExpJ-SievedData.csv", header=TRUE)
head(expJ_15min)
expJ_15min_sub <- expJ_15min[,c(1,4,5)]
head(expJ_15min_sub)
colnames(expJ_15min_sub) <- c("time","gravel_feed","gravel_rate")
head(expJ_15min_sub)
good <- complete.cases(expJ_15min_sub)
expJ_15min_sub <- expJ_15min_sub[good,]
head(expJ_15min_sub)

plot(expJ_15min_sub$time,expJ_15min_sub$gravel_rate,type="b")

# Scale the data by the RMS of the distribution:
expJ_dat <- expJ_15min_sub[,c(3)]
norm_expJ_dat <- sqrt(mean(expJ_dat^2))
expJ_dat_n <- expJ_dat/norm_expJ_dat
mean(expJ_dat_n)

# Create the timeseries:
expJ_ts <- ts(data = expJ_dat_n, start = 14,  frequency = 4)
head(expJ_ts)
tail(expJ_ts)
plot(expJ_ts)


# Basic SSA --------------------------------------------------------------------

s_expJ <- ssa(expJ_ts)
s_expJ #look inside the SSA object for relevant details
print(ssa.capabilities(s_expJ))
grouping.auto(s_expJ, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
plot(s_expJ) # These are the eignevalues
plot(s_expJ, type = "vectors", idx = 1:10) # These are the first eigenvectors given by the vector indicies -- the first panel is eigenvector 1, etc.
plot(s_expJ, type = "paired", idx = 1:9, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expJ)) # This produces a w-correlation matrix

# Extract the trend:
s_expJ_recon <- reconstruct(s_expJ, groups = list(1)) # List the eigenvectors that represent the trend -- "res1" contains the residuals
expJ_trend <- s_expJ_recon$F1 # Feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed 
plot(expJ_trend)

# Extract the residuals from the trend:
s_expJ_res <- residuals(s_expJ_recon)
plot(s_expJ_res) # Plot the results and look for cyclical behavoir in the data
period_expJ <- spec.pgram(s_expJ_res, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expJ #looks like a period at 12.5 years

# Stop and determine the best choice for L:
s2_expJ <- ssa(s_expJ_res, L=50) # For better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expJ)
plot(s2_expJ, type = "vectors", idx = 1:8)
plot(s2_expJ, type = "paired", idx =1:8, plot.contrib = FALSE)
plot(wcor(s2_expJ, groups = as.list(1:16)))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expJ <- reconstruct(s2_expJ, groups = list(1:2,3:4,5:6))
parestimate(s2_expJ, groups = list(1:2,3:4,5:6), method = "esprit")
plot(seasons_expJ)
saveRDS(seasons_expJ, "output/seasons_expJ.rds")

# For convience, repeat the reconstruction in one step:
recon_expJ <- reconstruct(s_expJ,
                   groups = list(Trend = 1, Seasonality = c(2:3,4:5,6:7)))

recon_expJ_final <- plot(recon_expJ, plot.method = "xyplot", superpose = FALSE,
                            auto.key = list(columns = 3),
                            col = c("blue", "green", "red", "violet"),
                            lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expJ_final)
saveRDS(recon_expJ_final, "output/recon_expJ_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expJ <- residuals(recon_expJ) #extract the residuals
env_expJ <- ssa(res_expJ^2, L=12)
rsd_expJ <- sqrt(reconstruct(env_expJ, groups=list(1))$F1)


noise.plot <- function(res,rsd) {
        plot(res, type="l",
             main = "Noise Envelope of Gravel-only Series",
             xlab = "Time (hrs)",
             ylab = "Scaled Gravel Rate (dimensionless)")
        lines(rsd, type="l",col="red")
        lines(-rsd, type="l",col="red")
        
}

noise.plot(res_expJ,rsd_expJ)
saveRDS(res_expJ, "output/res_expJ.rds")
saveRDS(rsd_expJ, "output/rsd_expJ.rds")


# Test for white noise
parestimate(s_expJ, groups = list(Trend = c(1),
                                   Seasonality = c(2:3,4:5,6:7,8:9)), method = "esprit")
Box.test(res_expJ,type="Ljung",lag=25)
# Reject the null hypothesis...The data are not independently distributed; they exhibit serial correlation

boot.mean <- function(x,i){boot.mean <- mean(x[i])} # Function to bootstrap the mean
trend_boot_expJ <- boot(expJ_trend, boot.mean, R = 10000)
print(trend_boot_expJ)

res_boot_expJ <- boot(abs(res_expJ), boot.mean, R = 10000)
print(res_boot_expJ)

# Precent of sigmal explained
1-mean(abs(res_expJ/expJ_dat_n))


# Transform data to orginal units ----------------------------------------------

signal_n <- data.frame(gravel_signal = expJ_trend + seasons_expJ$F1 + seasons_expJ$F2 + seasons_expJ$F3)
signal <- sweep(signal_n, 2, norm_expJ_dat, "*")

noise_n <- data.frame(gravel_noise = rsd_expJ)
noise <- sweep(noise_n, 2, norm_expJ_dat, "*")

trend_n <- data.frame(gravel_trend = expJ_trend)
trend <- sweep(trend_n, 2, norm_expJ_dat, "*")

expJ_15min_output <- merge.zoo(signal,noise,expJ_15min_sub$gravel_rate,expJ_15min_sub$gravel_feed)
names(expJ_15min_output) <- c("gravel_signal","gravel_noise","gravel_rate","gravel_feed")
expJ_15min_output$lower_gravel <- expJ_15min_output[,1] - expJ_15min_output[,2]
expJ_15min_output$upper_gravel <- expJ_15min_output[,1] + expJ_15min_output[,2]

plot_lims <- expJ_15min_output
plot(expJ_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp J")
lines(expJ_15min_output$upper_gravel,lty=2)
lines(expJ_15min_output$lower_gravel,lty=2)
points(expJ_15min_output$gravel_rate,pch=1,cex=.5)
lines(expJ_15min_output$gravel_feed,type="l",lty=3)

# Output data
write.csv(expJ_15min_output,"expJ_15min_out.csv")




























# Fit a sinusoid to the seasonal curves ----------------------------------------
plot(seasons_expJ)

sinu.stats <- function(timeseries,sinusoid){
        
        wave2analyze <- as.data.frame(sinusoid) # Isolate the cycle of interest
        time_steps <- as.vector(time(timeseries)) # Recreate the timeseries for convience
        wave2analyze <- cbind(time_steps,wave2analyze) 
        names(wave2analyze) <- c("t","y")
        plot(wave2analyze)
        
        ssp <- spectrum(wave2analyze$y)  
        per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
        reslm <- lm(wave2analyze$y ~ sin(2*pi/per*wave2analyze$t)+cos(2*pi/per*wave2analyze$t))
        summary(reslm)
        
        b0 <- coef(reslm)[1]
        alpha <- coef(reslm)[2]
        beta <- coef(reslm)[3]
        
        r <- sqrt(alpha^2 + beta^2) # amplitude
        phi <- atan2(beta, alpha)
        phase_shift <- 2*pi*phi # phase shift
        
        rg <- diff(range(wave2analyze$y))
        #plot(wave2analyze$y~wave2analyze$t,ylim=c(min(wave2analyze$y)-0.1*rg,max(wave2analyze$y)+0.1*rg))
        #lines(fitted(reslm)~wave2analyze$t,col=4,lty=2)   # dashed blue line is sin fit
        c(r,phase_shift)
        
}

sinu.stats(expJ_ts,seasons_expJ$F1)
sinu.stats(expJ_ts,seasons_expJ$F2)
sinu.stats(expJ_ts,seasons_expJ$F3)


# Evaluate the results ---------------------------------------------------------



cbind(time_steps,F3_expJ)





wave2analyze <- as.data.frame(seasons_expJ$F3) # Isolate the cycle of interest
time_steps <- as.vector(time(expJ_ts)) # Recreate the timeseries for convience
wave2analyze <- cbind(time_steps,wave2analyze) 
names(wave2analyze) <- c("t","y")
plot(wave2analyze)

ssp <- spectrum(wave2analyze$y)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(wave2analyze$y ~ sin(2*pi/per*wave2analyze$t)+cos(2*pi/per*wave2analyze$t))
summary(reslm)

b0 <- coef(reslm)[1]
alpha <- coef(reslm)[2]
beta <- coef(reslm)[3]

r <- sqrt(alpha^2 + beta^2) # amplitude
phi <- atan2(beta, alpha)
phase_shift <- 2*pi*phi # phase shift

rg <- diff(range(wave2analyze$y))
plot(wave2analyze$y~wave2analyze$t,ylim=c(min(wave2analyze$y)-0.1*rg,max(wave2analyze$y)+0.1*rg))
lines(fitted(reslm)~wave2analyze$t,col=4,lty=2)   # dashed blue line is sin fit
c(r,phase_shift)
par(mfrow=c(1,2))









