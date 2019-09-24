
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment D. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 77 hrs with a data point every 15 min

expD_15min <- read.csv("ExpD-15min.csv", header=TRUE)
head(expD_15min)
expD_15min_sub <- expD_15min[,c(2,6,9,8,10)]
head(expD_15min_sub)
colnames(expD_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expD_15min_sub)
tail(expD_15min_sub)
expD_15min_sub <- expD_15min_sub[1:48,]
expD_15min_sub$sand_rate <- ifelse(is.na(expD_15min_sub$sand_feed) & expD_15min_sub$sand_rate == 0, NA, expD_15min_sub$sand_rate)
#good <- complete.cases(expD_15min_sub)
#expD_15min_sub <- expD_15min_sub[good,]
#expD_15min_sub <- expD_15min_sub[16:263,]
#head(expD_15min_sub)

plot(expD_15min_sub$time,expD_15min_sub$gravel_rate,type="b")
plot(expD_15min_sub$time,expD_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expD_dat <- expD_15min_sub[14:48,c(2,4)]
norm_expD_dat <- sqrt(colMeans(expD_dat^2,na.rm=T))
expD_dat_n <- sweep(expD_dat, 2, norm_expD_dat, "/")

# Create the timeseries:
expD_ts <- ts(data = expD_dat_n, start = 14, frequency = 4)
head(expD_ts)
tail(expD_ts)
#expD_ts[,2] <- ifelse(expD_ts[,2]==0,NA,expD_ts[,2])
plot(expD_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expD <- length(expD_ts[,1])
s_expD <- length(expD_ts[1,])
window_length <- round(s_expD*(N_expD+1)/(s_expD+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expD <- ssa(expD_ts, L=window_length, kind="mssa") # Create an SSA object
s_expD # Look inside the SSA object for relevant details
plot(s_expD) # These are the eignevalues
plot(s_expD, type = "vectors", idx = 1:14) # These are the first eigenvectors given by the vector indicies
plot(s_expD, type = "paired", idx = 1:14, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expD)) # This produces a w-correlation matrix
grouping.auto(s_expD, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expD_recon <- reconstruct(s_expD, groups = list(1)) #list the eigenvectors that represent the trend -- "s_expD_recon" contains the residuals
expD_trend <- s_expD_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expD_recon$F1 + s_expD_recon$F2...)
plot(expD_trend) #plot the trend
plot(expD_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expD_res <- residuals(s_expD_recon) # Extract the residuals
plot(s_expD_res) # Plot the results and look for cyclical behavoir in the data
period_expD <- spec.pgram(s_expD_res, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expD # Looks like 3 hour period for gravel and 3.8 hours for sand

# Stop and determine the best choice for L:
s2_expD <- ssa(s_expD_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expD) #these are the eignevalues
plot(s2_expD, type = "vectors", idx = 1:12)
plot(s2_expD, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expD, groups = as.list(1:16)))

grouping.auto(s2_expD, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expD <- reconstruct(s2_expD, groups = list(c(5:6)))
parestimate(s2_expD, groups = list(c(5:6), method = "esprit"))
plot(seasons_expD)
saveRDS(seasons_expD, "output/seasons_expD.rds")

# For convience, repeat the reconstruction in one step:
recon_expD <- reconstruct(s_expD,
                     groups = list(Trend = 1, Seasonality = c(c(6:7))))

recon_expD_final <- plot(recon_expD, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expD_final)
saveRDS(recon_expD_final, "output/recon_expD_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expD <- residuals(recon_expD) #extract the residuals
env_expD <- ssa(res_expD^2, L=12)
rsd_expD <- sqrt(reconstruct(env_expD, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expD,rsd_expD)
s.noise.plot(res_expD,rsd_expD)

saveRDS(res_expD, "output/res_expD.rds")
saveRDS(rsd_expD, "output/rsd_expD.rds")


# Test for white noise
parestimate(s_expD, groups = list(Trend = 1:2, Seasonality = c(c(3,6),4:5)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series,
# use h = min(2m,T/5) where T = length of record; m = period of seasonality; h =
# no. lags to test

Box.test(rsd_expD[,1],type="Ljung",lag=6)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expD[,2],type="Ljung",lag=6)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expD <- boot(abs(rsd_expD[,1]), boot.mean, R = 10000)
print(gres_boot_expD)
sres_boot_expD <- boot(abs(rsd_expD[,2]), boot.mean, R = 10000)
print(sres_boot_expD)

# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expD_trend[,1] + seasons_expD$F1[,1])
gravel_signal <- gravel_signal_n * norm_expD_dat[1]

sand_signal_n <- data.frame(sand_signal = expD_trend[,2] + seasons_expD$F1[,2] )
sand_signal <- sand_signal_n * norm_expD_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expD[,1])
gravel_noise <- gravel_noise_n * norm_expD_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expD[,2])
sand_noise <- sand_noise_n * norm_expD_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expD_trend[,1])
gravel_trend <- gravel_trend_n * norm_expD_dat[1]

sand_trend_n <- data.frame(sand_trend = expD_trend[,2])
sand_trend <- sand_trend_n * norm_expD_dat[2]

expD_15min_sub <- expD_15min_sub[14:48,]

expD_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expD_15min_sub$gravel_rate,expD_15min_sub$sand_rate,expD_15min_sub$gravel_feed,expD_15min_sub$sand_feed)
names(expD_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expD_15min_output$lower_gravel <- expD_15min_output[,1] - expD_15min_output[,3]
expD_15min_output$upper_gravel <- expD_15min_output[,1] + expD_15min_output[,3]
expD_15min_output$lower_sand <- expD_15min_output[,2] - expD_15min_output[,4]
expD_15min_output$upper_sand <- expD_15min_output[,2] + expD_15min_output[,4]

plot_lims <- ifelse(is.na(expD_15min_output),0,expD_15min_output)
plot(expD_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp D")
lines(expD_15min_output$sand_signal,lwd=2,col="red")
lines(expD_15min_output$upper_gravel,lty=2)
lines(expD_15min_output$lower_gravel,lty=2)
points(expD_15min_output$gravel_rate,pch=1,cex=.5)
lines(expD_15min_output$gravel_feed,type="l",lty=3)
lines(expD_15min_output$upper_sand,lty=2,col="red")
lines(expD_15min_output$lower_sand,lty=2,col="red")
points(expD_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expD_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expD_15min_output,"expD_15min_out.csv")

mean(expD_15min_output$gravel_rate / expD_15min_output$gravel_noise)
