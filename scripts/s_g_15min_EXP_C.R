
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment C. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 77 hrs with a data point every 15 min

expC_15min <- read.csv("ExpC-15min.csv", header=TRUE)
head(expC_15min)
expC_15min_sub <- expC_15min[,c(2,5,7,6,8)]
head(expC_15min_sub)
colnames(expC_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expC_15min_sub)
good <- complete.cases(expC_15min_sub)
expC_15min_sub <- expC_15min_sub[good,]
expC_15min_sub <- expC_15min_sub[16:263,]
expC_15min_sub$sand_rate <- ifelse(expC_15min_sub$sand_feed == 0, NA, expC_15min_sub$sand_rate)
head(expC_15min_sub)

plot(expC_15min_sub$time,expC_15min_sub$gravel_rate,type="b")
plot(expC_15min_sub$time,expC_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expC_dat <- expC_15min_sub[,c(2,4)]
norm_expC_dat <- sqrt(colMeans(expC_dat^2,na.rm=T))
expC_dat_n <- sweep(expC_dat, 2, norm_expC_dat, "/")

# Create the timeseries:
expC_ts <- ts(data = expC_dat_n, start = 15.25, end = 77, frequency = 4)
head(expC_ts)
tail(expC_ts)
#expC_ts[,2] <- ifelse(expC_ts[,2]==0,NA,expC_ts[,2])
plot(expC_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expC <- length(expC_ts[,1])
s_expC <- length(expC_ts[1,])
window_length <- round(s_expC*(N_expC+1)/(s_expC+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expC <- ssa(expC_ts, L=window_length, kind="mssa") # Create an SSA object
s_expC # Look inside the SSA object for relevant details
plot(s_expC) # These are the eignevalues
plot(s_expC, type = "vectors", idx = 1:14) # These are the first eigenvectors given by the vector indicies
plot(s_expC, type = "paired", idx = 1:14, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expC)) # This produces a w-correlation matrix
grouping.auto(s_expC, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expC_recon <- reconstruct(s_expC, groups = list(1:2)) #list the eigenvectors that represent the trend -- "s_expC_recon" contains the residuals
expC_trend <- s_expC_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expC_recon$F1 + s_expC_recon$F2...)
plot(expC_trend) #plot the trend
plot(expC_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expC_res <- residuals(s_expC_recon) # Extract the residuals
plot(s_expC_res) # Plot the results and look for cyclical behavoir in the data
period_expC <- spec.pgram(s_expC_res[17:248,], detrend = FALSE, log = "no", spans = 3) # Plot a periodogram and look for periods in the data
period_expC # Looks like 34.3 hour period for gravel and 4.2 hours for sand

# Stop and determine the best choice for L:
s2_expC <- ssa(s_expC_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expC) #these are the eignevalues
plot(s2_expC, type = "vectors", idx = 1:12)
plot(s2_expC, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expC, groups = as.list(1:16)))

grouping.auto(s2_expC, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expC <- reconstruct(s2_expC, groups = list(1:2,5:6,8:9))
parestimate(s2_expC, groups = list(1:2,5:6,8:9), method = "esprit")
plot(seasons_expC)
saveRDS(seasons_expC, "output/seasons_expC.rds")

# For convience, repeat the reconstruction in one step:
recon_expC <- reconstruct(s_expC,
                     groups = list(Trend = 1:2, Seasonality = c(3:4,7:8,10:11)))

recon_expC_final <- plot(recon_expC, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expC_final)
saveRDS(recon_expC_final, "output/recon_expC_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expC <- residuals(recon_expC) #extract the residuals
env_expC <- ssa(res_expC^2, L=12)
rsd_expC <- sqrt(reconstruct(env_expC, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expC,rsd_expC)
s.noise.plot(res_expC,rsd_expC)

saveRDS(res_expC, "output/res_expC.rds")
saveRDS(rsd_expC, "output/rsd_expC.rds")


# Test for white noise
parestimate(s_expC, groups = list(Trend = 1:2, Seasonality = c(3:4,7:8,10:11)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(rsd_expC[,1],type="Ljung",lag=50)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expC[,2],type="Ljung",lag=50)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expC <- boot(abs(rsd_expC[,1]), boot.mean, R = 10000)
print(gres_boot_expC)
sres_boot_expC <- boot(abs(rsd_expC[38:192,2]), boot.mean, R = 10000)
print(sres_boot_expC)


# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expC_trend[,1] + seasons_expC$F1[,1] + seasons_expC$F2[,1] )
gravel_signal <- gravel_signal_n * norm_expC_dat[1]

sand_signal_n <- data.frame(sand_signal = expC_trend[,2] + seasons_expC$F1[,2] + seasons_expC$F2[,2] )
sand_signal <- sand_signal_n * norm_expC_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expC[,1])
gravel_noise <- gravel_noise_n * norm_expC_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expC[,2])
sand_noise <- sand_noise_n * norm_expC_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expC_trend[,1])
gravel_trend <- gravel_trend_n * norm_expC_dat[1]

sand_trend_n <- data.frame(sand_trend = expC_trend[,2])
sand_trend <- sand_trend_n * norm_expC_dat[2]

expC_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expC_15min_sub$gravel_rate,expC_15min_sub$sand_rate,expC_15min_sub$gravel_feed,expC_15min_sub$sand_feed)
names(expC_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expC_15min_output$lower_gravel <- expC_15min_output[,1] - expC_15min_output[,3]
expC_15min_output$upper_gravel <- expC_15min_output[,1] + expC_15min_output[,3]
expC_15min_output$lower_sand <- expC_15min_output[,2] - expC_15min_output[,4]
expC_15min_output$upper_sand <- expC_15min_output[,2] + expC_15min_output[,4]

plot_lims <- ifelse(is.na(expC_15min_output),0,expC_15min_output)
plot(expC_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp C")
lines(expC_15min_output$sand_signal,lwd=2,col="red")
lines(expC_15min_output$upper_gravel,lty=2)
lines(expC_15min_output$lower_gravel,lty=2)
points(expC_15min_output$gravel_rate,pch=1,cex=.5)
lines(expC_15min_output$gravel_feed,type="l",lty=3)
lines(expC_15min_output$upper_sand,lty=2,col="red")
lines(expC_15min_output$lower_sand,lty=2,col="red")
points(expC_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expC_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expC_15min_output,"expC_15min_out.csv")

mean(expC_15min_output$gravel_rate / expC_15min_output$gravel_noise)

