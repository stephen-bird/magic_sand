

# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment I. Gravel feed is set at 4.33 g/ms,
# while the sand feed starts at 0 g/ms and increases to 13.0 g/ms. The
# experiment was run for 95 hrs with a data point every 15 min.

expI_15min <- read.csv("ExpI-SievedData.csv", header=TRUE)
expI_15min_sub <- expI_15min[,c(3,10,13,12,14)]
colnames(expI_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expI_15min_sub)
good <- complete.cases(expI_15min_sub)
expI_15min_sub <- expI_15min_sub[good,]
expI_15min_sub$sand_rate <- ifelse(expI_15min_sub$sand_feed==0,NA,expI_15min_sub$sand_rate)
expI_15min_sub$sand_rate <- ifelse(expI_15min_sub$sand_feed == 0, NA, expI_15min_sub$sand_rate)
head(expI_15min_sub)

plot(expI_15min_sub$time,expI_15min_sub$gravel_rate,type="b")
plot(expI_15min_sub$time,expI_15min_sub$sand_rate,type="b")

# Standardize data by feed rate:
gravel_std <- expI_15min_sub$gravel_rate / expI_15min_sub$gravel_feed
gravel_std <- ifelse(is.finite(gravel_std)==TRUE,gravel_std,NA)
expI_15min_sub$gravel_std <- gravel_std

sand_std <- expI_15min_sub$sand_rate / expI_15min_sub$sand_feed
sand_std <- ifelse(is.nan(sand_std)==FALSE,sand_std,NA)
expI_15min_sub$sand_std <- sand_std

good <- complete.cases(expI_15min_sub)
expI_15min_sub <- expI_15min_sub[good,]

plot(expI_15min_sub$time,expI_15min_sub$gravel_std,type="b")
plot(expI_15min_sub$time,expI_15min_sub$sand_std,type="b")

expI_dat_n <- select(expI_15min_sub,gravel_std,sand_std)

# Create the timeseries:
expI_ts <- ts(data = expI_dat_n, start = 42.25, end = 95, frequency = 4)
#expI_ts <- ts(data = expI_dat_n, start = 0, end = 52.75, frequency = 4)
head(expI_ts)
tail(expI_ts)
plot(expI_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expI <- length(expI_ts[,1])
s_expI <- length(expI_ts[1,])
window_length <- round(s_expI*(N_expI+1)/(s_expI+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expI <- ssa(expI_ts, L=window_length, kind="mssa") # Create an SSA object
s_expI # Look inside the SSA object for relevant details
plot(s_expI) # These are the eignevalues
plot(s_expI, type = "vectors", idx = 1:14) # These are the first eigenvectors given by the vector indicies
plot(s_expI, type = "paired", idx = 1:14, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expI)) # This produces a w-correlation matrix
grouping.auto(s_expI, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expI_recon <- reconstruct(s_expI, groups = list(1)) #list the eigenvectors that represent the trend -- "s_expI_recon" contains the residuals
expI_trend <- s_expI_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expI_recon$F1 + s_expI_recon$F2...)
plot(expI_trend) #plot the trend
plot(expI_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expI_res <- residuals(s_expI_recon) # Extract the residuals
plot(s_expI_res) # Plot the results and look for cyclical behavoir in the data
period_expI <- spec.pgram(s_expI_res[21:212,], spans=5, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expI # Looks like 38.4 hour period for gravel and 7.4 hours for sand

# Stop and determine the best choice for L:
s2_expI <- ssa(s_expI_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expI) #these are the eignevalues
plot(s2_expI, type = "vectors", idx = 1:12)
plot(s2_expI, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expI, groups = as.list(1:16)))

grouping.auto(s2_expI, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expI <- reconstruct(s2_expI, groups = list(9:10,12:13))
parestimate(s2_expI, groups = list(9:10,12:13),method = "esprit")
plot(seasons_expI)
saveRDS(seasons_expI, "output/seasons_expI.rds")

# For convience, repeat the reconstruction in one step:
recon_expI <- reconstruct(s_expI,
                          groups = list(Trend = 1, Seasonality = c(10:11,13:14)))

recon_expI_final <- plot(recon_expI, plot.method = "xyplot", superpose = FALSE,
                         auto.key = list(columns = 3),
                         col = c("blue", "green", "red", "violet"),
                         lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expI_final)
saveRDS(recon_expI_final, "output/recon_expI_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expI <- residuals(recon_expI) #extract the residuals
env_expI <- ssa(res_expI^2, L=12)
rsd_expI <- sqrt(reconstruct(env_expI, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expI,rsd_expI)
s.noise.plot(res_expI,rsd_expI)

saveRDS(res_expI, "output/res_expI.rds")
saveRDS(rsd_expI, "output/rsd_expI.rds")


# Test for white noise
parestimate(s_expI, groups = list(Trend = 1:2, Seasonality = c(3:4,5:6,7:8)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(res_expI[,1],type="Ljung",lag=42)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(res_expI[,2],type="Ljung",lag=42)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i],na.rm=T)}

gres_boot_expI <- boot(abs(res_expI[,1]), boot.mean, R = 10000)
print(gres_boot_expI)
sres_boot_expI <- boot(abs(res_expI[,2]), boot.mean, R = 10000)
print(sres_boot_expI)


# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expI_trend[,1] + seasons_expI$F1[,1] + seasons_expI$F2[,1] + seasons_expI$F3[,1])
gravel_signal <- gravel_signal_n * norm_expI_dat[1]

sand_signal_n <- data.frame(sand_signal = expI_trend[,2] + seasons_expI$F1[,2] + seasons_expI$F2[,2] + seasons_expI$F3[,2])
sand_signal <- sand_signal_n * norm_expI_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expI[,1])
gravel_noise <- gravel_noise_n * norm_expI_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expI[,2])
sand_noise <- sand_noise_n * norm_expI_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expI_trend[,1])
gravel_trend <- gravel_trend_n * norm_expI_dat[1]

sand_trend_n <- data.frame(sand_trend = expI_trend[,2])
sand_trend <- sand_trend_n * norm_expI_dat[2]

expI_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expI_15min_sub$gravel_rate,expI_15min_sub$sand_rate,expI_15min_sub$gravel_feed,expI_15min_sub$sand_feed)
names(expI_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expI_15min_output$lower_gravel <- expI_15min_output[,1] - expI_15min_output[,3]
expI_15min_output$upper_gravel <- expI_15min_output[,1] + expI_15min_output[,3]
expI_15min_output$lower_sand <- expI_15min_output[,2] - expI_15min_output[,4]
expI_15min_output$upper_sand <- expI_15min_output[,2] + expI_15min_output[,4]

plot_lims <- ifelse(is.na(expI_15min_output),0,expI_15min_output)
plot(expI_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp I")
lines(expI_15min_output$sand_signal,lwd=2,col="red")
lines(expI_15min_output$upper_gravel,lty=2)
lines(expI_15min_output$lower_gravel,lty=2)
points(expI_15min_output$gravel_rate,pch=1,cex=.5)
lines(expI_15min_output$gravel_feed,type="l",lty=3)
lines(expI_15min_output$upper_sand,lty=2,col="red")
lines(expI_15min_output$lower_sand,lty=2,col="red")
points(expI_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expI_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expI_15min_output,"expI_15min_out.csv")

mean(expI_15min_output$gravel_rate / expI_15min_output$gravel_noise)



