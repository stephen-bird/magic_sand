
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment B. Gravel feed begins at 86.58 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 86.58 g/ms. The
# experiment was run for 51 hrs with a data point every 15 min.

expB_15min <- read.csv("ExpB-15min.csv", header=TRUE)
head(expB_15min)
expB_15min_sub <- expB_15min[,c(2,9,12,11,13)]
head(expB_15min_sub)
colnames(expB_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expB_15min_sub)
good <- complete.cases(expB_15min_sub)
expB_15min_sub <- expB_15min_sub[good,]
expB_15min_sub <- expB_15min_sub[3:195,]
expB_15min_sub$sand_rate <- ifelse(expB_15min_sub$sand_feed == 0, NA, expB_15min_sub$sand_rate)
head(expB_15min_sub)

plot(expB_15min_sub$time,expB_15min_sub$gravel_rate,type="b")
plot(expB_15min_sub$time,expB_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expB_dat <- expB_15min_sub[,c(2,4)]
norm_expB_dat <- sqrt(colMeans(expB_dat^2,na.rm=T))
expB_dat_n <- sweep(expB_dat, 2, norm_expB_dat, "/")

# Create the timeseries:
expB_ts <- ts(data = expB_dat_n, start = 3, frequency = 4)
head(expB_ts)
tail(expB_ts)
plot(expB_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expB <- length(expB_ts[,1])
s_expB <- length(expB_ts[1,])
window_length <- round(s_expB*(N_expB+1)/(s_expB+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expB <- ssa(expB_ts, L=window_length, kind="mssa") # Create an SSA object
s_expB # Look inside the SSA object for relevant details
plot(s_expB) # These are the eignevalues
plot(s_expB, type = "vectors", idx = 1:12) # These are the first eigenvectors given by the vector indicies
plot(s_expB, type = "paired", idx = 1:12, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expB)) # This produces a w-correlation matrix
grouping.auto(s_expB, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expB_recon <- reconstruct(s_expB, groups = list(1:2)) #list the eigenvectors that represent the trend -- "s_expB_recon" contains the residuals
expB_trend <- s_expB_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expB_recon$F1 + s_expB_recon$F2...)
plot(expB_trend) #plot the trend
plot(expB_ts)

# Extract the residuals from the trend:
s_expB_res <- residuals(s_expB_recon) # Extract the residuals
plot(s_expB_res) # Plot the results and look for cyclical behavoir in the data
period_expB <- spec.pgram(s_expB_res, spans = 3, detrend = FALSE, log = "no", na.action = na.omit) # Plot a periodogram and look for periods in the data
period_expB # Looks like a period at 15 hours

# Stop and determine the best choice for L:
s2_expB <- ssa(s_expB_res, L=window_length, kind="mssa") # For better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expB) #these are the eignevalues
plot(s2_expB, type = "vectors", idx = 1:12)
plot(s2_expB, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expB, groups = as.list(1:16)))
grouping.auto(s2_expB, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expB <- reconstruct(s2_expB, groups = list(2:3,4:5))
parestimate(s2_expB, groups = list(2:3,4:5), method = "esprit")
plot(seasons_expB)
saveRDS(seasons_expB, "output/seasons_expB.rds")

# For convience, repeat the reconstruction in one step:
recon_expB <- reconstruct(s_expB,
                     groups = list(Trend = 1:2, Seasonality = c(4:5,6:7)))

recon_expB_final <- plot(recon_expB, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expB_final)
saveRDS(recon_expB_final, "output/recon_expB_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expB <- residuals(recon_expB) #extract the residuals
env_expB <- ssa(res_expB^2, L=12)
rsd_expB <- sqrt(reconstruct(env_expB, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expB,rsd_expB)
s.noise.plot(res_expB,rsd_expB)

saveRDS(res_expB, "output/res_expB.rds")
saveRDS(rsd_expB, "output/rsd_expB.rds")


# Test for white noise
parestimate(s_expB, groups = list(Trend = 1:2, Seasonality = c(4:5,6:7,8:11)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(res_expB[,1],type="Ljung",lag=12)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(res_expB[,2],type="Ljung",lag=12)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expB <- boot(abs(res_expB[,1]), boot.mean, R = 10000)
print(gres_boot_expB)
sres_boot_expB <- boot(abs(res_expB[38:192,2]), boot.mean, R = 10000)
print(sres_boot_expB)


# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expB_trend[,1] + seasons_expB$F1[,1] + seasons_expB$F2[,1] )
gravel_signal <- gravel_signal_n * norm_expB_dat[1]

sand_signal_n <- data.frame(sand_signal = expB_trend[,2] + seasons_expB$F1[,2] + seasons_expB$F2[,2])
sand_signal <- sand_signal_n * norm_expB_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expB[,1])
gravel_noise <- gravel_noise_n * norm_expB_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expB[,2])
sand_noise <- sand_noise_n * norm_expB_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expB_trend[,1])
gravel_trend <- gravel_trend_n * norm_expB_dat[1]

sand_trend_n <- data.frame(sand_trend = expB_trend[,2])
sand_trend <- sand_trend_n * norm_expB_dat[2]

expB_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expB_15min_sub$gravel_rate,expB_15min_sub$sand_rate,expB_15min_sub$gravel_feed,expB_15min_sub$sand_feed)
names(expB_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expB_15min_output$lower_gravel <- expB_15min_output[,1] - expB_15min_output[,3]
expB_15min_output$upper_gravel <- expB_15min_output[,1] + expB_15min_output[,3]
expB_15min_output$lower_sand <- expB_15min_output[,2] - expB_15min_output[,4]
expB_15min_output$upper_sand <- expB_15min_output[,2] + expB_15min_output[,4]

plot_lims <- ifelse(is.na(expB_15min_output),0,expB_15min_output)
plot(expB_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp B")
lines(expB_15min_output$sand_signal,lwd=2,col="red")
lines(expB_15min_output$upper_gravel,lty=2)
lines(expB_15min_output$lower_gravel,lty=2)
points(expB_15min_output$gravel_rate,pch=1,cex=.5)
lines(expB_15min_output$gravel_feed,type="l",lty=3)
lines(expB_15min_output$upper_sand,lty=2,col="red")
lines(expB_15min_output$lower_sand,lty=2,col="red")
points(expB_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expB_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expB_15min_output,"expB_15min_out.csv")






# SNR
plot(expB_15min_output$gravel_rate / expB_15min_output$gravel_noise, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp B")
plot()
lines(expB_15min_output$sand_rate / noise[,2],col="red")
lines(expB_15min_sub$time,expB_15min_sub$gravel_feed,type="l",lty=3)
lines(expB_15min_sub$time,expB_15min_sub$sand_feed,type="l",lty=3,col="red")

plot(expB_15min_output$gravel_rate, expB_15min_output$gravel_noise)
summary(lm(expB_15min_output$gravel_rate ~ expB_15min_output$gravel_noise))
