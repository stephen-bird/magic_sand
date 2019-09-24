
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment E. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 77 hrs with a data point every 15 min

expE_15min <- read.csv("expE-15min.csv", header=TRUE)
head(expE_15min)
expE_15min_sub <- expE_15min[,c(2,7,10,9,10)]
head(expE_15min_sub)
colnames(expE_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expE_15min_sub)
tail(expE_15min_sub)
#good <- complete.cases(expE_15min_sub)
#expE_15min_sub <- expE_15min_sub[good,]
#expE_15min_sub <- expE_15min_sub[16:263,]
expE_15min_sub$sand_rate <- ifelse(expE_15min_sub$sand_feed == 0, NA, expE_15min_sub$sand_rate)
#head(expE_15min_sub)

plot(expE_15min_sub$time,expE_15min_sub$gravel_rate,type="b")
plot(expE_15min_sub$time,expE_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expE_dat <- expE_15min_sub[4:28,c(2,4)]
norm_expE_dat <- sqrt(colMeans(expE_dat^2,na.rm=T))
expE_dat_n <- sweep(expE_dat, 2, norm_expE_dat, "/")

# Create the timeseries:
expE_ts <- ts(data = expE_dat_n, start = 1, end = 7, frequency = 4)
head(expE_ts)
tail(expE_ts)
#expE_ts[,2] <- ifelse(expE_ts[,2]==0,NA,expE_ts[,2])
plot(expE_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expE <- length(expE_ts[,1])
s_expE <- length(expE_ts[1,])
window_length <- round(s_expE*(N_expE+1)/(s_expE+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expE <- ssa(expE_ts, L=window_length, kind="mssa") # Create an SSA object
s_expE # Look inside the SSA object for relevant details
plot(s_expE) # These are the eignevalues
plot(s_expE, type = "vectors", idx = 1:14) # These are the first eigenvectors given by the vector indicies
plot(s_expE, type = "paired", idx = 1:14, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expE)) # This produces a w-correlation matrix
grouping.auto(s_expE, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expE_recon <- reconstruct(s_expE, groups = list(1:2)) #list the eigenvectors that represent the trend -- "s_expE_recon" contains the residuals
expE_trend <- s_expE_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expE_recon$F1 + s_expE_recon$F2...)
plot(expE_trend) #plot the trend
plot(expE_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expE_res <- residuals(s_expE_recon) # Extract the residuals
plot(s_expE_res) # Plot the results and look for cyclical behavoir in the data
period_expE <- spec.pgram(s_expE_res, spans=3,detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expE # Looks like 6.25 hour period for gravel and 4.2 hours for sand

# Stop and determine the best choice for L:
s2_expE <- ssa(s_expE_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expE) #these are the eignevalues
plot(s2_expE, type = "vectors", idx = 1:12)
plot(s2_expE, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expE, groups = as.list(1:16)))

grouping.auto(s2_expE, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expE <- reconstruct(s2_expE, groups = list(1:2,3:4))
parestimate(s2_expE, groups = list(1:2,3:4), method = "esprit")
plot(seasons_expE)
saveRDS(seasons_expE, "output/seasons_expE.rds")

# For convience, repeat the reconstruction in one step:
recon_expE <- reconstruct(s_expE,
                     groups = list(Trend = 1:2, Seasonality = c(c(3:4),5:6)))

recon_expE_final <- plot(recon_expE, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expE_final)
saveRDS(recon_expE_final, "output/recon_expE_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expE <- residuals(recon_expE) #extract the residuals
env_expE <- ssa(res_expE^2, L=12)
rsd_expE <- sqrt(reconstruct(env_expE, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expE,rsd_expE)
s.noise.plot(res_expE,rsd_expE)

saveRDS(res_expE, "output/res_expE.rds")
saveRDS(rsd_expE, "output/rsd_expE.rds")


# Test for white noise
parestimate(s_expE, groups = list(Trend = 1:2, Seasonality = c(3,4,5:6)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series,
# use h = min(2m,T/5) where T = length of record; m = period of seasonality; h =
# no. lags to test

Box.test(rsd_expE[,1],type="Ljung",lag=5)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expE[,2],type="Ljung",lag=5)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expE <- boot(abs(rsd_expE[,1]), boot.mean, R = 10000)
print(gres_boot_expE)
sres_boot_expE <- boot(abs(rsd_expE[,2]), boot.mean, R = 10000)
print(sres_boot_expE)


# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expE_trend[,1] + seasons_expE$F1[,1] + seasons_expE$F2[,1] )
gravel_signal <- gravel_signal_n * norm_expE_dat[1]

sand_signal_n <- data.frame(sand_signal = expE_trend[,2] + seasons_expE$F1[,2] + seasons_expE$F2[,2] )
sand_signal <- sand_signal_n * norm_expE_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expE[,1])
gravel_noise <- gravel_noise_n * norm_expE_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expE[,2])
sand_noise <- sand_noise_n * norm_expE_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expE_trend[,1])
gravel_trend <- gravel_trend_n * norm_expE_dat[1]

sand_trend_n <- data.frame(sand_trend = expE_trend[,2])
sand_trend <- sand_trend_n * norm_expE_dat[2]

expE_15min_sub <- expE_15min_sub[4:28,]

expE_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expE_15min_sub$gravel_rate,expE_15min_sub$sand_rate,expE_15min_sub$gravel_feed,expE_15min_sub$sand_feed)
names(expE_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expE_15min_output$lower_gravel <- expE_15min_output[,1] - expE_15min_output[,3]
expE_15min_output$upper_gravel <- expE_15min_output[,1] + expE_15min_output[,3]
expE_15min_output$lower_sand <- expE_15min_output[,2] - expE_15min_output[,4]
expE_15min_output$upper_sand <- expE_15min_output[,2] + expE_15min_output[,4]

plot_lims <- ifelse(is.na(expE_15min_output),0,expE_15min_output)
plot(expE_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp E")
lines(expE_15min_output$sand_signal,lwd=2,col="red")
lines(expE_15min_output$upper_gravel,lty=2)
lines(expE_15min_output$lower_gravel,lty=2)
points(expE_15min_output$gravel_rate,pch=1,cex=.5)
lines(expE_15min_output$gravel_feed,type="l",lty=3)
lines(expE_15min_output$upper_sand,lty=2,col="red")
lines(expE_15min_output$lower_sand,lty=2,col="red")
points(expE_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expE_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expE_15min_output,"expE_15min_out.csv")

mean(expE_15min_output$gravel_rate / expE_15min_output$gravel_noise)



