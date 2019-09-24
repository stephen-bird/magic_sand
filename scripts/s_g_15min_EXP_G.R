
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment C. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 77 hrs with a data point every 15 min

expG_15min <- read.csv("expG-15min.csv", header=TRUE)
head(expG_15min)
expG_15min_sub <- expG_15min[,c(2,8,11,10,12)]
head(expG_15min_sub)
colnames(expG_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expG_15min_sub)
#good <- complete.cases(expG_15min_sub)
#expG_15min_sub <- expG_15min_sub[good,]
expG_15min_sub <- expG_15min_sub[7:24,]
expG_15min_sub$sand_rate <- ifelse(expG_15min_sub$sand_feed == 0, NA, expG_15min_sub$sand_rate)
head(expG_15min_sub)

plot(expG_15min_sub$time,expG_15min_sub$gravel_rate,type="b")
plot(expG_15min_sub$time,expG_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expG_dat <- expG_15min_sub[,c(2,4)]
norm_expG_dat <- sqrt(colMeans(expG_dat^2,na.rm=T))
expG_dat_n <- sweep(expG_dat, 2, norm_expG_dat, "/")

# Create the timeseries:
expG_ts <- ts(data = expG_dat_n, start = 7, frequency = 4)
head(expG_ts)
tail(expG_ts)
#expG_ts[,2] <- ifelse(expG_ts[,2]==0,NA,expG_ts[,2])
plot(expG_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expG <- length(expG_ts[,1])
s_expG <- length(expG_ts[1,])
window_length <- round(s_expG*(N_expG+1)/(s_expG+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expG <- ssa(expG_ts,  L=12, kind="mssa") # Create an SSA object
s_expG # Look inside the SSA object for relevant details
plot(s_expG) # These are the eignevalues
plot(s_expG, type = "vectors", idx = 1:6) # These are the first eigenvectors given by the vector indicies
plot(s_expG, type = "paired", idx = 1:6, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expG)) # This produces a w-correlation matrix
grouping.auto(s_expG, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expG_recon <- reconstruct(s_expG, groups = list(1)) #list the eigenvectors that represent the trend -- "s_expG_recon" contains the residuals
expG_trend <- s_expG_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expG_recon$F1 + s_expG_recon$F2...)
plot(expG_trend) #plot the trend
plot(expG_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expG_res <- residuals(s_expG_recon) # Extract the residuals
plot(s_expG_res) # Plot the results and look for cyclical behavoir in the data
period_expG <- spec.pgram(s_expG_res, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expG # Looks like 6 hour period for gravel and 4.8 hours for sand

# Stop and determine the best choice for L:
s2_expG <- ssa(s_expG_res, L=12, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expG) #these are the eignevalues
plot(s2_expG, type = "vectors", idx = 1:6)
plot(s2_expG, type = "paired", idx =1:6, plot.contrib = FALSE)
plot(wcor(s2_expG, groups = as.list(1:6)))

grouping.auto(s2_expG, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expG <- reconstruct(s2_expG, groups = list(1:2))
parestimate(s2_expG, groups = list(1:2), method = "esprit")
plot(seasons_expG)
saveRDS(seasons_expG, "output/seasons_expG.rds")

# For convience, repeat the reconstruction in one step:
recon_expG <- reconstruct(s_expG,
                     groups = list(Trend = 1, Seasonality = c(2:3)))

recon_expG_final <- plot(recon_expG, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expG_final)
saveRDS(recon_expG_final, "output/recon_expG_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expG <- residuals(recon_expG) #extract the residuals
env_expG <- ssa(res_expG^2, L=12)
rsd_expG <- sqrt(reconstruct(env_expG, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expG,rsd_expG)
s.noise.plot(res_expG,rsd_expG)

saveRDS(res_expG, "output/res_expG.rds")
saveRDS(rsd_expG, "output/rsd_expG.rds")


# Test for white noise
parestimate(s_expG, groups = list(Trend = 1, Seasonality = c(2:3)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(rsd_expG[,1],type="Ljung",lag=3)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expG[,2],type="Ljung",lag=3)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expG <- boot(abs(rsd_expG[,1]), boot.mean, R = 10000)
print(gres_boot_expG)
sres_boot_expG <- boot(abs(rsd_expG[,2]), boot.mean, R = 10000)
print(sres_boot_expG)

# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expG_trend[,1] + seasons_expG$F1[,1] )
gravel_signal <- gravel_signal_n * norm_expG_dat[1]

sand_signal_n <- data.frame(sand_signal = expG_trend[,2] + seasons_expG$F1[,2]  )
sand_signal <- sand_signal_n * norm_expG_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expG[,1])
gravel_noise <- gravel_noise_n * norm_expG_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expG[,2])
sand_noise <- sand_noise_n * norm_expG_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expG_trend[,1])
gravel_trend <- gravel_trend_n * norm_expG_dat[1]

sand_trend_n <- data.frame(sand_trend = expG_trend[,2])
sand_trend <- sand_trend_n * norm_expG_dat[2]

expG_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expG_15min_sub$gravel_rate,expG_15min_sub$sand_rate,expG_15min_sub$gravel_feed,expG_15min_sub$sand_feed)
names(expG_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expG_15min_output$lower_gravel <- expG_15min_output[,1] - expG_15min_output[,3]
expG_15min_output$upper_gravel <- expG_15min_output[,1] + expG_15min_output[,3]
expG_15min_output$lower_sand <- expG_15min_output[,2] - expG_15min_output[,4]
expG_15min_output$upper_sand <- expG_15min_output[,2] + expG_15min_output[,4]

plot_lims <- ifelse(is.na(expG_15min_output),0,expG_15min_output)
plot(expG_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp G")
lines(expG_15min_output$sand_signal,lwd=2,col="red")
lines(expG_15min_output$upper_gravel,lty=2)
lines(expG_15min_output$lower_gravel,lty=2)
points(expG_15min_output$gravel_rate,pch=1,cex=.5)
lines(expG_15min_output$gravel_feed,type="l",lty=3)
lines(expG_15min_output$upper_sand,lty=2,col="red")
lines(expG_15min_output$lower_sand,lty=2,col="red")
points(expG_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expG_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expG_15min_output,"expG_15min_out.csv")

mean(expG_15min_output$gravel_rate / expG_15min_output$gravel_noise)


