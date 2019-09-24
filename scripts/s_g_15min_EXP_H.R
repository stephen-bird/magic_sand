
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment C. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 77 hrs with a data point every 15 min

expH_15min <- read.csv("expH-15min.csv", header=TRUE)
head(expH_15min)
expH_15min_sub <- expH_15min[,c(2,7,10,9,11)]
head(expH_15min)
colnames(expH_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expH_15min_sub)
#good <- complete.cases(expH_15min_sub)
#expH_15min_sub <- expH_15min_sub[good,]
expH_15min_sub <- expH_15min_sub[7:29,]
head(expH_15min_sub)
expH_15min_sub$sand_rate <- ifelse(expH_15min_sub$sand_feed == 0, NA, expH_15min_sub$sand_rate)

plot(expH_15min_sub$time,expH_15min_sub$gravel_rate,type="b")
plot(expH_15min_sub$time,expH_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expH_dat <- expH_15min_sub[,c(2,4)]
norm_expH_dat <- sqrt(colMeans(expH_dat^2,na.rm=T))
expH_dat_n <- sweep(expH_dat, 2, norm_expH_dat, "/")

# Create the timeseries:
expH_ts <- ts(data = expH_dat_n, start = 7, frequency = 4)
head(expH_ts)
tail(expH_ts)
#expH_ts[,2] <- ifelse(expH_ts[,2]==0,NA,expH_ts[,2])
plot(expH_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expH <- length(expH_ts[,1])
s_expH <- length(expH_ts[1,])
window_length <- round(s_expH*(N_expH+1)/(s_expH+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expH <- ssa(expH_ts, L=window_length, kind="mssa") # Create an SSA object
s_expH # Look inside the SSA object for relevant details
plot(s_expH) # These are the eignevalues
plot(s_expH, type = "vectors", idx = 1:8) # These are the first eigenvectors given by the vector indicies
plot(s_expH, type = "paired", idx = 1:8, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expH)) # This produces a w-correlation matrix
grouping.auto(s_expH, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expH_recon <- reconstruct(s_expH, groups = list(1)) #list the eigenvectors that represent the trend -- "s_expH_recon" contains the residuals
expH_trend <- s_expH_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expH_recon$F1 + s_expH_recon$F2...)
plot(expH_trend) #plot the trend
plot(expH_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expH_res <- residuals(s_expH_recon) # Extract the residuals
plot(s_expH_res) # Plot the results and look for cyclical behavoir in the data
period_expH <- spec.pgram(s_expH_res, spans=3, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expH # Looks like 4 hour period for gravel and 3.3 hours for sand

# Stop and determine the best choice for L:
s2_expH <- ssa(s_expH_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expH) #these are the eignevalues
plot(s2_expH, type = "vectors", idx = 1:12)
plot(s2_expH, type = "paired", idx =1:12, plot.contrib = FALSE)
plot(wcor(s2_expH, groups = as.list(1:12)))

grouping.auto(s2_expH, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expH <- reconstruct(s2_expH, groups = list(2:4))
parestimate(s2_expH, groups = list(2:4), method = "esprit")
plot(seasons_expH)
saveRDS(seasons_expH, "output/seasons_expH.rds")

# For convience, repeat the reconstruction in one step:
recon_expH <- reconstruct(s_expH,
                     groups = list(Trend = 1, Seasonality = c(3:4)))

recon_expH_final <- plot(recon_expH, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expH_final)
saveRDS(recon_expH_final, "output/recon_expH_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expH <- residuals(recon_expH) #extract the residuals
env_expH <- ssa(res_expH^2, L=12)
rsd_expH <- sqrt(reconstruct(env_expH, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expH,rsd_expH)
s.noise.plot(res_expH,rsd_expH)

saveRDS(res_expH, "output/res_expH.rds")
saveRDS(rsd_expH, "output/rsd_expH.rds")


# Test for white noise
parestimate(s_expH, groups = list(Trend = 1, Seasonality = c(2:5)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(rsd_expH[,1],type="Ljung",lag=4)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expH[,2],type="Ljung",lag=4)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expH <- boot(abs(rsd_expH[,1]), boot.mean, R = 10000)
print(gres_boot_expH)
sres_boot_expH <- boot(abs(rsd_expH[,2]), boot.mean, R = 10000)
print(sres_boot_expH)

# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expH_trend[,1] + seasons_expH$F1[,1] )
gravel_signal <- gravel_signal_n * norm_expH_dat[1]

sand_signal_n <- data.frame(sand_signal = expH_trend[,2] + seasons_expH$F1[,2]  )
sand_signal <- sand_signal_n * norm_expH_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expH[,1])
gravel_noise <- gravel_noise_n * norm_expH_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expH[,2])
sand_noise <- sand_noise_n * norm_expH_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expH_trend[,1])
gravel_trend <- gravel_trend_n * norm_expH_dat[1]

sand_trend_n <- data.frame(sand_trend = expH_trend[,2])
sand_trend <- sand_trend_n * norm_expH_dat[2]

expH_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expH_15min_sub$gravel_rate,expH_15min_sub$sand_rate,expH_15min_sub$gravel_feed,expH_15min_sub$sand_feed)
names(expH_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expH_15min_output$lower_gravel <- expH_15min_output[,1] - expH_15min_output[,3]
expH_15min_output$upper_gravel <- expH_15min_output[,1] + expH_15min_output[,3]
expH_15min_output$lower_sand <- expH_15min_output[,2] - expH_15min_output[,4]
expH_15min_output$upper_sand <- expH_15min_output[,2] + expH_15min_output[,4]

plot_lims <- ifelse(is.na(expH_15min_output),0,expH_15min_output)
plot(expH_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp H")
lines(expH_15min_output$sand_signal,lwd=2,col="red")
lines(expH_15min_output$upper_gravel,lty=2)
lines(expH_15min_output$lower_gravel,lty=2)
points(expH_15min_output$gravel_rate,pch=1,cex=.5)
lines(expH_15min_output$gravel_feed,type="l",lty=3)
lines(expH_15min_output$upper_sand,lty=2,col="red")
lines(expH_15min_output$lower_sand,lty=2,col="red")
points(expH_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expH_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expH_15min_output,"expH_15min_out.csv")

mean(expH_15min_output$gravel_rate / expH_15min_output$gravel_noise)



