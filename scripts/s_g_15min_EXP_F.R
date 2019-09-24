
# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load expFriment D. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# expFriment was run for 77 hrs with a data point every 15 min

expF_15min <- read.csv("expF-15min.csv", header=TRUE)
head(expF_15min)
expF_15min_sub <- expF_15min[3:20,c(2,7,10,9,10)]
head(expF_15min_sub)
colnames(expF_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expF_15min_sub)
tail(expF_15min_sub)
#good <- complete.cases(expF_15min_sub)
#expF_15min_sub <- expF_15min_sub[good,]
#expF_15min_sub <- expF_15min_sub[16:263,]
#head(expF_15min_sub)

plot(expF_15min_sub$time,expF_15min_sub$gravel_rate,type="b")
plot(expF_15min_sub$time,expF_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expF_dat <- expF_15min_sub[,c(2,4)]
norm_expF_dat <- sqrt(colMeans(expF_dat^2,na.rm=T))
expF_dat_n <- sweep(expF_dat, 2, norm_expF_dat, "/")

# Create the timeseries:
expF_ts <- ts(data = expF_dat_n, start = 1, frequency = 4)
head(expF_ts)
tail(expF_ts)
#expF_ts[,2] <- ifelse(expF_ts[,2]==0,NA,expF_ts[,2])
plot(expF_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expF <- length(expF_ts[,1])
s_expF <- length(expF_ts[1,])
window_length <- round(s_expF*(N_expF+1)/(s_expF+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expF <- ssa(expF_ts, L=window_length, kind="mssa") # Create an SSA object
s_expF # Look inside the SSA object for relevant details
plot(s_expF) # These are the eignevalues
plot(s_expF, type = "vectors", idx = 1:8) # These are the first eigenvectors given by the vector indicies
plot(s_expF, type = "paired", idx = 1:8, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expF)) # This produces a w-correlation matrix
grouping.auto(s_expF, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expF_recon <- reconstruct(s_expF, groups = list(1)) #list the eigenvectors that represent the trend -- "s_expF_recon" contains the residuals
expF_trend <- s_expF_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expF_recon$F1 + s_expF_recon$F2...)
plot(expF_trend) #plot the trend
plot(expF_ts)

# Stop an determine the eigenvectors that represent the trend:
s_expF_res <- residuals(s_expF_recon) # Extract the residuals
plot(s_expF_res) # Plot the results and look for cyclical behavoir in the data
period_expF <- spec.pgram(s_expF_res, detrend = FALSE, log = "no") # Plot a periodogram and look for periods in the data
period_expF # Looks like 4.5 hour period for gravel and NA hours for sand

# Stop and determine the best choice for L:
s2_expF <- ssa(s_expF_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expF) #these are the eignevalues
plot(s2_expF, type = "vectors", idx = 1:8)
plot(s2_expF, type = "paired", idx =1:8, plot.contrib = FALSE)
plot(wcor(s2_expF, groups = as.list(1:8)))

grouping.auto(s2_expF, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expF <- reconstruct(s2_expF, groups = list(2:3))
parestimate(s2_expF, groups = list(2:3), method = "esprit")
plot(seasons_expF)
saveRDS(seasons_expF, "output/seasons_expF.rds")

# For convience, repeat the reconstruction in one step:
recon_expF <- reconstruct(s_expF,
                     groups = list(Trend = 1, Seasonality = c(4:5)))

recon_expF_final <- plot(recon_expF, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expF_final)
saveRDS(recon_expF_final, "output/recon_expF_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expF <- residuals(recon_expF) #extract the residuals
env_expF <- ssa(res_expF^2, L=12)
rsd_expF <- sqrt(reconstruct(env_expF, groups=list(1))$F1)

par(mfrow=c(2,1))
g.noise.plot(res_expF,rsd_expF)
s.noise.plot(res_expF,rsd_expF)

saveRDS(res_expF, "output/res_expF.rds")
saveRDS(rsd_expF, "output/rsd_expF.rds")


# Test for white noise
parestimate(s_expF, groups = list(Trend = 1, Seasonality = c(2:5)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series,
# use h = min(2m,T/5) where T = length of record; m = period of seasonality; h =
# no. lags to test

Box.test(rsd_expF[,1],type="Ljung",lag=4)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(rsd_expF[,2],type="Ljung",lag=4)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expF <- boot(abs(rsd_expF[,1]), boot.mean, R = 10000)
print(gres_boot_expF)
sres_boot_expF <- boot(abs(rsd_expF[,2]), boot.mean, R = 10000)
print(sres_boot_expF)

# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expF_trend[,1] + seasons_expF$F1[,1] )
gravel_signal <- gravel_signal_n * norm_expF_dat[1]

sand_signal_n <- data.frame(sand_signal = expF_trend[,2] + seasons_expF$F1[,2]  )
sand_signal <- sand_signal_n * norm_expF_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expF[,1])
gravel_noise <- gravel_noise_n * norm_expF_dat[1]

sand_noise_n <- data.frame(sand_noise = rsd_expF[,2])
sand_noise <- sand_noise_n * norm_expF_dat[2]

gravel_trend_n <- data.frame(gravel_trend = expF_trend[,1])
gravel_trend <- gravel_trend_n * norm_expF_dat[1]

sand_trend_n <- data.frame(sand_trend = expF_trend[,2])
sand_trend <- sand_trend_n * norm_expF_dat[2]

expF_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expF_15min_sub$gravel_rate,expF_15min_sub$sand_rate,expF_15min_sub$gravel_feed,expF_15min_sub$sand_feed)
names(expF_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expF_15min_output$lower_gravel <- expF_15min_output[,1] - expF_15min_output[,3]
expF_15min_output$upper_gravel <- expF_15min_output[,1] + expF_15min_output[,3]
expF_15min_output$lower_sand <- expF_15min_output[,2] - expF_15min_output[,4]
expF_15min_output$upper_sand <- expF_15min_output[,2] + expF_15min_output[,4]

plot_lims <- ifelse(is.na(expF_15min_output),0,expF_15min_output)
plot(expF_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp F")
lines(expF_15min_output$sand_signal,lwd=2,col="red")
lines(expF_15min_output$upper_gravel,lty=2)
lines(expF_15min_output$lower_gravel,lty=2)
points(expF_15min_output$gravel_rate,pch=1,cex=.5)
lines(expF_15min_output$gravel_feed,type="l",lty=3)
lines(expF_15min_output$upper_sand,lty=2,col="red")
lines(expF_15min_output$lower_sand,lty=2,col="red")
points(expF_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expF_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expF_15min_output,"expF_15min_out.csv")

mean(expF_15min_output$gravel_rate / expF_15min_output$gravel_noise)
