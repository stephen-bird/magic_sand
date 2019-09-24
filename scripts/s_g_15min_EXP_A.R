

# Read Marwan's data, subset and rename 5 columns ------------------------------

# Load Experiment A. Gravel feed begins at 10.83 g/ms and decreases to 0 g/ms,
# while the sand feed starts at 0 g/ms and increases to 10.83 g/ms. The
# experiment was run for 88 hrs with a data point every 15 min.

expA_15min <- read.csv("ExpA-15min.csv", header=TRUE)
expA_15min_sub <- expA_15min[,c(3,8:11)]
colnames(expA_15min_sub) <- c("time","gravel_rate","gravel_feed","sand_rate","sand_feed")
head(expA_15min_sub)
expA_15min_sub <- expA_15min_sub[25:277,]

plot(expA_15min_sub$time,expA_15min_sub$gravel_rate,type="b")
plot(expA_15min_sub$time,expA_15min_sub$sand_rate,type="b")

# Scale the data by the RMS of the distribution:
expA_dat <- expA_15min_sub[,c(2,4)]
norm_expA_dat <- sqrt(colMeans(expA_dat^2,na.rm=T))
expA_dat_n <- sweep(expA_dat, 2, norm_expA_dat, "/")

# Create the timeseries:
expA_ts <- ts(data = expA_dat_n, start = 25, frequency = 4)
head(expA_ts)
tail(expA_ts)
plot(expA_ts)
#expA_ts[1,2] <- NA 
head(expA_ts)

# Set window length (see Golyandina, N., Korobeynikov, A., Shlemov, A., Usevich,
# K., 2013. Multivariate and 2D Extensions of Singular Spectrum Analysis with the
# Rssa Package. arXiv. doi:10.18637/jss.v067.i02):
N_expA <- length(expA_ts[,1])
s_expA <- length(expA_ts[1,])
window_length <- round(s_expA*(N_expA+1)/(s_expA+1))
print(window_length)


# Basic MSSA -------------------------------------------------------------------

s_expA <- ssa(expA_ts, L=window_length, kind="mssa") # Create an SSA object
s_expA # Look inside the SSA object for relevant details
plot(s_expA) # These are the eignevalues
plot(s_expA, type = "vectors", idx = 1:12) # These are the first eigenvectors given by the vector indicies
plot(s_expA, type = "paired", idx = 1:12, plot.contrib = FALSE) # These are the pairs of eigenvectors -- the 'groups' arguement specifies the grouping plot
plot(wcor(s_expA)) # This produces a w-correlation matrix
grouping.auto(s_expA, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the trend:
s_expA_recon <- reconstruct(s_expA, groups = list(1:2)) #list the eigenvectors that represent the trend -- "s_expA_recon" contains the residuals
expA_trend <- s_expA_recon$F1 #feature vector 1 (F1) represents eignevector 1, while F2 represents eigenvector 2... add additional vectors as needed if the trend is represented by multiple eigenvectos (e.g. s_expA_recon$F1 + s_expA_recon$F2...)
plot(expA_trend) #plot the trend
plot(expA_ts)

# Extract the residuals from the trend:
s_expA_res <- residuals(s_expA_recon) #extract the residuals
plot(s_expA_res) #plot the results and look for cyclical behavoir in the data
period_expA <- spec.pgram(s_expA_res, detrend = FALSE, log = "no") #plot a periodogram and look for periods in the data
period_expA # Looks like a period at 25.6 years for gravel and 7.1 years for sand

# Stop and determine the best choice for L:
s2_expA <- ssa(s_expA_res, L=window_length, kind="mssa") # for better seperatebility, set "L" to the maximium value of "N/2" and evenly divisible by the period
plot(s2_expA) #these are the eignevalues
plot(s2_expA, type = "vectors", idx = 1:12)
plot(s2_expA, type = "paired", idx =1:12, plot.contrib = FALSE) #these are the first six eigenvectors given by the vector indicies
plot(wcor(s2_expA, groups = as.list(1:16))) #this produces a w-correlation matrix
grouping.auto(s2_expA, grouping.method = c("wcor"))

# Stop an determine the eigenvectors that represent the seasonal component:
seasons_expA <- reconstruct(s2_expA, groups = list(1:2,5:6,7:8))
parestimate(s2_expA, groups = list(1:2,5:6,7:8), method = "esprit")
plot(seasons_expA)
saveRDS(seasons_expA, "output/seasons_expA.rds")

# For convience, repeat the reconstruction in one step:
recon_expA <- reconstruct(s_expA,
                     groups = list(Trend = 1:2, Seasonality = c(3:4,7:8,9:10)))

recon_expA_final <- plot(recon_expA, plot.method = "xyplot", superpose = FALSE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

print(recon_expA_final)
saveRDS(recon_expA_final, "output/recon_expA_final.rds")


# Examine the noise ------------------------------------------------------------

# Noise Envelope:
res_expA <- residuals(recon_expA) #extract the residuals
env_expA <- ssa(res_expA^2, L=12)
rsd_expA <- sqrt(reconstruct(env_expA, groups=list(1))$F1)

g.noise.plot <- function(res,rsd){
        plot(res[,1], type="l",
             main = "Noise Envelope of Gravel Series",
             xlab = "Time (hrs)",
             ylab = "Scaled Gravel Rate (dimensionless)")
        lines(rsd[,1], type="l",col="red")
        lines(-rsd[,1], type="l",col="red")
}

s.noise.plot <- function(res,rsd){
        plot(res[,2], type="l",
             main = "Noise Envelope of Sand Series",
             xlab = "Time (hrs)",
             ylab = "Scaled Sand Rate (dimensionless)")
        lines(rsd[,2], type="l",col="red")
        lines(-rsd[,2], type="l",col="red")
}


par(mfrow=c(2,1))
g.noise.plot(res_expA,rsd_expA)
s.noise.plot(res_expA,rsd_expA)

saveRDS(res_expA, "output/res_expA.rds")
saveRDS(rsd_expA, "output/rsd_expA.rds")

# Test for white noise
parestimate(s_expA, groups = list(Trend = 1:2, Seasonality = c(3:4,7:8,9:10)),method = "esprit")

# https://robjhyndman.com/hyndsight/ljung-box-test/ For seasonal time series, use h = min(2m,T/5) where T = length of record; m = period of seasonality; h = no. lags to test

Box.test(res_expA[,1],type="Ljung",lag=51)
# Reject the null hypothesis...The gravel data are not independently distributed; they exhibit serial correlation

Box.test(res_expA[,2],type="Ljung",lag=51)
# Reject the null hypothesis...The sand data are not independently distributed; they exhibit serial correlation


boot.mean <- function(x,i){boot.mean <- mean(x[i])}

gres_boot_expA <- boot(abs(res_expA[,1]), boot.mean, R = 10000)
print(gres_boot_expA)
sres_boot_expA <- boot(abs(res_expA[,2]), boot.mean, R = 10000)
print(sres_boot_expA)

# Precent of sigmal explained
1-mean(abs(res_expA[,1]/expA_dat_n[,1]))
1-mean(abs(res_expA[,2]/expA_dat_n[,2]))



# Transform data to orginal units ----------------------------------------------

gravel_signal_n <- data.frame(gravel_signal = expA_trend[,1] + seasons_expA$F1[,1] + seasons_expA$F2[,1] + seasons_expA$F3[,1])
gravel_signal <- gravel_signal_n * norm_expA_dat[1]

sand_signal_n <- data.frame(sand_signal = expA_trend[,2] + seasons_expA$F1[,2] + seasons_expA$F2[,2] + seasons_expA$F3[,2])
sand_signal <- sand_signal_n * norm_expA_dat[2]

gravel_noise_n <- data.frame(gravel_noise = rsd_expA[,1])
gravel_noise <- sweep(gravel_noise_n, 2, norm_expA_dat[1], "*")

sand_noise_n <- data.frame(sand_noise = rsd_expA[,2])
sand_noise <- sweep(sand_noise_n, 2, norm_expA_dat[2], "*")

gravel_trend_n <- data.frame(gravel_trend = expA_trend[,1])
gravel_trend <- sweep(gravel_trend_n, 2, norm_expA_dat[1], "*")

sand_trend_n <- data.frame(sand_trend = expA_trend[,2])
sand_trend <- sweep(sand_trend_n, 2, norm_expA_dat[2], "*")

expA_15min_output <- merge.zoo(gravel_signal,sand_signal,gravel_noise,sand_noise,expA_15min_sub$gravel_rate,expA_15min_sub$sand_rate,expA_15min_sub$gravel_feed,expA_15min_sub$sand_feed)
names(expA_15min_output) <- c("gravel_signal","sand_signal","gravel_noise","sand_noise","gravel_rate","sand_rate","gravel_feed","sand_feed")
expA_15min_output$lower_gravel <- expA_15min_output[,1] - expA_15min_output[,3]
expA_15min_output$upper_gravel <- expA_15min_output[,1] + expA_15min_output[,3]
expA_15min_output$lower_sand <- expA_15min_output[,2] - expA_15min_output[,4]
expA_15min_output$upper_sand <- expA_15min_output[,2] + expA_15min_output[,4]

plot_lims <- expA_15min_output
plot(expA_15min_output$gravel_signal, ylim = c(min(plot_lims),max(plot_lims)), ylab = "Rate (g/ms)", xlab = "Time (hrs)",type = "l",lwd=2,main="Exp A")
lines(expA_15min_output$sand_signal,lwd=2,col="red")
lines(expA_15min_output$upper_gravel,lty=2)
lines(expA_15min_output$lower_gravel,lty=2)
points(expA_15min_output$gravel_rate,pch=1,cex=.5)
lines(expA_15min_output$gravel_feed,type="l",lty=3)
lines(expA_15min_output$upper_sand,lty=2,col="red")
lines(expA_15min_output$lower_sand,lty=2,col="red")
points(expA_15min_output$sand_rate,pch=1,cex=.5,col="red")
lines(expA_15min_output$sand_feed,type="l",lty=3,col="red")

# Output data
write.csv(expA_15min_output,"expA_15min_out.csv")




