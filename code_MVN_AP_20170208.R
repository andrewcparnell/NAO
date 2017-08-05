# Code for Guiomar model written 8/2/17
# This model fits a different type of MDP model
# It's a joint polynomial model rather than an independent one
# Let's see if it makes a difference

setwd("~/GDrive/Students/Guiomar/paper_20170211")
#setwd("/Volumes/MacintoshHD2/GDrive/Students/Guiomar/paper_20170211")

# Clear the workspace
rm(list=ls())

# Load in packages
pkgs = c("R2jags", "fitdistrplus", "mvtnorm", "tidyverse")
lapply(pkgs, library, character.only = TRUE)

# Get today's date for saving
#today = format(Sys.time(),"%Y%m%d_%H%M")
today = '20170210'

# Start a plot
#pdf(file = paste0('plots_MVN_polynomial_',today,'.pdf'), width = 8, height = 5)

# Load data ---------------------------------------------------------------

load("DATA_NAO_Model.RData")
# We now have:
# NAO_m - 69 modern NAO values
# time_m - 69 modern time values (1934 to 2012 AD)
# xrf_m - 69 by 13 matrix of modern xrf values - all large integers
# xrf_f - 578 by 13 matrix of fossil xrf values
# time_f - 578 fossil time values 173.5BC to 1921.213 AD

# Set up the modern data for use in JAGS
xrf_m_resc = scale(xrf_m)
xrf_m_means = attributes(xrf_m_resc)$`scaled:center`
xrf_m_sds = attributes(xrf_m_resc)$`scaled:scale`
xrf_f_resc = scale(xrf_f, center = xrf_m_means, scale = xrf_m_sds)

# Extra priors for Wishart prior
k = 13
R = diag(k)

data = list(N.rows_m = nrow(xrf_m),
            N.cols_m = ncol(xrf_m),
            xrf_m = xrf_m_resc,
            NAO_m = NAO_m,
            time_m = time_m,
            k = k,
            R = R)

# Fit a polynomial JAGS model to the modern data --------------------------

# Which parameters are we interested in
pars_to_save = c("beta0",
                 "beta1",
                 "beta2",
                 "Sigma_inv",
                 "xrf_m_pred",
                 "sd_rw")

## RUN THE MODEL
if(!file.exists(paste0('run_model1_MVN_polynomial_',today,'.RData'))) {
  run1 = jags(data = data,
              parameters.to.save = pars_to_save,
              model.file = 'modern_MVN_polynomial_model.jags',
              n.chains = 4,
              n.iter = 100000,
              n.burnin = 40000,
              n.thin = 24)
  save(run1,
       file = paste0('run_model1_MVN_polynomial_',today,'.RData'))
} else {
  load(file = paste0('run_model1_MVN_polynomial_',today,'.RData'))
}

plot(run1) # Convergence very good
print(run1)

# Get the predicted mean values
xrf_m_pred = run1$BUGSoutput$mean$xrf_m_pred

# Now for each column of xrf_m plot them against the true values
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
par(mfrow = c(7,2))
for(i in 1:13) {
  plot(data$xrf_m[,i], xrf_m_pred[,i],
       pch = 19,
       xlab = 'True',
       ylab = 'Predicted',
       ylim = range(data$xrf_m[,i]))
  abline(a = 0, b = 1, col = 'red', lty = 'dotted')
}
par(mfrow = c(1,1))
# All of these are useless - no predictive power

# Plot the relationships with uncertainties
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01,las=1)
par(mfrow = c(7,2))
par_means = run1$BUGSoutput$mean
NAO_grid = seq(-3,3,length=50)
for(i in 1:13) {
  plot(data$NAO_m, data$xrf_m[,i],
       pch = 19,
       xlab = 'NAO',
       ylab = 'XRF',
       ylim = range(c(data$xrf_m[,i], xrf_f_resc[,i])))
  mean_pred = with(par_means, beta0[i] + beta1[i] * NAO_grid + beta2[i] * (NAO_grid^2))
  lines(NAO_grid, mean_pred)
  low_pred = with(par_means, beta0[i] + beta1[i] * NAO_grid + beta2[i] * (NAO_grid^2) - 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
  lines(NAO_grid, low_pred, lty = 2)
  high_pred = with(par_means, beta0[i] + beta1[i] * NAO_grid + beta2[i] * (NAO_grid^2) + 2 * sqrt(solve(par_means$Sigma_inv)[i,i]) )
  lines(NAO_grid, high_pred, lty = 2)
  #abline(h = xrf_f_resc[510,i], col = 'red')
}
par(mfrow = c(1,1))

# Have a look at the median of Sigma
sigma_inv_med = run1$BUGSoutput$median$Sigma_inv
sigma_med = solve(sigma_inv_med)
round(sigma_med, 3) # Just PD
plot(eigen(sigma_med)$values) # Only really 2 effective dimensions

# Finally fit a log-normal distribution to sd_rw for later use

## Fit a log-normal distribution to the 50 random data set
results = run1$BUGSoutput$sims.list
lnorm_pars = fitdist(as.vector(results$sd_rw),'lnorm')

# Create MDPs -------------------------------------------------------------

# Create MDPs of NAO_f from NAO_grid
# Run importance sampling on a grid of NAO_f values to get MDPs of NAO_f, one for each layer
# Store these as just the mean and prec for use in the next jags model
pow = function(x, p) return(x^p)
NAO_grid = seq(-3,3,length=50)
# beta0 = apply(results$beta0,2,'median')
# beta1 = apply(results$beta1,2,'median')
# beta2 = apply(results$beta2,2,'median')
# mu0 = apply(results$Sigma_inv,2,'median')
N.rows_f = nrow(xrf_f)

if(!file.exists(paste0('MDP_final_MVN_polynomial_',today,'.Rda'))) {

  MDP_raw = MDP = matrix(0,nrow=N.rows_f,ncol=length(NAO_grid))
  # Create MDPs
  for(i in 1:N.rows_f) {
    print(i)
    for(k in 1:length(NAO_grid)) {
      for(j in 1:nrow(results$beta0)) {
        # Get parameters
        beta0 = results$beta0[j,]
        beta1 = results$beta1[j,]
        beta2 = results$beta2[j,]
        Sigma_inv = results$Sigma_inv[j,,]
        Sigma = solve(Sigma_inv)
        # Create mean
        mean = beta0 + beta1*NAO_grid[k] + beta2*pow(NAO_grid[k],2)
        # Create MDP density (log scale) for current parameters
        MDP_raw[i,k] = MDP_raw[i,k] +
          dmvnorm(xrf_f_resc[i,],
                  mean = mean,
                  sigma = Sigma,
                  log = TRUE)
      }
      # Rescale back into a probabilitiy
      MDP[i,k] = exp(MDP_raw[i,k]/nrow(results$beta0))
    }
  }
  MDP_final = sweep(MDP,1,apply(MDP,1,'sum'),'/')
  save(MDP_final,file=paste0('MDP_final_MVN_polynomial_',today,'.Rda'))
} else {
  load(paste0('MDP_final_MVN_polynomial_',today,'.Rda'))
}

MDP_mean = round(MDP_final,5)%*%NAO_grid # Rounding required otherwise numerical errors
MDP_sd = sqrt(round(MDP_final,5)%*%(NAO_grid^2)-MDP_mean^2)
MDP_prec = 1/(MDP_sd^2)
MDP_high = MDP_mean+2*MDP_sd
MDP_low = MDP_mean-2*MDP_sd

# If any of the MDP_prec values are infinite make them something finite
MDP_prec[is.infinite(MDP_prec)] = 1e10

# Plot MDPs
plot(time_f,MDP_mean,xlim=range(c(time_f,time_m)),ylim=range(c(MDP_low,MDP_high,NAO_m)),col='red',pch=19)
points(time_m,NAO_m,pch=19)
for(i in 1:N.rows_f) {
  lines(c(time_f[i],time_f[i]),c(MDP_low[i],MDP_high[i]),col='red')
}

# Create MDPs on modern data ----------------------------------------------

# Checking process - do the MDPs on the modern data match the true climates
pow = function(x, p) return(x^p)
NAO_grid = seq(-3,3,length=50)
N.rows_m = nrow(xrf_m)

if(!file.exists(paste0('MDP_modern_MVN_',today,'.Rda'))) {

  MDP_modern_raw = MDP_modern = matrix(0,nrow=N.rows_m,ncol=length(NAO_grid))
  # Create MDPs
  for(i in 1:N.rows_m) {
    print(i)
    for(k in 1:length(NAO_grid)) {
      for(j in 1:nrow(results$beta0)) {
        # Get parameters
        beta0 = results$beta0[j,]
        beta1 = results$beta1[j,]
        beta2 = results$beta2[j,]
        Sigma_inv = results$Sigma_inv[j,,]
        Sigma = solve(Sigma_inv)
        # Create mean
        mean = beta0 + beta1*NAO_grid[k] + beta2*pow(NAO_grid[k],2)
        # Create MDP density (log scale) for current parameters
        MDP_modern_raw[i,k] = MDP_modern_raw[i,k] +
          dmvnorm(xrf_m_resc[i,],
                  mean = mean,
                  sigma = Sigma,
                  log = TRUE)
      }
      # Rescale back into a probabilitiy
      MDP_modern[i,k] = exp(MDP_modern_raw[i,k]/nrow(results$beta0))
    }
  }
  MDP_modern_final = sweep(MDP_modern,1,apply(MDP_modern,1,'sum'),'/')
  save(MDP_modern_final,file=paste0('MDP_modern_MVN_',today,'.Rda'))
} else {
  load(paste0('MDP_modern_MVN_',today,'.Rda'))
}

MDP_modern_mean = round(MDP_modern_final,5)%*%NAO_grid # Rounding required otherwise numerical errors
MDP_modern_sd = sqrt(round(MDP_modern_final,5)%*%(NAO_grid^2)-MDP_modern_mean^2)
MDP_modern_prec = 1/(MDP_modern_sd^2)
MDP_modern_high_95 = MDP_modern_mean+qnorm(0.975)*MDP_modern_sd
MDP_modern_low_95 = MDP_modern_mean-qnorm(0.975)*MDP_modern_sd
MDP_modern_high_50 = MDP_modern_mean+qnorm(0.75)*MDP_modern_sd
MDP_modern_low_50 = MDP_modern_mean-qnorm(0.75)*MDP_modern_sd

# If any of the MDP_prec values are infinite make them something finite
MDP_modern_prec[is.infinite(MDP_modern_prec)] = 1e10

# Plot MDPs
#pdf(file = paste0('MDP_modern_plot_',today,'.pdf'), width = 8, height = 5)
cairo_ps(file = paste0('MDP_modern_plot_',today,'.eps'), width = 8, height = 5)
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01, las=1)
plot(time_m,MDP_modern_mean,
     xlim=range(c(time_m,time_m)),
     ylim=range(c(MDP_modern_low_95,MDP_modern_high_95+0.2,NAO_m)),
     xlab = 'Year',
     ylab = 'NAO',
     col='red',
     pch=19)
grid()
points(time_m,NAO_m,pch=19)
for(i in 1:N.rows_m) {
  lines(c(time_m[i],time_m[i]),c(MDP_modern_low_95[i],MDP_modern_high_95[i]),col='red')
  lines(c(time_m[i],time_m[i]),c(MDP_modern_low_50[i],MDP_modern_high_50[i]),col='red', lwd = 3)
}
legend('topleft',legend = c('True NAO', 'Predicted NAO', '95% interval', '50% interval'),
       pch = c(19, 19, -1, -1),
       lty = c(-1, -1, 1, 1),
       lwd = c(-1, -1, 1, 3),
       col = c('black', rep('red', 3)),
       bty = 'n',
       horiz = TRUE)
dev.off()

(sum(NAO_m < MDP_modern_low_95) + sum(NAO_m > MDP_modern_high_95))/ length(NAO_m) # 2.9%
(sum(NAO_m < MDP_modern_low_50) + sum(NAO_m > MDP_modern_high_50))/ length(NAO_m) # 44.9%

# Fit the fossil model ----------------------------------------------------

# Create a time grid
time_grid = seq(-150.51, 2012.51, by = 10) # Don't want any matching times
# Put all the times together
time_fm = sort(c(time_f, time_m))
time_all = sort(c(time_fm,time_grid))
which(diff(time_all)==0) # Should be integer(0)
pick = which(time_all %in% time_fm)
MDP_mean_all = c(MDP_mean[,1], NAO_m)
MDP_prec_all = c(MDP_prec[,1], rep(1e5, length(NAO_m)))
N.rows_fm = length(MDP_mean_all)

data = list(N.rows_fm = N.rows_fm,
            N.rows_all = length(time_all),
            pick = pick,
            MDP_mean = MDP_mean_all,
            MDP_prec = MDP_prec_all,
            time_all = time_all,
            a_rw = lnorm_pars$estimate[1],
            b_rw = 1/pow(lnorm_pars$estimate[2],2))

pars=c("NAO_all","sd_rw")

## RUN THE MODEL
if(!file.exists(paste0('run_model2_MVN_polynomial_',today,'.RData'))) {
  run2 = jags(data=data,
              parameters.to.save=pars,
              model.file='fossil_MDP_model.jags',
              n.chains = 4,
              n.iter = 6000,
              n.burnin = 4000,
              n.thin = 4)
  save(run2,
       file = paste0('run_model2_MVN_polynomial_',today,'.RData'))
} else {
  load(paste0('run_model2_MVN_polynomial_',today,'.RData'))
}

plot(run2)

# Extract out the gridded values
NAO_all = run2$BUGSoutput$sims.list$NAO_all
pick_grid = which(time_all %in% time_grid)
NAO_grid = NAO_all[,pick_grid]

# Create a better plot of interpolated climate
#pdf(file = paste0('NAO_plot_',today,'.pdf'), width = 8, height = 5)
stop()
cairo_ps(file = paste0('NAO_plot_',today,'.eps'), width = 8, height = 5)
par(mfrow = c(1, 1))
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01, las=1)

NAO_all_med = apply(NAO_grid, 2, median)
NAO_all_low_95 = apply(NAO_grid, 2, quantile, 0.025)
NAO_all_high_95 = apply(NAO_grid, 2, quantile, 0.975)
NAO_all_low_50 = apply(NAO_grid, 2, quantile, 0.25)
NAO_all_high_50 = apply(NAO_grid, 2, quantile, 0.75)

plot(time_grid, NAO_all_med,
     type="l",
     xaxt='n',
     ylab = 'Reconstructed NAO',
     xlab = 'Year AD',
     ylim=c(-3,2))

#grid()

axis(1,at=seq(-200,2000,by=200))

polygon(c(time_grid, rev(time_grid)),
        c(NAO_all_low_95, rev(NAO_all_high_95)),
        col = rgb(0.75, 0.75, 0.75, 0.2),
        border = NA)
polygon(c(time_grid, rev(time_grid)),
        c(NAO_all_low_50, rev(NAO_all_high_50)),
        col = rgb(0.75, 0.75, 0.75, 0.5),
        border = NA)

# Add in modern
lines(time_m, NAO_m, col = 'blue')

# Add in horizontal
abline(h=0,col='red')

# Legend
legend('bottomleft', c('Median', '95% interval', '50% interval', 'Modern'),
       lty = c(1, -1, -1, 1),
       col = c(rgb(0, 0, 0),
               rgb(0.75, 0.75, 0.75, 0.2),
               rgb(0.75, 0.75, 0.75, 0.5),
               rgb(0, 0, 1)),
       pch = c(-1, 15, 15, -1),
       bty = 'n')

dev.off()
