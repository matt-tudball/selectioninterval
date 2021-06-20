rm(list=ls())
library(devtools); library(foreign); library(ggplot2)
setwd("C:/Users/ow18301/OneDrive - University of Bristol/MyFiles-Migrated/Documents/IEASBOS_FILES_")

#p <- c(0.02, 0.2, 2) # 1
#p <- c(0.01, 0.5, 2) # 2
#p <- c(0.02, 0.2, 2.5) # 3
p <- c(0.01, 0.5, 2.5) # 4
#p <- c(0.02, 0.2, 1.75) # 5
#p <- c(0.01, 0.5, 1.75) # 6

# Load function
load_all(path='CODE/selectioninterval')

# Load data
dat <- read.dta("DATA/data_for_example.dta")

# Format data
x <- dat[,colnames(dat)[startsWith(colnames(dat),'x')]]
z <- dat[,c('z_rosla', colnames(x)[colnames(x) != 'x_edu'])]
w <- dat[,colnames(dat)[startsWith(colnames(dat),'w')]]
y <- dat[,'y_inc']

# Interaction term
w$int <- w$w_edu*w$w_inc

# Constraints
mycons1 <- list(
  list('RESP', 0.055)
)

mycons2 <- list(
  list('RESP', 0.055),
  list('COVMEAN', w[,2], 0.495)
)

mycons3 <- list(
  list('RESP', 0.055),
  list('COVMEAN', w[,2], 0.495),
  list('COVMEAN', w[,3], 0.21)
)

mycons4 <- list(
  list('RESP', 0.055),
  list('COVMEAN', w[,2], 0.495),
  list('COVMEAN', w[,3], 0.21),
  list('COVMEAN', w[,4], -0.297)
)

#mycons5 <- list(
#  list('RESP', 0.055),
#  list('COVMEAN', w[,2], 0.495),
#  list('COVMEAN', w[,3], 0.21),
#  list('COVMEAN', w[,4], -0.297),
#  list('DIREC', 1, '+')
#)

# Intervals
out0 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=NULL)

out1 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=mycons1)

out2 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=mycons2)

out3 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=mycons3)

out4 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=mycons4)

#out5 <- selection_bound(y=y, x=x, z=z, w=w, L0l=p[1], L0u=p[2], L1=p[3], cons=mycons5)

# Statistics implied by weights
u <- as.matrix(w); u <- apply(u, 2, function(v) (v - mean(v))/sd(v))
inv_wgt_min <- 1+exp(as.matrix(-cbind(1,u))%*%out4$theta_min)
inv_wgt_max <- 1+exp(as.matrix(-cbind(1,u))%*%out4$theta_max)

# For min
print(paste("Population mean of males:", round(weighted.mean(x=w$w_sex, w=inv_wgt_min),4)))
print(paste("Population mean of income:", round(weighted.mean(x=w$w_inc, w=inv_wgt_min),4)))
print(paste("Population mean of age:", round(weighted.mean(x=w$w_rosla, w=inv_wgt_min),4)))

# For max
print(paste("Population mean of males:", round(weighted.mean(x=w$w_sex, w=inv_wgt_max),4)))
print(paste("Population mean of income:", round(weighted.mean(x=w$w_inc, w=inv_wgt_max),4)))
print(paste("Population mean of age:", round(weighted.mean(x=w$w_rosla, w=inv_wgt_max),4)))

# Plot
main <- data.frame(type=c("0","1","2","3","4"),
                   ci_lower=c(out0$ci[1], out1$ci[1], out2$ci[1], out3$ci[1], out4$ci[1]),
                   ci_upper=c(out0$ci[2], out1$ci[2], out2$ci[2], out3$ci[2], out4$ci[2]),
                   interval_lower=c(out0$interval[1], out1$interval[1], out2$interval[1], out3$interval[1], out4$interval[1]),
                   interval_upper=c(out0$interval[2], out1$interval[2], out2$interval[2], out3$interval[2], out4$interval[2]))

saveRDS(main, file="SIMULATIONS/DATA/applied_example_4.rds")

plot <- ggplot(main, aes(type)) +
  geom_hline(aes(yintercept=0), color = "grey", linetype="dashed") +
  geom_linerange(
    aes(ymin = interval_lower,
        ymax = interval_upper,
        color = as.factor(type)),
    position = position_dodge(0.3), size=3
  ) +
  geom_errorbar(
    aes(ymin = ci_lower,
        ymax = ci_upper,
        color = as.factor(type)),
    position = position_dodge(0.3), width=0.2
  ) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ylab("Intervals for the effect estimate") +
  scale_colour_discrete(name="Constraints", labels=c("None","1","2","3","4")) +
  coord_flip()

plot
ggsave(filename=paste("FIGURES/applied_example_4",type=".pdf",sep=""),plot=plot)

