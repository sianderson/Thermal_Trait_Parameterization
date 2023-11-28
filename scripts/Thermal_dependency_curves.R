# Code for Anderson et al. 2023 
# "Phytoplankton thermal trait parameterization alters community structure and biogeochemical processes in a modeled ocean"
# Stephanie I. Anderson updated 11/21/23

# Contains: 
# Table 1: Model Parameters
# Table S1: a_PCmax
# Figure 1: Thermal Dependency curves for each model and PFT
# Figure S2 & S3: exponential curve fits for Kremer (S2) and Anderson (S3)
# Figure S4: allometrically scaled maximum growth curves

###########################################################################
# Load packages
library(ggplot2)
library(quantreg)
library(lme4)
library(dplyr)
library(data.table)
library(cowplot)
library(bbmle)

# data from Anderson et al. (2021) + Diazotrophs and Green algae from Kremer et al (2017)
setwd("/Users/Stephanie/Desktop/MIT/Q10_variability/Code & Datasets/")
isolates <- read.csv("data/derived_traits_diazogreen.csv") # information on isolates
rates <- read.csv("data/growth_rates_diazogreen.csv") # growth rates
group = list("diatoms","coccolithophores", "dinoflagellates","cyanobacteria", "greens", "diazotrophs")

###########################################################################
##### Fitting Exponential Thermal Dependencies ####
###########################################################################

##### Kremer Method
# Weighted Q10s by number of measurements within each group 
rate2 <- rates
rate2$n <- ifelse(rates$group == "coccolithophores", 202, 
                  ifelse(rates$group == "cyanobacteria", 502,
                         ifelse(rates$group == "diatoms", 1794,
                                ifelse(rates$group == "diazotrophs", 144,
                                       ifelse(rates$group == "dinoflagellates", 748,175)))))
rate2$weight = 1/rate2$n
biss3 <- rq(ln.r~temperature+group, data=rate2,tau=0.99,ci=T, weights=weight) #weighted by group
cf.bd3 <- coef(biss3)
qr<- boot.rq(rate2$temperature,
             rate2$ln.r,tau=0.99, R=10000, method="mcmb")
ci_k <- t(apply(qr$B, 2, quantile, c(0.025,0.975)))
se = summary(biss3, se = "boot", bsmethod = "mcmb")$coefficients


######## Eppley
# weighting by number of growth rates for each strain in the database so no strain is overrepresented 
rates <-ddply(rates,. (isolate.code), mutate, n=length(isolate.code), weight=1/n)
biss6 <- rq(ln.r~temperature, data=rates,tau=0.99,ci=T, weights=weight) 
cf.bd6 <- coef(biss6)
qr<- boot.rq(rates$temperature,
             rates$ln.r,tau=0.99, R=10000, method="mcmb")
ci_e <- t(apply(qr$B, 2, quantile, c(0.025,0.975)))

###########################################################################
##### Table 1 ####
###########################################################################
# Parameter confidence intervals are calculated here, though they are not included in Table 1

######## Kremer
group <- c("coccolithophores", "cyanobacteria", "diatoms", "diazotrophs","dinoflagellates","greens")
{model <- c(rep("Kremer", 6))
table_sameQ10 = as.data.frame(model)
table_sameQ10$group = group
table_sameQ10$n <- rbind(length(subset(isolates, group == "coccolithophores")$isolate.code), 
                         length(subset(isolates, group == "cyanobacteria")$isolate.code),
                         length(subset(isolates, group == "diatoms")$isolate.code),
                         length(subset(isolates, group == "diazotrophs")$isolate.code),
                         length(subset(isolates, group == "dinoflagellates" )$isolate.code),
                         length(subset(isolates, group == "greens")$isolate.code))
table_sameQ10$N <- rbind(length(subset(rates, group == "coccolithophores")$r), 
                         length(subset(rates, group == "cyanobacteria")$r),
                         length(subset(rates, group == "diatoms")$r),
                         length(subset(rates, group == "diazotrophs")$r),
                         length(subset(rates, group == "dinoflagellates" )$r),
                         length(subset(rates, group == "greens")$r))
table_sameQ10$a <- rbind(round(cf.bd3[[1]],4), round(cf.bd3[[1]]+cf.bd3[[3]],4), round(cf.bd3[[1]]+cf.bd3[[4]],4), 
                         round(cf.bd3[[1]]+cf.bd3[[5]],4), round(cf.bd3[[1]]+cf.bd3[[6]],4),round(cf.bd3[[1]]+cf.bd3[[7]],4))
## Confidence Intervals (a)
# table_sameQ10$a_cil <- rbind(round(cf.bd3[[8]],4), round(cf.bd3[[8]]+CI(cyano, c=0.95)[[3]],4),
#                              round(CI(ints, c=0.95)[[3]]+CI(diatoms, c=0.95)[[3]],4), round(CI(ints, c=0.95)[[3]]+CI(diazos, c=0.95)[[3]],4),
#                              round(CI(ints, c=0.95)[[3]]+CI(dinos, c=0.95)[[3]],4),round(CI(ints, c=0.95)[[3]]+CI(greens, c=0.95)[[3]],4))
# table_sameQ10$a_ciu <- rbind(round(CI(ints, c=0.95)[[1]],4), round(CI(ints, c=0.95)[[1]]+CI(cyano, c=0.95)[[1]],4),
#                              round(CI(ints, c=0.95)[[1]]+CI(diatoms, c=0.95)[[1]],4), round(CI(ints, c=0.95)[[1]]+CI(diazos, c=0.95)[[1]],4),
#                              round(CI(ints, c=0.95)[[1]]+CI(dinos, c=0.95)[[1]],4),round(CI(ints, c=0.95)[[1]]+CI(greens, c=0.95)[[1]],4))
table_sameQ10$b <- round(cf.bd3[[2]], 4)
## Confidence Intervals (b)
#table_sameQ10$b_cil <- round(CI(coeff, ci=0.95)[[3]], 4)
#table_sameQ10$b_ciu <- round(CI(coeff, ci=0.95)[[1]], 4)
table_sameQ10$Q10 <- round(exp(cf.bd3[[2]]*10),4)
table_sameQ10$intercept <- round(exp(table_sameQ10$a),4)

######## Eppley
model <- c(rep("Eppley", 6))
table_sameumax = as.data.frame(model)
table_sameumax$group = group
table_sameumax$n = length(isolates$isolate.code)
table_sameumax$N = length(rates$isolate.code)
table_sameumax$a <- round(cf.bd6[[1]],4)
#table_sameumax$a_ci <- paste0("[",round(cf.bd6[[3]],3),", ",round(cf.bd6[[5]],3),"]") # confidence intervals (a)
table_sameumax$b <- round(cf.bd6[[2]], 4)
#table_sameumax$b_ci <- paste0("[",round(cf.bd6[[4]],3),", ",round(cf.bd6[[6]],3),"]") # confidence intervals (b)
table_sameumax$Q10 <- round(exp(cf.bd6[[2]]*10),4)
table_sameumax$intercept <- round(exp(table_sameumax$a),4)
table <- rbind(table_sameumax, table_sameQ10) # removed table_sameQ10_un
table$umax20 = exp(table$a+table$b*20) }


######## Table from Anderson et al. 2021
{
  table_natcomm <- read.csv("data/Anderson_etal_2021_table1.csv")
  table_natcomm$model = "Anderson"
  table_natcomm$weighted = "no"
  table1<-rbindlist(list(table, table_natcomm), fill=TRUE,use.names=TRUE)
  table1 = table1[c(1:12, 14:17),1:9]
  table1 <- rbind(table1, table1[10,]) # copying the diazotrophs from Kremer
  table1[17, 1] = "Anderson"
}

# Adding green algae for Anderson model
{
  ratesg = rates
  bissd_g<-rq(ln.r~temperature, data=subset(ratesg, group=="greens"),tau=0.99,ci=T)
  cf.bg<-coef(bissd_g)
  QR.g <- boot.rq(cbind(1,subset(ratesg, group=="greens")$temperature),
                (subset(ratesg, group=="greens")$ln.r),tau=0.99, R=10000, method="mcmb")
  ci_g<-t(apply(QR.g$B, 2, quantile, c(0.025,0.975)))
  model <- "Anderson"
  table_green = as.data.frame(model)
  table_green$group = "greens"
  table_green$n = length(subset(isolates, group == "greens")$isolate.code)
  table_green$N = length(subset(rates, group == "greens")$r)
  table_green$a = round(cf.bg[[1]],4)
  table_green$b = round(cf.bg[[2]],4)
  table_green$Q10 = round(exp(cf.bg[[2]]*10),4)
  table_green$intercept <- round(exp(table_green$a),4)
  table_green$umax20 = exp(table_green$a+table_green$b*20)

  table1 = rbind(table1, table_green)
  table1$intercept = round(table1$intercept, 4)
  table1$umax20 = round(table1$umax20, 4)
}

#write.csv(table1, "output/table1.csv")
table1 <- read.csv("output/table1.csv")


################################################################################
#### a_PCmax as in Table S1
################################################################################

# biovolumes of all tracers (i.e. plankton) in the MIT Darwin model
biovol= c(0.125892541179417     ,  0.410204102986607     ,   1.33659551654644     ,   4.35511873685569     ,
          14.1905752168909     ,   46.2381021399260     ,   150.660706618674     ,   490.907876152603     ,
          1599.55802861467     , 14.1905752168909     ,   46.2381021399260     ,   150.660706618674     ,
          490.907876152603     ,   1599.55802861467     , 14.1905752168909     ,   46.2381021399260     ,
          150.660706618674     ,   490.907876152603     ,   1599.55802861467     , 5211.94711105080     ,
          16982.4365246174     ,   55335.0109215737     ,   180301.774085957     ,   150.660706618674     ,
          490.907876152603     ,   1599.55802861467     ,   5211.94711105080     ,   16982.4365246174     ,
          55335.0109215737     , 180301.774085957     ,   587489.352529777     ,   46.2381021399261     ,
          150.660706618674     ,   490.907876152603     , 1599.55802861467     ,   5211.94711105081     ,
          16982.4365246175     ,   55335.0109215737     ,   180301.774085957     , 587489.352529777     ,
          1914255.92502109     ,   6237348.35482419     ,   20323570.1093622     ,   66221650.3701762     ,
          215774440.915267     ,   703072319.883834     ,   2290867652.76777     ,  3.863669770540691E-002,
          0.125892541179417     , 0.410204102986607)

# Numbers corresponding to the biovolume of each PFT phenotype
# 1-2: cyano
# 3-4: greens
# 5-9: coccolithophores
# 10-14: diazotroph
# 15-23: diatoms
# 24-31: dinos
table1$b_coeff =  ifelse(table1$group %in% c("cyanobacteria", "greens"), 0.28, -0.1) # if 
table1$mumax_biovolnum = ifelse(table1$group == "cyanobacteria", 2, 
                                ifelse(table1$group == "greens", 4,
                                       ifelse(table1$group == "coccolithophores", 5,
                                              ifelse(table1$group == "diazotrophs", 10,                                                     ifelse(table1$group == "diatoms", 15, 24)))))

table1$a_PCmax = 0
for(i in 1:length(table1$a_PCmax)){
  table1$a_PCmax[i] = round(table1$umax20[i]/(biovol[table1$mumax_biovolnum[i]]^table1$b_coeff[i]),4)
}
table1 <- table1[ , -which(names(table1) %in% c("mumax_biovolnum","b_coeff", "notes"))]

#write.csv(table1, "output/tableS1.csv")


################################################################################
#### Figure 1
################################################################################

x<-seq(-2, 35, by=0.1)
eq = function(a,b,x){exp(a+b*x)}

{pdf("figures/Figure1.pdf", width = 7.2, height = 3)
par(mfrow = c(1, 3),mar=c(1,1, 0.5, 0.5), oma=c(3,3,1,1), bg = "transparent")

# Eppley
plot(eq(cf.bd6[[1]],	cf.bd6[[2]], x), type='l', x=x,lwd=2.5, ylim=c(0,4),xlab='',ylab = "")
points(rates$temperature, rates$r, col = alpha("black", 0.1), pch=16)
text(0,3.75,quote(bold("(a)") * " Eppley"), line=1, adj=0.05, cex=1.2)

# Kremer
par(mar=c(1,1,0.5,0.5))
plot(eq(cf.bd3[[1]],cf.bd3[[2]], x), type='l', x=x,lwd=2.5, ylim=c(0,4),yaxt="n",
     xlab='', ylab='', col="orange")
points(subset(rates, group == "diatoms")$temperature, 
       subset(rates, group == "diatoms")$r, col = alpha("#026cb1", 0.1), pch=16)
points(subset(rates, group == "coccolithophores")$temperature, 
       subset(rates, group == "coccolithophores")$r, col = alpha("orange", 0.1), pch=17)
points(subset(rates, group == "cyanobacteria")$temperature, 
       subset(rates, group == "cyanobacteria")$r, col = alpha("#ec3a25", 0.1), pch=16)
points(subset(rates, group == "diazotrophs")$temperature, 
       subset(rates, group == "diazotrophs")$r, col = alpha("darkorchid3", 0.1), pch=17)
points(subset(rates, group == "dinoflagellates")$temperature, 
       subset(rates, group == "dinoflagellates")$r, col = alpha("brown", 0.1), pch=18)
points(subset(rates, group == "greens")$temperature, 
       subset(rates, group == "greens")$r, col = alpha("#3ea127", 0.1), pch=18)
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[3]]),add=T,col="#ec3a25",lwd=2.5,lty=1) 
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[4]]),add=T,col="#026cb1",lwd=2.5,lty=1) 
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[5]]),add=T,col="darkorchid3",lwd=2.5,lty=1) # diazos
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[6]]),add=T,col="brown",lwd=2.5,lty=1) 
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[7]]),add=T,col="#3ea127",lwd=2.5,lty=1) 
text(0,3.75,quote(bold("(b)")* " Kremer"), line=1, adj=0.05, cex=1.2)
legend(-2.5, 3.65, legend=c(group[[1]], group[[2]],group[[3]],group[[4]], group[[5]],"green algae"),
       col=c("orange", "#ec3a25","#026cb1", "darkorchid3", "brown", "#3ea127"),
       lwd=2, lty=1, cex=1,box.lty=0,bg="transparent",seg.len=1)

# Anderson
par(mar=c(1,1,0.5,0.5))
plot(eq(-0.3,	0.035, x), type='l', x=x,lwd=2.5, ylim=c(0,4),yaxt="n",
     xlab='', ylab='', col="orange")
points(subset(rates, group == "diatoms")$temperature, 
       subset(rates, group == "diatoms")$r, col = alpha("#026cb1", 0.1), pch=16)
points(subset(rates, group == "coccolithophores")$temperature, 
       subset(rates, group == "coccolithophores")$r, col = alpha("orange", 0.1), pch=17)
points(subset(rates, group == "cyanobacteria")$temperature, 
       subset(rates, group == "cyanobacteria")$r, col = alpha("#ec3a25", 0.1), pch=16)
points(subset(rates, group == "diazotrophs")$temperature, 
       subset(rates, group == "diazotrophs")$r, col = alpha("darkorchid3", 0.1), pch=17)
points(subset(rates, group == "dinoflagellates")$temperature, 
       subset(rates, group == "dinoflagellates")$r, col = alpha("brown", 0.1), pch=18)
points(subset(rates, group == "greens")$temperature, 
       subset(rates, group == "greens")$r, col = alpha("#3ea127", 0.1), pch=18)
curve(eq(-1.661,	0.076, x),add=T,col="#ec3a25",lwd=2.5)
curve(eq(-0.226,	0.044, x),add=T,col="#026cb1",lwd=2.5)
curve(eq(-1.262,	0.0514, x),add=T,col="brown",lwd=2.5)
curve(eq(cf.bg[[1]],	cf.bg[[2]], x),add=T,col="#3ea127",lwd=2.5) # greens
curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[5]]),add=T,col="darkorchid3",lwd=2.5,lty=1)
text(0,3.75,quote(bold("(c)")* " Anderson"), line=1, adj=0.05, cex=1.2)

mtext(expression("Temperature (ºC)"), side = 1, outer = TRUE, line = 1.3, cex=0.8)
mtext(expression("Specific Growth Rate (d"^"-1" *")"), side = 2, outer = TRUE, line = 1.3, cex=0.8)
dev.off()
}


################################################################################
#### Figure S2: Fits for Kremer method
################################################################################

colors <- c("orange","#ec3a25", "#026cb1","darkorchid3","brown","#3ea127")
{pdf("figures/FigureS2.pdf", width = 7.2, height = 4.5)
  par(mfrow=c(3,2), mar=c(0,3.5, 3, 0.5))
  x =seq(min(subset(rates, group=="coccolithophores")$temperature),
         max(subset(rates, group=="coccolithophores")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "coccolithophores"),ylim=c(0,3),xlim=c(-2,40),xaxt="n",
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(a) Coccolithophores")), line=-1, adj=0.05, cex=1)
  axis(side=1,labels=F)
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x),min(x), max(x),add=T,col=colors[[1]],lwd=2.5,lty=1)
  y1_co <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x))
  y2_co <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x))
  polygon(c(x, rev(x)),c(y1_co, rev(y2_co)),col=alpha(colors[1],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "coccolithophores")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="coccolithophores")$isolate.code)), adj=c(0,0))
  
  par(mar=c(0,0.5,3,3.5))
  x =seq(min(subset(rates, group=="cyanobacteria")$temperature),
         max(subset(rates, group=="cyanobacteria")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "cyanobacteria"),ylim=c(0,3),xlim=c(-2,40),xaxt="n",yaxt="n",
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(b) Cyanobacteria")), line=-1, adj=0.05, cex=1)
  axis(side=1,labels=F)
  y1_cy <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x + (cf.bd3[3,1]-se[3,2]*1.96)))
  y2_cy <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x + (cf.bd3[3,1]+se[3,2]*1.96)))
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[3]]),min(x), max(x),add=T,col=colors[[2]],lwd=2.5,lty=1)
  polygon(c(x, rev(x)),c(y1_cy, rev(y2_cy)),col=alpha(colors[2],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "cyanobacteria")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="cyanobacteria")$isolate.code)), adj=c(0,0))
  
  par(mar=c(1.5, 3.5,1.5,0.5))
  x =seq(min(subset(rates, group=="diatoms")$temperature),
         max(subset(rates, group=="diatoms")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "diatoms"),ylim=c(0,3),xlim=c(-2,40),xaxt="n",
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(c) Diatoms")), line=-1, adj=0.05, cex=1)
  y1_d <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x + (cf.bd3[4,1]-se[4,2]*1.96)))
  y2_d <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x + (cf.bd3[4,1]+se[4,2]*1.96)))
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[4]]),min(x), max(x),add=T,col=colors[[3]],lwd=2.5,lty=1) #weighted by mean mass 
  polygon(c(x, rev(x)),c(y1_d, rev(y2_d)),col=alpha(colors[3],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "diatoms")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="diatoms")$isolate.code)), adj=c(0,0))
  
  par(mar=c(1.5, 0.5,1.5,3.5))
  x =seq(min(subset(rates, group=="diazotrophs")$temperature),
         max(subset(rates, group=="diazotrophs")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "diazotrophs"),ylim=c(0,3),xlim=c(-2,40),
       yaxt="n",xaxt="n",
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(d) Diazotrophs")), line=-1, adj=0.05, cex=1)
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[5]]),min(x), max(x),add=T,col=colors[[4]],lwd=2.5,lty=1) #weighted by mean mass 
  y1_di <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x +(cf.bd3[5,1]-se[5,2]*1.96)))
  y2_di <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x +(cf.bd3[5,1]+se[5,2]*1.96)))
  polygon(c(x, rev(x)),c(y1_di, rev(y2_di)),col=alpha(colors[4],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "diazotrophs")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="diazotrophs")$isolate.code)), adj=c(0,0))
  
  par(mar=c(3.5, 3.5,0,0.5))
  x =seq(min(subset(rates, group=="dinoflagellates")$temperature),
         max(subset(rates, group=="dinoflagellates")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "dinoflagellates"),ylim=c(0,3),xlim=c(-2,40),
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(e) Dinoflagellates")), line=-1, adj=0.05, cex=1)
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[6]]),min(x), max(x),add=T,col=colors[[5]],lwd=2.5,lty=1) #weighted by mean mas
  y1_dn <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x +(cf.bd3[6,1]-se[6,2]*1.96)))
  y2_dn <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x +(cf.bd3[6,1]+se[6,2]*1.96)))
  polygon(c(x, rev(x)),c(y1_dn, rev(y2_dn)),col=alpha(colors[5],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "dinoflagellates")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="dinoflagellates")$isolate.code)), adj=c(0,0))
  
  par(mar=c(3.5, 0.5,0,3.5))
  x =seq(min(subset(rates, group=="greens")$temperature),
         max(subset(rates, group=="greens")$temperature), by=0.1)
  plot(r~temperature,data=subset(rates, group == "greens"),ylim=c(0,3),xlim=c(-2,40),yaxt="n",
       xlab='',ylab='',pch=20, col=alpha("black", 0.4))
  title(main=expression(bold("(f) Green Algae")), line=-1, adj=0.05, cex=1)
  curve(exp(cf.bd3[[1]]+cf.bd3[[2]]*x+cf.bd3[[7]]),min(x), max(x),add=T,col=colors[[6]],lwd=2.5,lty=1) #weighted by mean mass 
  y1_g <- c(exp((cf.bd3[1,1]-se[1,2]*1.96)+(cf.bd3[2,1]-se[2,2]*1.96)*x +(cf.bd3[7,1]-se[7,2]*1.96)))
  y2_g <- c(exp((cf.bd3[1,1]+se[1,2]*1.96)+(cf.bd3[2,1]+se[2,2]*1.96)*x +(cf.bd3[7,1]+se[7,2]*1.96)))
  polygon(c(x, rev(x)),c(y1_g, rev(y2_g)),col=alpha(colors[6],alpha=0.2), border=FALSE)
  text(-1.6, 2.2, paste0("n=", length(unique(subset(isolates, group == "greens")$isolate.code))), adj=c(0,0))
  text(-1.6, 1.8, paste0("N=", length(subset(ratesg, group=="greens")$isolate.code)), adj=c(0,0))
  
  mtext(expression(bold("Temperature (ºC)")), side = 1, outer = TRUE, line = -1, cex=0.8)
  mtext(expression(bold("Specific Growth Rate (d"^"-1" *")")), side = 2, outer = TRUE, line = -1.5, cex=0.8)
  dev.off()}

################################################################################
#### Figure S3: Fits for Anderson method
################################################################################
{
colors <- c("orange","#ec3a25", "#026cb1","darkorchid3","brown","#3ea127")
cocco<-subset(rates, group =='coccolithophores')
cyano<-subset(rates, group =='cyanobacteria')
diatoms<-subset(rates, group =='diatoms')
dinos<-subset(rates, group =='dinoflagellates')

# Coefficient values are all taken from Anderson et al. 2021
pdf("figures/FigureS3.pdf", width = 7.2, height = 4.5)
x=seq(6,30, by=0.1)
y1_co <- c(exp(-0.2493228 + 0.04086702 * x))
y2_co <- c(exp(-0.3970879 + 0.03192307 * x))
par(mfrow = c(3,2), mar=c(0.5,0.5,0.4,0.4), oma=c(3.5,3.5,1,1), bg = "transparent")
plot(cocco$temperature, cocco$r, xlim=c(-2, 40), ylim=c(0,3), xaxt='n',xlab='', ylab='', pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(-0.3004739 + 0.03533314 * x),min(x), max(x),add=T,col=colors[[1]],lwd=2.5)
polygon(c(x, rev(x)),c(y1_co, rev(y2_co)),col=alpha(colors[[1]],alpha=0.2), border=FALSE)
title(main=expression(bold("(a) Coccolithophores")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(cocco$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="coccolithophores")$isolate.code)), adj=c(0,0))

x=seq(9,35, by=0.1)
y1_cy <- c(exp(-1.383528 + 0.08958329 * x))
y2_cy <- c(exp(-2.043423 + 0.06136036 * x))
plot(cyano$temperature, cyano$r, xlim=c(-2, 40), ylim=c(0,3), xaxt='n', yaxt='n',xlab='', ylab='', pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(-1.661365 + 0.0757542 * x),min(x), max(x),add=T,col=colors[2], lwd=2.5)
polygon(c(x, rev(x)),c(y1_cy, rev(y2_cy)),col=alpha(colors[2], alpha=0.2), border=FALSE)
title(main=expression(bold("(b) Cyanobacteria")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(cyano$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates,group=="cyanobacteria")$isolate.code)), adj=c(0,0))

tempd = seq(min(diatoms$temperature),
               max(diatoms$temperature), by=0.1)
y1 <- c(exp(-0.1837165 + 0.04833889 * tempd))
y2 <- c(exp(-0.2765283 + 0.04049121 * tempd))
plot(diatoms$temperature, diatoms$r, xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='', xaxt='n',pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(-0.226313 + 0.04382048 * x),min(x), max(x),add=T,col=colors[3],lwd=2.5)
polygon(c(tempd, rev(tempd)),c(y1, rev(y2)),col=alpha(colors[3],alpha=0.2), border=FALSE)
title(main=expression(bold("(c) Diatoms")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(diatoms$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="diatoms")$isolate.code)), adj=c(0,0))

x=seq(18, 34, by=0.1)
y1_diaz <- c(exp(-1.510912 + 0.2310916 * x))
y2_diaz <- c(exp(-6.546415 + 0.03031033 * x))
diazo.rates <- read.csv("data/diazo_growth_abbr.csv")
diazo.rates$group = "diazotrophs"
plot(diazo.rates$temperature, diazo.rates$r, xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='', yaxt='n',xaxt='n',pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(-4.445751+0.1494887*x),min(x), max(x),add=T,col=colors[4],lwd=2.5)
polygon(c(x, rev(x)),c(y1_diaz, rev(y2_diaz)),col=alpha(colors[4],alpha=0.2), border=FALSE)
title(main=expression(bold("(d) Diazotrophs")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=7"), adj=c(0,0))
text(-1.6, 1.9, paste0("N=144"), adj=c(0,0))

tempdi = seq(min(dinos$temperature),max(dinos$temperature), by=0.1)
x= tempdi
y1_df <- c(exp(-0.8729579 + 0.05962429 * tempdi))
y2_df <- c(exp(-1.358066 + 0.02777962 * tempdi))
plot(dinos$temperature, dinos$r, xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='',pch=20, col=alpha("black", 0.4))
curve(exp(-1.23558 + 0.05120331 *x),min(x), max(x),add=T,col=colors[5],lwd=2.5)
polygon(c(tempdi, rev(tempdi)),c(y1_df, rev(y2_df)),col=alpha(colors[5],alpha=0.2), border=FALSE)
title(main=expression(bold("(e) Dinoflagellates")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(dinos$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="dinoflagellates")$isolate.code)), adj=c(0,0))

# Green algae are new to this study
tempg=seq(min(ratesg$temperature),max(ratesg$temperature), by=0.1)
x=tempg
y1_g <- c(exp(ci_g[1,2]+ci_g[2,2]*tempg))
y2_g <- c(exp(ci_g[1,1]+ci_g[2,1]*tempg))
plot(subset(ratesg, group == "greens")$temperature, subset(ratesg, group == "greens")$r,  xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='', yaxt='n', pch=20, col=alpha("black", 0.4))
curve(exp(cf.bg[[1]]+cf.bg[[2]]*x),min(x), max(x),add=T,col=colors[6],lwd=2.5)
polygon(c(tempg, rev(tempg)),c(y1_g, rev(y2_g)),col=alpha(colors[6],alpha=0.2), border=FALSE)
title(main=expression(bold("(f) Green Algae")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(subset(isolates, group == "greens")$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(ratesg, group=="greens")$isolate.code)), adj=c(0,0))

mtext(expression(bold("Temperature (ºC)")), side = 1, outer = TRUE, line = 1.75)
mtext(expression(bold("Specific Growth Rate (d"^"-1" *")")), side = 2, outer = TRUE, line = 1.75)
dev.off()}

################################################################################
#### Figure S4: Allometrically scaled maximum growth curves
################################################################################
data <- read.csv("output/tableS1.csv")
target <- c("cyanobacteria", "greens", "coccolithophores", "diazotrophs", "diatoms", "dinoflagellates")
data$group = factor(data$group, levels=target)
data = data[order(data$group),]
data$num = ifelse(data$model == "Eppley", 1, 
                  ifelse(data$model == "Kremer", 2, 3))

temp = seq(-2, 40, by=0.1)
ipnum1 = c(1, 3, 5, 10, 15, 24) 
ipnum2 = c(2, 4, 9, 14, 23, 31)

# maximum growth rates either ascending or descending depending on PFT
data$b_coeff =  ifelse(data$group %in% c("cyanobacteria", "greens"), 0.28, -0.1)

# Numbers corresponding to the biovolume of each PFT phenotype
data$mumax_biovolnum = ifelse(data$group == "cyanobacteria", 2, 
                                ifelse(data$group == "greens", 4,
                                       ifelse(data$group == "coccolithophores", 5,
                                              ifelse(data$group == "diazotrophs", 10, 
                                                     ifelse(data$group == "diatoms", 15, 24)))))

# will need biovolume values (line 1155)
colors <- c("#ec3a25","#3ea127","orange", "darkorchid3","#026cb1","brown")
colors2 <- c("orange","#ec3a25", "#026cb1","darkorchid3","brown","#3ea127")
{pdf("figures/FigureS4.pdf", width = 6, height = 4)
  par(mfrow = c(2,3), mar=c(0.5,0.5,0.4,0.4), oma=c(2.5,2.5,1,1), bg = "transparent")
  for(k in 1:3){
    model = subset(data, num == k)
    target <- c("cyanobacteria", "greens", "coccolithophores", "diazotrophs", "diatoms", "dinoflagellates")
    model$group = factor(model$group, levels=target)
    model = model[order(model$group),]
    plot.new()
    plot.window(c(-2,30),c(0,3))
    box()
    for(i in 1:length(model$n)){
      curve(exp(model$a[i]+model$b[i]*x),add=T,col=colors[i], lwd=2)
    }
    if(k == 1){
      legend(13, 1.5, legend=c("coccolithophores","cyanobacteria",  "diatoms","diazotrophs",
                               "dinoflagellates","green algae"),
             col=c(colors2[[1]],colors2[[2]],colors2[[3]],colors2[[4]],colors2[[5]],colors2[[6]]),
             lwd=2, lty=1, cex=0.8,box.lty=0,bg="transparent",seg.len=1)
      axis(2, 0.5*(0:6), mgp=c(1,0.5,0), cex.axis=0.9)
      text(-1.6, 2.71, substitute(paste(bold((a)))), adj=c(0,0))
    }
    if(k == 2){
      text(-1.6, 2.71, substitute(paste(bold((b)))), adj=c(0,0))
    }
    if(k == 3){
      text(-1.6, 2.71, substitute(paste(bold((c)))), adj=c(0,0))
    }
    text(1.6, 2.75, paste0(model$model[1]), adj=c(0,0))
  }
  
  for(k in 4:6){
    model = subset(data, num == k-3)
    target <- c("cyanobacteria", "greens", "coccolithophores", "diazotrophs", "diatoms", "dinoflagellates")
    model$group = factor(model$group, levels=target)
    model = model[order(model$group),]
    mu20sa_vol = array()
    exposa <- model$b
    plot.new()
    plot.window(c(-2,30),c(0,3))
    if(k==4){
      axis(2, 0.5*(0:6), mgp=c(5,0.5,0), cex.axis=0.9)
    }
    axis(1, 10*(-2:30), mgp=c(5,0.5,0), cex.axis=0.9)
    box()
    for (i in 1:6){
      for(j in ipnum1[i]:ipnum2[i]){
        mu20sa_vol[j]=model$a_PCmax[i]*biovol[j]^model$b_coeff[i]
        x=temp
        curve((mu20sa_vol[j])*exp(exposa[i]*(x-20)),add=T,col=colors[i], lwd=2)
      }
    }
    mtext(expression("Temperature (ºC)"), side = 1, outer = TRUE, line = 1, cex=0.7)
    mtext(expression("Specific Growth Rate (d"^"-1" *")"), side = 2, outer = TRUE, line = 1, cex=0.7)
  }
  dev.off()
}

