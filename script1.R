rm(list=ls())

### Utils ###

library(timeDate)
library(timeSeries)
library(tseries)
library(zoo)
library(astsa)
library(urca)
library(lmtest)
library(fUnitRoots)
library(fBasics)
library(fImport)
library(stats)
library(fAssets)
library(forecast)
library(fPortfolio)
library(fRegression)
library(fTrading)
library(stabledist)
library(fGarch)
library(rugarch)
library(rmgarch)
library(FinTS)
library(ggplot2)
library(EnvStats)
library(pROC)
library(boot)
library(rvest)
# 
library(mfGARCH)

head(df_mfgarch)
head(df_financial)
tail(df_mfgarch$rv)
cor(v$return, v$dhousing, v$dindpro, v$nai, v$nfci)
v<-data.frame(df_mfgarch[,1:11])
head(v)
v[,c(1,11)]
cor(v[,c(2, 7:10)])
str(v)
v[v$date=="1990-01-02",]
v[v$date=="2017-12-29",]
v<-v[1:11644,]
head(v)
tail(v)

ggplot(data=v, aes(x=date, y=return)) +
  geom_line(color="steelblue") +
  scale_x_date(date_labels="%Y") +
  labs(title = "", x = "time", y = "index")

sp500_stat<-round(basicStats(v$return), 2)

xtable(round(basicStats(v$return),3))

plot1<-ggplot(data=v, aes(x=date, y=dhousing)) +
  geom_line(color="steelblue") +
  scale_x_date(date_labels="%Y") +
  labs(title = expression(Delta ~ "Housing"), x = "time", y = "index")

plot2<-ggplot(data=v, aes(x=date, y=dindpro)) +
  geom_line(color="steelblue") +
  scale_x_date(date_labels="%Y") +
  labs(title = expression(Delta ~ "IP"), x = "time", y = "index")

plot3<-ggplot(data=v, aes(x=date, y=nai)) +
  geom_line(color="steelblue") +
  scale_x_date(date_labels="%Y") +
  labs(title = "NAI", x = "time", y = "index")

plot4<-ggplot(data=v, aes(x=date, y=nfci)) +
  geom_line(color="steelblue") +
  scale_x_date(date_labels="%Y") +
  labs(title = "NFCI", x = "time", y = "index")

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)

# statistiche principali
sp500_stat<-round(basicStats(v$return), 2)
H_stat<-round(basicStats(v$dhousing),2)
IP_stat<-round(basicStats(v$dindpro),2)
nai_stat<-round(basicStats(v$nai),2)
nfci_stat<-round(basicStats(v$nfci),2)

prin_stat<-round(basicStats(v[,c(2,7:10)]),2)
xtable(prin_stat)

# stima del modello 
fitsp500_dhousing$tau
fitsp500_dhousing<-fit_mfgarch(data=v, y="return", x="dhousing", low.freq="year_month", K=36)
fitsp500_dindpro<-fit_mfgarch(data=v, y="return", x="dindpro", low.freq="year_month", K=36)
fitsp500_nai<-fit_mfgarch(data = v, y = "return", x = "nai", low.freq = "year_month", K = 36)
fitsp500_nfci<-fit_mfgarch(data = v, y = "return", x = "nfci", low.freq = "year_week", K = 52)

fitsp500_dhousing$broom.mgarch
summary(fitsp500_dhousing$std.err)

par_dhousing<-fitsp500_dhousing$broom.mgarch[c(2,3)]
par_dinpro<-fitsp500_dindpro$broom.mgarch[c(2,3)]
par_nai<-fitsp500_nai$broom.mgarch[c(2,3)]
par_nfci<-fitsp500_nfci$broom.mgarch[c(2,3)]

par_fit<-t(cbind(par_dhousing, par_dinpro, par_nai, par_nfci))
nomi_fit<-t(c(".","dhousing","dhousing", "dinpro", "dinpro", "NAI", "NAI", "NFCI", "NFCI"))
xtable(cbind(t(nomi_fit), round(par_fit, 3)))


# dataset per stime di dhousing 
tau_dhousing<-fitsp500_dhousing$tau
g_dhousing<-fitsp500_dhousing$g
v[v$date=="2000-01-03",]

length(tau_dhousing)
data_dhousing<-data.frame(date=v$date[7329:11644], rv_h=v$rv[7329:11644], 
                          tau_h=tau_dhousing[7329:11644], g_h=g_dhousing[7329:11644])
str(data_dhousing)

plot1<-ggplot () +
  geom_line(aes(x = data_dhousing$date, y = 12*sqrt(data_dhousing$rv_h), color = "steelblue"), size = 0.35, alpha = 1.1) +
  geom_line(aes(x = data_dhousing$date, y = 12*sqrt(data_dhousing$tau_h), color = "red"), size = 1, alpha = 0.5) +
  geom_line(aes(x = data_dhousing$date, y = 12*sqrt(data_dhousing$tau_h * data_dhousing$g_h), color = "green"), size = 0.8, alpha = 0.5) +
  scale_color_manual(name = " ", breaks = c("steelblue", "red", "green"), values = c("steelblue" = "steelblue", "red" = "red", "green"="green"),
                     labels = c("rv", "tau", "g*tau")) +
  labs(title = "1", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y")

# dataset per stime dinpro

tau_dinpro<-fitsp500_dindpro$tau
g_dinpro<-fitsp500_dindpro$g
data_dinpro<-data.frame(date=v$date[7329:11644], rv_p=v$rv[7329:11644],
                        tau_p=tau_dinpro[7329:11644], g_p=g_dinpro[7329:11644])

plot2<-ggplot () +
  geom_line(aes(x = data_dinpro$date, y = 12*sqrt(data_dinpro$rv_p), color = "steelblue"), size = 0.35, alpha = 1.1) +
  geom_line(aes(x = data_dinpro$date, y = 12*sqrt(data_dinpro$tau_p), color = "red"), size = 1, alpha = 0.5) +
  geom_line(aes(x = data_dinpro$date, y = 12*sqrt(data_dinpro$tau_p * data_dinpro$g_p), color = "green"), size = 0.8, alpha = 0.5) +
  scale_color_manual(name = " ", breaks = c("steelblue", "red", "green"), values = c("steelblue" = "steelblue", "red" = "red", "green"="green"),
                     labels = c("rv", "tau", "g*tau")) +
  labs(title = "2", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y")

# dataset for NAI estimation
tau_nai<-fitsp500_nai$tau
g_nai<-fitsp500_nai$g
data_nai<-data.frame(date=v$date[7329:11644], rv_nai=v$rv[7329:11644],
                     tau_n=tau_nai[7329:11644], g_n=g_nai[7329:11644])

plot3<-ggplot () +
  geom_line(aes(x = data_nai$date, y = 12*sqrt(data_nai$rv_nai), color = "steelblue"), size = 0.35, alpha = 1.1) +
  geom_line(aes(x = data_nai$date, y = 12*sqrt(data_nai$tau_n), color = "red"), size = 1, alpha = 0.5) +
  geom_line(aes(x = data_nai$date, y = 12*sqrt(data_nai$tau_n * data_nai$g_n), color = "green"), size = 0.8, alpha = 0.5) +
  scale_color_manual(name = " ", breaks = c("steelblue", "red", "green"), values = c("steelblue" = "steelblue", "red" = "red", "green"="green"),
                     labels = c("rv", "tau", "g*tau")) +
  labs(title = "3", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y")

# dataset for NFCI estimation
tau_nfci<-fitsp500_nfci$tau
g_nfci<-fitsp500_nfci$g
data_nfci<-data.frame(date=v$date[7329:11644], rv_f=v$rv[7329:11644],
                 tau_f=tau_nfci[7329:11644], g_f=g_nfci[7329:11644])
plot4<-ggplot () +
  geom_line(aes(x = data_nfci$date, y = 12*sqrt(data_nfci$rv_f), color = "steelblue"), size = 0.35, alpha = 1.1) +
  geom_line(aes(x = data_nfci$date, y = 12*sqrt(data_nfci$tau_f), color = "red"), size = 1, alpha = 0.5) +
  geom_line(aes(x = data_nfci$date, y = 12*sqrt(data_nfci$tau_f * data_nfci$g_f), color = "green"), size = 0.8, alpha = 0.5) +
  scale_color_manual(name = " ", breaks = c("steelblue", "red", "green"), values = c("steelblue" = "steelblue", "red" = "red", "green"="green"),
                     labels = c("rv", "tau", "g*tau")) +
  labs(title = "4", x = " ", y = " ") +
  scale_x_date(date_labels = "%Y")

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)


# plot of weights
w_dhousing<-data.frame(k=c(1:fitsp500_dhousing$K), w2=fitsp500_dhousing$est.weighting)
plot1<-ggplot(data=w_dhousing, aes(x=k, y=w2)) +
  geom_line(color="steelblue") +
  labs(title = expression(Delta ~ "Housing"), x = "k", y = "weight")
w_dinpro<-data.frame(k=c(1:fitsp500_dindpro$K), w2=fitsp500_dindpro$est.weighting)
plot2<-ggplot(data=w_dinpro, aes(x=k, y=w2)) +
  geom_line(color="steelblue") +
  labs(title = expression(Delta ~ "IP"), x = "k", y = "index")
w_nai<-data.frame(k=c(1:fitsp500_nai$K), w2=fitsp500_nai$est.weighting)
plot3<-ggplot(data=w_nai, aes(x=k, y=w2)) +
  geom_line(color="steelblue") +
  labs(title = expression("NAI"), x = "k", y = "index")
w_nfci<-data.frame(k=c(1:fitsp500_nfci$K), w2=fitsp500_nfci$est.weighting)
plot4<-ggplot(data=w_nfci, aes(x=k, y=w2)) +
  geom_line(color="steelblue") +
  labs(title = expression("NFCI"), x = "k", y = "index")

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)

llh_val<-c(fitsp500_dhousing$llh, fitsp500_dindpro$llh, 
           fitsp500_nai$llh, fitsp500_nfci$llh)
bic_val<-c(fitsp500_dhousing$bic, fitsp500_dindpro$bic,
           fitsp500_nai$bic, fitsp500_nfci$bic)
vr_val<-c(fitsp500_dhousing$variance.ratio, fitsp500_dindpro$variance.ratio,
          fitsp500_nai$variance.ratio, fitsp500_nfci$variance.ratio)
val<-rbind(llh_val, bic_val, vr_val)
nome_val<-c("llh", "BIC", "VR")
nomi_fit1<-c("dhousing", "dinpro", "NAI", "NFCI")
val1<-cbind(nomi_fit1, round(t(val), 3))
xtable(val1)

llh_val_n<-c(fitsp500_dhousing$llh, fitsp500_dindpro$llh, 
           fitsp500_nai$llh, fitsp500_nfci$llh)
bic_val_n<-c(fitsp500_dhousing$bic, fitsp500_dindpro$bic,
           fitsp500_nai$bic, fitsp500_nfci$bic)
vr_val_n<-c(fitsp500_dhousing$variance.ratio, fitsp500_dindpro$variance.ratio,
          fitsp500_nai$variance.ratio, fitsp500_nfci$variance.ratio)
aic_val_n<-c(AIC_dhousing_n, AIC_dinpro_n, AIC_nai_n, AIC_nfci_n)

llh_val_t<-c(fit_prova$llh, fit_t_dinpro$llh,
             fit_t_nai$llh, fit_t_nfci$llh)
bic_val_t<-c(fit_prova$bic, fit_t_dinpro$bic,
             fit_t_nai$bic, fit_t_nfci$bic)
vr_val_t<-c(fit_prova$vr, fit_t_dinpro$vr,
             fit_t_nai$vr, fit_t_nfci$vr)
aic_val_t<-c(AIC_dhousing_t, AIC_dinpro_t, AIC_nai_t, AIC_nfci_t)

llh_val_GED<-c(fit_prova_GED$llh, fit_GED_dinpro$llh,
               fit_GED_nai$llh, fit_GED_nfci$llh)
bic_val_GED<-c(fit_prova_GED$bic, fit_GED_dinpro$bic,
               fit_GED_nai$bic, fit_GED_nfci$bic)
rv_val_GED<-c(fit_prova_GED$rv, fit_GED_dinpro$rv,
               fit_GED_nai$rv, fit_GED_nfci$rv)
aic_val_GED<-c(AIC_dhousing_GED, AIC_dinpro_GED, AIC_nai_GED, AIC_nfci_GED)

dis_n<-c("Normale", "Normale", "Normale", "Normale")
dis_t<-c("t-Student", "t-Student", "t-Student", "t-Student")
dis_GED<-c("GED", "GED", "GED", "GED")
mod<-c("dhousing", "dindpro", "nai", "nfci")

normv<-cbind(dis_n, mod, round(cbind(llh_val_n, bic_val_n, aic_val_n), 3))
tstud<-cbind(dis_t, mod, round(cbind(llh_val_t, bic_val_t, aic_val_t), 3))
gedv<-cbind(dis_GED, mod, round(cbind(llh_val_GED, bic_val_GED, aic_val_GED), 3))
val_fin<-rbind(normv, tstud, gedv)

xtable(val_fin)



# AIC estimation
AIC_dhousing_n<-2*7 - 2*fitsp500_dhousing$llh
AIC_dhousing_t<-2*8 - 2*fit_prova$llh
AIC_dhousing_GED<-2*8 - 2*fit_prova_GED$llh

AIC_dinpro_n<-2*7 - 2*fitsp500_dindpro$llh
AIC_dinpro_t<-2*8 - 2*fit_t_dinpro$llh
AIC_dinpro_GED<-2*8 - 2*fit_GED_dinpro$llh

AIC_nai_n<-2*7 - 2*fitsp500_nai$llh
AIC_nai_t<-2*8 - 2*fit_t_nai$llh
AIC_nai_GED<-2*8 - 2*fit_GED_nai$llh

AIC_nfci_n<-2*7 - 2*fitsp500_nfci$llh
AIC_nfci_t<-2*8 - 2*fit_t_nfci$llh
AIC_nfci_GED<-2*8 - 2*fit_GED_nfci$llh



# Validation
qqnormPlot(v$return)

# Normality test
ksnormTest(v$return)
jarqueberaTest(v$return)
dagoTest(v$return)



fit.norm<-nFit(v$return)
fit.std<-stdFit(v$return)
fit.ged<-gedFit(v$return)

N<-length(v$return)
ggplot(data=v, aes(x=sort(v$return))) + 
  geom_point(aes(x=sort(v$return),y=(1:N/N)), col="steelblue", alpha = 1.4) +
  stat_function(fun=pnorm,args=list(mean=fit.norm@fit$estimate[1],sd=fit.norm@fit$estimate[2]),col=2,alpha=0.7,lwd=1)+
  labs(title="Confronto tra CDF empirica e teorica (std)", x="rendimenti", y="frequenze cumulate ") +
  theme(axis.text.x=element_text(angle=0, hjust=1, size=9), axis.text.y=element_text(size = 9), title = element_text(size=15))

plot1<-ggplot(data=v, aes(x=v$return)) + 
  geom_histogram(aes(y=..density..,),bins = 100, col="white", fill="steelblue") +
  geom_density(col=3,lwd=1) + theme(axis.text.x=element_text(angle=0, hjust=1, size=10), axis.text.y=element_text(size = 10), title = element_text(size=15)) +
  labs(title="", x="Rendimenti", y="Frequenza") +
  theme(axis.text.x=element_text(angle=0, hjust=1, size=9), axis.text.y=element_text(size = 9), title = element_text(size=10))+
  stat_function(fun=dstd,args=list(mean=fit.std$par[1],sd=fit.norm@fit$estimate[2]),col=2,alpha=0.7,lwd=1)
plot2<-ggplot(data = v, aes(sample = v$return)) + stat_qq(col="steelblue") + 
  stat_qq_line(col=2) + theme(axis.text.x = element_text(angle = 0, hjust = 1), title = element_text(size = 10)) + 
  labs(title = " ", x = "theoretical quantile", y = "sample quantile")
grid.arrange(plot1, plot2, nrow=1, ncol=2)

plot1<-ggplot(data=v, aes(x=v$return)) + 
  geom_histogram(aes(y=..density..,),bins = 100, col="white", fill="steelblue") +
  geom_density(col=3,lwd=1) + theme(axis.text.x=element_text(angle=0, hjust=1, size=10), axis.text.y=element_text(size = 10), title = element_text(size=10)) +
  labs(title="", x="rendimenti", y="Frequenza") +
  stat_function(fun=dstd,args=list(mean=fit.std$par[1],sd=fit.std$par[2], nu=fit.std$par[3]),col=2,alpha=0.7,lwd=1)

plot2<-ggplot(data=v, aes(x=v$return)) + 
  geom_histogram(aes(y=..density..,),bins = 100, col="white", fill="steelblue") +
  geom_density(col=3,lwd=1) + theme(axis.text.x=element_text(angle=0, hjust=1, size=10), axis.text.y=element_text(size = 10), title = element_text(size=10)) +
  labs(title="", x="rendimenti", y="Frequenza") +
  stat_function(fun=dged,args=list(mean=fit.ged$par[1],sd=fit.ged$par[2], nu=fit.ged$par[3]),col=2,alpha=0.7,lwd=1)

grid.arrange(plot1, plot2, nrow=1, ncol=2)


