rm(list=ls())

# Estimate MIDAS-GARCH models
fit_n_dhousing<-fit_mfgarch(data = v, y = "return", x = "dhousing", K = 36, low.freq = "year_month")
fit_t_dhousing<-fit_mfgarch_t(data = v, y = "return", x = "dhousing", K = 36, low.freq = "year_month")
fit_GED_dhousing<-fit_mfgarch_GED(data = v, y = "return", x = "dhousing", K = 36, low.freq = "year_month")

## Estimate MIDAS-GARCH with three different models
fit_n_dindpro<-fit_mfgarch(data = v, y = "return", x = "dindpro", K = 36, low.freq = "year_month")
fit_n_nai<-fit_mfgarch(data = v, y = "return", x = "nai", K = 36, low.freq = "year_month")
fit_n_nfci<-fit_mfgarch(data = v, y = "return", x = "nfci", K = 52, low.freq = "year_week")
n<-length(v$return)
alfa<-c(0.01,0.05)
alfa.n<-qnorm(1-alfa)
alfa.t<-qstd(1-alfa, nu = fit_t_dhousing$par[8])
alfa.GED<-qged(1-alfa, nu = fit_GED_dhousing$par[8])

# Estimate Normal VaR with dhousing

# VaR01
Var1.n<-fit_n_dhousing$par[1]-
        sqrt(fit_n_dhousing$df.fitted$g*fit_n_dhousing$df.fitted$tau)*alfa.n[1]

# VaR05
Var2.n<-fit_n_dhousing$par[1]-
        sqrt(fit_n_dhousing$df.fitted$g*fit_n_dhousing$df.fitted$tau)*alfa.n[2]

# data.frame for ggplot visualization
Var_data.n<-data.frame(v$date, v$return, Var1.n, Var2.n)

# Estimate t-Student VaR with dhousing

# VaR01
Var1.t<-fit_t_dhousing$par[1]-
  sqrt(fit_t_dhousing$df.fitted$g*fit_t_dhousing$df.fitted$tau)*alfa.t[1]

# VaR05
Var2.t<-fit_t_dhousing$par[1]-
  sqrt(fit_t_dhousing$df.fitted$g*fit_t_dhousing$df.fitted$tau)*alfa.t[2]

# data.frame for ggplot visualization
Var_data.t<-data.frame(v$date, v$return, Var1.t, Var2.t)

# Estimate GED VaR with dhousing

# VaR01
Var1.GED<-fit_GED_dhousing$par[1]-
  sqrt(fit_GED_dhousing$df.fitted$g*fit_GED_dhousing$df.fitted$tau)*alfa.GED[1]

# VaR05
Var2.GED<-fit_GED_dhousing$par[1]-
  sqrt(fit_GED_dhousing$df.fitted$g*fit_GED_dhousing$df.fitted$tau)*alfa.GED[2]

# data.frame for ggplot visualization
Var_data.GED<-data.frame(v$date, v$return, Var1.GED, Var2.GED)

#
ggplot()+
  geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$v.return[757:n], color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$Var1.GED[757:n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$Var2.GED[757:n], color = "yellow"), alpha = 1.4 , size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 10), axis.text.y=element_text(size = 10), title=element_text(size=10)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow"),
                     labels=c("log--returns ", "VaR01", "VaR05")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(x = "time", y = "index") +
  scale_x_date(date_labels = "%Y")

#
ggplot()+
  geom_line(aes(x=Var_data.t$v.date[757:n],y=Var_data.t$v.return[757:n], color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=Var_data.t$v.date[757:n],y=Var_data.t$Var1.t[757:n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=Var_data.t$v.date[757:n],y=Var_data.t$Var2.t[757:n], color = "yellow"), alpha = 1.4 , size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 10), axis.text.y=element_text(size = 10), title=element_text(size=10)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow"),
                     labels=c("log--returns ", "VaR01", "VaR05")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(x = "time", y = "index") +
  scale_x_date(date_labels = "%Y")

# BackTest
BackTest.n.01<-BacktestVaR(data = Var_data.n$v.return[757:n], 
                           VaR = Var_data.n$Var1.n[757:n], alpha = alfa[1], Lags = 4)
BackTest.t.01<-BacktestVaR(data = Var_data.t$v.return[757:n], 
                           VaR = Var_data.t$Var1.t[757:n], alpha = alfa[1], Lags = 4)
BackTest.GED.01<-BacktestVaR(data = Var_data.GED$v.return[757:n], 
                           VaR = Var_data.GED$Var1.GED[757:n], alpha = alfa[1], Lags = 4)
BackTest.n.01$AE
BackTest.t.01$AE
BackTest.GED.01$AE

BackTest.n.05<-BacktestVaR(data = Var_data.n$v.return[757:n], 
                           VaR = Var_data.n$Var2.n[757:n], alpha = alfa[2], Lags = 4)
BackTest.t.05<-BacktestVaR(data = Var_data.t$v.return[757:n], 
                           VaR = Var_data.t$Var2.t[757:n], alpha = alfa[2], Lags = 4)
BackTest.GED.05<-BacktestVaR(data = Var_data.GED$v.return[757:n], 
                             VaR = Var_data.GED$Var2.GED[757:n], alpha = alfa[2], Lags = 4)
BackTest.n.05$AE
BackTest.t.05$AE
BackTest.GED.05$AE


ris_BackTest.n.01<-c(BackTest.n.01$AE, BackTest.n.01$AD[1], BackTest.n.01$AD[2], 
                     sum(Var_data.n$v.return[757:n]<Var_data.n$Var1.n[757:n]), 
                     BackTest.n.01$LRuc, BackTest.n.01$LRcc)
ris_BackTest.n.05<-c(BackTest.n.05$AE, BackTest.n.05$AD[1], BackTest.n.05$AD[2], 
                     sum(Var_data.n$v.return[757:n]<Var_data.n$Var1.n[757:n]), 
                     BackTest.n.05$LRuc, BackTest.n.05$LRcc)

ris_BackTest.t.01<-c(BackTest.t.01$AE, BackTest.t.01$AD[1], BackTest.t.01$AD[2], 
                     sum(Var_data.t$v.return[757:n]<Var_data.t$Var1.t[757:n]), 
                     BackTest.t.01$LRuc, BackTest.t.01$LRcc)
ris_BackTest.t.05<-c(BackTest.t.05$AE, BackTest.t.05$AD[1], BackTest.t.05$AD[2], 
                     sum(Var_data.t$v.return[757:n]<Var_data.t$Var1.t[757:n]), 
                     BackTest.t.05$LRuc, BackTest.t.05$LRcc)

ris_BackTest.GED.01<-c(BackTest.GED.01$AE, BackTest.GED.01$AD[1], BackTest.GED.01$AD[2], 
                     sum(Var_data.GED$v.return[757:n]<Var_data.GED$Var1.n[757:n]), 
                     BackTest.GED.01$LRuc, BackTest.GED.01$LRcc)
ris_BackTest.GED.05<-c(BackTest.GED.05$AE, BackTest.GED.05$AD[1], BackTest.GED.05$AD[2], 
                     sum(Var_data.GED$v.return[757:n]<Var_data.GED$Var1.GED[757:n]), 
                     BackTest.GED.05$LRuc, BackTest.GED.05$LRcc)

ris_back<-rbind(ris_BackTest.n.01, ris_BackTest.n.05, ris_BackTest.t.01, ris_BackTest.t.05,
                ris_BackTest.GED.01, ris_BackTest.GED.05)
# xtable(round(ris_back, 4))

# 
ggplot()+
  geom_line(aes(x=Var_data.n$v.date[757:n],y=Var_data.n$v.return[757:n], color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=Var_data.n$v.date[757:n],y=Var_data.n$Var1.g.n[757:n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=Var_data.n$v.date[757:n],y=Var_data.n$Var2.g.n[757:n], color = "yellow"), alpha = 1.4 , size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 16), axis.text.y=element_text(size = 16), title=element_text(size=17)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow"),
                     labels=c("log--returns ", "VaR01", "VaR05")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(x = "index", y = "")
scale_x_date(date_labels = "%Y")


# Estimate t-Student VaR for dhousing
Var1.g<-fit_t_dhousing$par[1]-
        sqrt(fit_t_dhousing$g*fit_t_dhousing$tau)*alfa.q[1]

Var2.g<-fit_t_dhousing$par[1]-
  sqrt(fit_t_dhousing$df.fitted$g*fit_t_dhousing$df.fitted$tau)*alfa.q[2]

Var_data<-data.frame(v$date, v$return, Var1.g, Var2.g)

ggplot()+
  geom_line(aes(x=Var_data$v.date[757:n],y=Var_data$v.return[757:n], color="blue"),size=0.15,alpha=0.3) +
  geom_line(aes(x=Var_data$v.date[757:n],y=Var_data$Var1.g[757:n], color="red"), alpha=1.4, size=0.35) +
  geom_line(aes(x=Var_data$v.date[757:n],y=Var_data$Var2.g[757:n], color = "yellow"), alpha = 1.4 , size = 0.35) +
  theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 16), axis.text.y=element_text(size = 16), title=element_text(size=17)) +
  scale_color_manual(name="Series", breaks=c("blue", "red", "yellow"),
                     values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow"),
                     labels=c("log--returns ", "VaR01", "VaR05")) +
  theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
  labs(x = "index", y = "")+
  scale_x_date(date_labels = "%Y")


# Estimate GED VaR for dhousing
Var1.g.GED<-fit_GED_dhousing$par[1]-
    sqrt(fit_GED_dhousing$g*fit_GED_dhousing$tau)*alfa.q.GED[1]
  
Var2.g.GED<-fit_GED_dhousing$par[1]-
    sqrt(fit_GED_dhousing$df.fitted$g*fit_GED_dhousing$df.fitted$tau)*alfa.q.GED[2]
  
Var_data.GED<-data.frame(v$date, v$return, Var1.g.GED, Var2.g.GED)
  
ggplot()+
    geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$v.return[757:n], color="blue"),size=0.15,alpha=0.3) +
    geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$Var1.g.GED[757:n], color="red"), alpha=1.4, size=0.35) +
    geom_line(aes(x=Var_data.GED$v.date[757:n],y=Var_data.GED$Var2.g.GED[757:n], color = "yellow"), alpha = 1.4 , size = 0.35) +
    theme(axis.text.x=element_text(angle = 0, hjust = 1, size = 16), axis.text.y=element_text(size = 16), title=element_text(size=17)) +
    scale_color_manual(name="Series", breaks=c("blue", "red", "yellow"),
                       values=c("blue" = "steelblue", "red" = "red", "yellow"= "yellow"),
                       labels=c("log--returns ", "VaR01", "VaR05")) +
    theme(legend.title=element_text(color = "black", size = 14),legend.text=element_text(color="black", size=14)) +
    labs(x = "index", y = "") +
    scale_x_date(date_labels = "%Y")

# Backtest
library(GAS)
BackTest.n.01<-BacktestVaR(data = Var_data$v.return[757:n], 
                      VaR = Var_data.n$Var1.g.n[757:n], alpha = alpha[1], Lags = 4)
BackTest.t.01<-BacktestVaR(data = Var_data$v.return[757:n], 
                           VaR = Var_data$Var1.g[757:n], alpha = alpha[1], Lags = 4)
BackTest.GED.01<-BacktestVaR(data = Var_data.GED$v.return[757:n], 
                           VaR = Var_data.GED$Var1.g.GED[757:n], alpha = alpha[1], Lags = 4)

BackTest.GED.01$AE
BackTest.n.01$AE
BackTest.t.01$AE


