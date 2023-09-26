rm(list = ls())

# datasets
d_apple<-read.csv2("AAPL.csv", header = T, sep = ",")
d_microsoft<-read.csv2("MSFT.csv", header = T, sep = ",")

d_pro<-read.csv("indpro.csv", header = T, sep = ",")
indpro<-diff(log(d_pro$INDPRO))
ind_pro<-rep(indpro, num)
ind_pro<-ind_pro*100
length(ind_pro)

d_unrate<-read.csv("unrate.csv", header = T, sep = ",")
un_rate<-diff(log(d_unrate$UNRATE))
unrate<-rep(un_rate, num)
unrate<-unrate*100
length(unrate)


length(d_pro$INDPRO_20201016)
head(d_apple)
head(d_microsoft)
str(d_apple)
d_apple$Date<-as.Date(d_apple$Date)
d_apple$Adj.Close<-as.numeric(d_apple$Adj.Close)
d_microsoft$Date<-as.Date(d_microsoft$Date)
d_microsoft$Adj.Close<-as.numeric(d_microsoft$Adj.Close)

da_apple<-data.frame(d_apple[,c(1,6)])
da_microsoft<-data.frame(d_microsoft[,c(1,6)])
head(da_apple)
str(da_apple)

# estimation of log-return for apple
return_a<-diff(log(da_apple$Adj.Close))
return_a100<-diff(log(da_apple$Adj.Close))*100
dat_apple<-data.frame(date=da_apple[-1,1], ret_a=return_a, ret_a100=return_a100)
head(dat_apple)
str(dat_apple)

# estimation of log-return of microsoft
return_m<-diff(log(da_microsoft$Adj.Close))
return_m100<-diff(log(da_microsoft$Adj.Close))*100
dat_microsoft<-data.frame(date=da_microsoft[-1,1], ret_m=return_m, ret_m100=return_m100)
head(dat_microsoft)
# ggplot 
plot1<-ggplot(data=dat_apple, aes(x=date, y=ret_a)) +
  geom_line(color="Indianred2") +
  scale_x_date(date_labels="%Y") +
  labs(title = "", x = "time", y = "index")
plot2<-ggplot(data=dat_microsoft, aes(x=date, y=ret_m)) +
  geom_line(color="royalblue4") +
  scale_x_date(date_labels="%Y") +
  labs(title = "", x = "time", y = "index")
grid.arrange(plot1, plot2, nrow=1, ncol=2)

basicStats(dat_apple$ret_a)
basicStats(dat_microsoft$ret_m)
cbind(basicStats(dat_apple$ret_a),
      basicStats(dat_microsoft$ret_m))
xtable(cbind(basicStats(dat_apple$ret_a),
             basicStats(dat_microsoft$ret_m)))


ind_unrate<-read.csv("un_rate.csv", header=T)
head(ind_unrate)
str(ind_unrate)
ind_unrate$DATE<-as.Date(ind_unrate$DATE)
ind_unrate<-data.frame(month=ind_unrate$DATE , unrate=ind_unrate$UNRATE_20201002)
head(ind_unrate)
str(ind_unrate)
n<-length(dat_apple$date)
l<-length(ind_unrate$month)


date_month<-format(dat_apple$date, "%m")
date_year<-format(dat_apple$date, "%Y")

year_month<-format(dat_apple$date, "%Y-%m-01")
str(year_month)
year_month<-as.Date(year_month)

num<-table(year_month)

length(num)

ind_un<-rep(ind_unrate$unrate, num)

length(ind_un)

data_apple<-data.frame(dat_apple, year_month, unrate=ind_un)
head(data_apple)

apple_fit_n<-fit_mfgarch(data = data_apple, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_t<-fit_mfgarch_t(data = data_apple, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_GED<-fit_mfgarch_GED(data = data_apple, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_n$bic
apple_fit_t$bic
apple_fit_GED$bic

tail(data_apple)

dati<-data.frame(date=data_apple$date, ret_a=data_apple$ret_a, ret_a100=data_apple$ret_a100,
                 ret_m=dat_microsoft$ret_m, ret_m100=dat_microsoft$ret_m100,
                 year_month, unrate, ind_pro)
head(dati)

apple_fit_n<-fit_mfgarch(data = dati, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_t<-fit_mfgarch_t(data = data_apple, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_GED<-fit_mfgarch_GED(data = data_apple, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
apple_fit_n$bic
apple_fit_t$bic
apple_fit_GED$bic

str(dati)
basicStats(dati[,c(3,5,7,8)])
xtable(basicStats(dati[,c(3,5,7,8)]))

jarque.bera.test(dati$ret_m100)
plot1<-ggplot(data = dati, aes(sample = ret_a)) + stat_qq(col="Indianred2") + 
  stat_qq_line(col="royalblue4") + theme(axis.text.x = element_text(angle = 0, hjust = 1), title = element_text(size = 10)) + 
  labs(title = "", x = "theoretical quantile", y = "sample quantile")
plot2<-ggplot(data = dati, aes(sample = ret_m)) + stat_qq(col="royalblue4") + 
  stat_qq_line(col=2) + theme(axis.text.x = element_text(angle = 0, hjust = 1), title = element_text(size = 10)) + 
  labs(title = "", x = "theoretical quantile", y = "sample quantile")
grid.arrange(plot1, plot2, nrow=1, ncol=2)

### MODELS ESTIMATION

# unrate apple
fit_apunr_n<-fit_mfgarch(data = dati, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
fit_apunr_t<-fit_mfgarch_t(data = dati, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
fit_apunr_GED<-fit_mfgarch_GED(data = dati, y = "ret_a100", x = "unrate", K = 36, low.freq = "year_month")
# unrate microsoft
fit_miunr_n<-fit_mfgarch(data = dati, y = "ret_m100", x = "unrate", K = 36, low.freq = "year_month")
fit_miunr_t<-fit_mfgarch_t(data = dati, y = "ret_m100", x = "unrate", K = 36, low.freq = "year_month")
fit_miunr_GED<-fit_mfgarch_GED(data = dati, y = "ret_m100", x = "unrate", K = 36, low.freq = "year_month")

# ind pro apple 
fit_appro_n<-fit_mfgarch(data = dati, y = "ret_a100", x = "ind_pro", K = 36, low.freq = "year_month")
fit_appro_t<-fit_mfgarch_t(data = dati, y = "ret_a100", x = "ind_pro", K = 36, low.freq = "year_month")
fit_appro_GED<-fit_mfgarch_GED(data = dati, y = "ret_a100", x = "ind_pro", K = 36, low.freq = "year_month")

# ind pro microsoft 
fit_mipro_n<-fit_mfgarch(data = dati, y = "ret_m100", x = "ind_pro", K = 36, low.freq = "year_month")
fit_mipro_t<-fit_mfgarch_t(data = dati, y = "ret_m100", x = "ind_pro", K = 36, low.freq = "year_month")
fit_mipro_GED<-fit_mfgarch_GED(data = dati, y = "ret_m100", x = "ind_pro", K = 36, low.freq = "year_month")

par_aunr_n<-fit_apunr_n$broom.mgarch[2:3]
par_aunr_n<-rbind(par_aunr_n, c(0, 0))
par_aunr_t<-fit_apunr_t$broom.mgarch[2:3]
par_aunr_GED<-fit_apunr_GED$broom.mgarch[2:3]

par_aunr<-cbind(par_aunr_n, par_aunr_t, par_aunr_GED)
bic<-c(fit_apunr_n$bic, 0, fit_apunr_t$bic, 0, fit_apunr_GED$bic, 0)
llh_val<-c(fit_apunr_n$llh, 0, fit_apunr_t$llh, 0, fit_apunr_GED$llh, 0)

xtable(t(round(rbind(par_aunr, llh_val, bic), 3)))

par_miunr_n<-fit_miunr_n$broom.mgarch[2:3]
par_miunr_n<-rbind(par_miunr_n, c(0, 0))
par_miunr_t<-fit_miunr_t$broom.mgarch[2:3]
par_miunr_GED<-fit_miunr_GED$broom.mgarch[2:3]

par_miunr<-cbind(par_miunr_n, par_miunr_t, par_miunr_GED)
bic<-c(fit_miunr_n$bic, 0, fit_miunr_t$bic, 0, fit_miunr_GED$bic, 0)
llh_val<-c(fit_miunr_n$llh, 0, fit_miunr_t$llh, 0, fit_miunr_GED$llh, 0)
xtable(t(round(rbind(par_miunr, llh_val, bic), 3)))

par_apro_n<-fit_appro_n$broom.mgarch[2:3]
par_apro_n<-rbind(par_apro_n, c(0, 0))
par_apro_t<-fit_appro_t$broom.mgarch[2:3]
par_apro_GED<-fit_appro_GED$broom.mgarch[2:3]

par_miunr<-cbind(par_apro_n, par_apro_t, par_apro_GED)
bic<-c(fit_appro_n$bic, 0, fit_appro_t$bic, 0, fit_appro_GED$bic, 0)
llh_val<-c(fit_appro_n$llh, 0, fit_appro_t$llh, 0, fit_appro_GED$llh, 0)
xtable(t(round(rbind(par_miunr, llh_val, bic), 3)))



par_mipro_n<-fit_mipro_n$broom.mgarch[2:3]
par_mipro_n<-rbind(par_mipro_n, c(0, 0))
par_mipro_t<-fit_mipro_t$broom.mgarch[2:3]
par_mipro_GED<-fit_mipro_GED$broom.mgarch[2:3]

par_mipro<-cbind(par_mipro_n, par_mipro_t, par_mipro_GED)
bic<-c(fit_mipro_n$bic, 0, fit_mipro_t$bic, 0, fit_mipro_GED$bic, 0)
llh_val<-c(fit_mipro_n$llh, 0, fit_mipro_t$llh, 0, fit_mipro_GED$llh, 0)
xtable(t(round(rbind(par_mipro, llh_val, bic), 3)))





