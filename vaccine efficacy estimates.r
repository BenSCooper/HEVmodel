  # vaccine efficacy of Hep E vaccine  after 1,2, and 3 doses

#  Numbers based on Zhu et al Lancet 2010 (including key numbers in webappendix)
#  Numbers are are slightly revised from earlier version to exclude those in the reactogenicity subset who were not followed up after the first 6 months.
#  Checked on August 10th 2015.

# Efficacy based on 1 dose

# 3817 in the vaccine arm received one dose, 3874 in the placebo arm received one dose
# one HEV case in vaccine arm in this group, and one in the placebo arm
# Therefore no evidence of of vaccine efficacy for a single dose.
#
# Efficacy based on 2 doses

#n_v<- 3542  # number vaccinated with doses 1 and 2 only

# Use these numbers..this corresponds to number receiving the first 2 doses who were followed up between months 1.5 and 6 for HEV outcomes
n_p<- 52167-1329  # number getting placebo with doses 1 and 2 only (see fig 1. Excluces those who got 1 dose, dose 1 and 3 and reactogencity subset)
n_v<- 52235-1316  # number vaccinated (see fig 1 52235)

# Don't use numbers below here br - these represent total number receiving  at least 2 doses
 
#n_p<- 54972  # number getting placebo with doses 1 and 2 only
#n_v<- 54986  # number vaccinated

# Don't use numbers below for n_p and n_v - as the paper does not report HEV outcomes corresponding to these
# people alone (Ruby is writing to ask for these data)
#n_p<- 3504   # number getting placebo with doses 1 and 2 only
#n_v<- 3542  # number vaccinated

# Use these numbers..these correspond to HEV outcomes in those receiving the first 2 doses who were followed up between months 1.5 and 6 
# (note though that there may be longer follow-up of those receiving only 2-doses which would give us more data with which to estimate effect 
# of vaccine - this woluld follow-up for HEV outcomes of those who received only 2 doses beyond 6 months)
c_v<- 0  # HEV cases in those vaccinated with doses 1 and 2 only
c_p<- 5  # HEV cases in those getting placebo with doses 1 and 2 only

ARV<-c_v/n_v # attack rate vaccinated
ARU<-c_p/n_p # attack rate unvaccinated

VE<-100*(ARU-ARV)/ARU

# various ways to calculate CIs for this. We use Bayesain approach

#  Let p_v be the prob of being infected in the vaccine group
#  Let p_p be the prob of being infected in the placebo group
# our estimate of VE is 100(p_p - p_v)/p_p

# If we give the p_v and p_p beta(1,1) priors were can use beta-binomial conjugacy to estimate posterior distributions of both
#  If the prior for the probability of success (in this case reporting of a case) is beta(a, b), 
#  and we have n trials with x successes, then the posterior for the probability of success will 
# be beta(a+x, b+n-x). 
# so posterior for p_v will be beta(1+ c_v, 1+ n_v-c_v)  - so we create a sample form this - similarly for p_p
p_v.post.sample<-rbeta(100000, 1+ c_v, 1+ n_v-c_v)
p_p.post.sample<-rbeta(100000,1+ c_p, 1+ n_p-c_p)

VE2.post.sample<-100*(p_p.post.sample -p_v.post.sample)/p_p.post.sample

VE2.mean<-mean(VE2.post.sample)
VE2.qs<-quantile(VE2.post.sample,c(.025,.5,.975))
# gives a mean of 60.5, VE calc is 80.2, 95% CrIs are -34.2, 96.2

# Efficacy based on 3 doses

n_v<- 48693  # number vaccinated with doses 1  2  and 3 
n_p<- 48663  # number getting placebo with doses 1  2  and 3
c_v<- 0  # HEV cases in those vaccinated with doses 1  2  and 3
c_p<- 15  # HEV cases in those getting placebo with doses 1  2  and 3

ARV<-c_v/n_v # attack rate vaccinated
ARU<-c_p/n_p # attack rate unvaccinated

VE<-100*(ARU-ARV)/ARU

p_v.post.sample<-rbeta(100000, 1+ c_v, 1+ n_v-c_v)
p_p.post.sample<-rbeta(100000,1+ c_p, 1+ n_p-c_p)

VE3.post.sample<-100*(p_p.post.sample -p_v.post.sample)/p_p.post.sample

VE3.mean<-mean(VE3.post.sample)
VE3.qs<-quantile(VE3.post.sample,c(.025,.5,.975))

pdf("Histogram of two and three dose vaccine efficacy.pdf")
par(mfrow=c(2,1))
hist(VE2.post.sample[VE2.post.sample>0],breaks=c(0,1:100),freq=FALSE, main = "Posterior distribution of 2-dose vaccine efficacy", xlab= "Vaccine efficacy",ylim=c(0,.15))
hist(VE3.post.sample[VE3.post.sample>0],breaks=c(0,1:100),freq=FALSE, main = "Posterior distribution of 3-dose vaccine efficacy", xlab= "Vaccine efficacy",ylim=c(0,.15))
dev.off()


pdf("Posterior density of two and three dose vaccine efficacy.pdf")
plot(density(VE3.post.sample,bw="nrd",adjust=1),type='l', xlab="Vaccine efficacy", ylab="Density",main="",xlim=c(0,100),ylim=c(0,1.0))
par(new=T)
plot(density(VE2.post.sample,bw="nrd",adjust=1),type='l', xlab="Vaccine efficacy", ylab="Density",main="",xlim=c(0,100),ylim=c(0,1.0),col="blue",)
legend("topleft", legend=c("2 doses", "3 doses"),col=c("blue","black"),lty=1)

VE3subset.post.sample<-100*(p_p.post.sample -p_v.post.sample)/p_p.post.sample
VE3subset.mean<-mean(VE2subset.post.sample)
VE3subset.qs<-quantile(VE2subset.post.sample,c(.025,.5,.975))

dev.off()




# the above gives 
#  2.5%      50%    97.5% 
#  74.21090 95.61477 99.84041 
# which is very close reported CIs in the paper (72.1, 100)

# now add a line for what Zhu call the first two doses subset 
n_v<- 54986  # number vaccinated  in "first doses subset"
n_p<- 54973  # number getting placebo  in "first doses subset"
c_v<- 0  # HEV cases in those vaccinated in "first doses subset"
c_p<- 5  # HEV cases in those getting placebo  in "first doses subset"

ARV<-c_v/n_v # attack rate vaccinated
ARU<-c_p/n_p # attack rate unvaccinated

VE<-100*(ARU-ARV)/ARU

p_v.post.sample<-rbeta(100000, 1+ c_v, 1+ n_v-c_v)
p_p.post.sample<-rbeta(100000,1+ c_p, 1+ n_p-c_p)

VE2subset.post.sample<-100*(p_p.post.sample -p_v.post.sample)/p_p.post.sample
VE2subset.mean<-mean(VE2subset.post.sample)
VE2subset.qs<-quantile(VE2subset.post.sample,c(.025,.5,.975))
#> VE2subset.qs
#2.5%      50%    97.5% 
# 14.67174 87.66518 99.56908 




pdf("Posterior density of two and three dose vaccine efficacy.pdf")
plot(density(VE3.post.sample),type='l', xlab="Vaccine efficacy", ylab="Density",main="",xlim=c(-50,100),ylim=c(0,0.15))
par(new=T)
plot(density(VE2.post.sample),type='l', xlab="Vaccine efficacy", ylab="Density",main="",xlim=c(-50,100),ylim=c(0,0.15),col="blue",)
par(new=T)
plot(density(VE2subset.post.sample),type='l', xlab="Vaccine efficacy", ylab="Density",main="",xlim=c(-50,100),ylim=c(0,0.15),col="green",)
legend("topleft", legend=c("Exactly 2 doses","First 2 doses subset", "3 doses"),col=c("blue","green", "black"),lty=1)

dev.off()



