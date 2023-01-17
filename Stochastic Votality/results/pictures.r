#-*- coding:utf-8 -*-
alpha_91_ori = read.csv("mse_all_sys_alpha91_ori.csv",header = F)
alpha_91_ori = as.matrix(alpha_91_ori)
boxplot(alpha_91_ori)
a = which.max(alpha_91_ori[,4])
alpha_91_ori_m = colMeans(alpha_91_ori[-a,])
ess = seq(0.1,1,by=0.1)
plot(ess,alpha_91_ori[6,],xlab = "ESS Threshold", ylab = "MSE")

alpha_91 = read.csv("mse_all_sys_alpha91.csv",header = F)
alpha_91 = as.matrix(alpha_91)
boxplot(alpha_91)
alpha_91_m = colMeans(alpha_91)
plot(ess, alpha_91[2,], xlab = "ESS Threshold", ylab = "MSE")

alpha_08 = read.csv("mse_all_alpha_08.csv", header=F)
alpha_08 = as.matrix(alpha_08)
boxplot(alpha_08)
alpha_08_m = colMeans(alpha_08)
plot(alpha_08_m)
for (i in 1:length(alpha_08[,1])) {
  plot(alpha_08[i,])
}
plot(ess,alpha_08[4,],xlab = "ESS Threshold", ylab = "MSE")

alpha_08_ori = read.csv("mse_all_alpha_08_ori.csv", header=F)
alpha_08_ori = as.matrix(alpha_08_ori)
boxplot(alpha_08_ori)
alpha_08_ori_m = colMeans(alpha_08_ori)
plot(alpha_08_ori_m)
for (i in 1:length(alpha_08_ori[,1])) {
  plot(alpha_08_ori[i,])
}
plot(ess, alpha_08_ori[6,], xlab = "ESS Threshold", ylab="MSE")
