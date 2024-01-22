library(survival)
# install.packages('dplyr')
rm(list = ls())
library(dplyr)
# setwd("C:/Users/Lingkai/Desktop/KM R")

# ===================================================== settings ==================================================
nwindows = 12
Path_data = "HB/sliding windows/"
Path_out = "HB/KMR/"
sliding_windows = seq(1, nwindows, 1)
for(window in sliding_windows){
  # ================= load data
  file = paste(Path_data, window, '.csv', sep="")
  data<-read.csv(file,header = TRUE)
  
  # ================= KM for W
  y1 = data$y1
  y2 = data$y2
  y3 = data$y3
  tags = data$tag
  
  stime = y1
  etime = y1+y2+y3
  observation = tags
  
  ID_delete = which(stime == etime)
  if(length(which(stime == etime)) > 0){
    stime = stime[-ID_delete]
    etime = etime[-ID_delete]
    observation = observation[-ID_delete]    
  } 
  observation[observation==5 | observation==6 | observation==8] = 0
  observation[observation==1 | observation==2 | observation==3 | observation==4 | observation==7] = 1
  
  fit <- survfit(Surv(stime, etime, observation) ~ 1, data=data)
  Xtime = matrix(fit$time)
  Ytime = matrix(1 - fit$surv)
  # plot(Xtime, Ytime)
  
  # cbind()  and  rbind()
  out <- cbind(Xtime, Ytime)
  df <- data.frame(out)
  file_out = paste(Path_out, window, '_W.csv', sep="")
  write.table(df, file_out, sep = ",", row.names = FALSE, col.names =FALSE)
  
  # ================= KM for W1
  term1 = matrix(data$tag)
  term2 = matrix(c(1,2,3,4,5,6))
  data1 = data[is.element(term1, term2),]
  
  y1 = data1$y1
  y2 = data1$y2
  y3 = data1$y3
  tags = data1$tag
  
  stime = y1
  etime = y1+y2
  observation = tags
  
  ID_delete = which(stime == etime)
  if(length(which(stime == etime)) > 0){
    stime = stime[-ID_delete]
    etime = etime[-ID_delete]
    observation = observation[-ID_delete]    
  } 
  observation[observation==3 | observation==4 | observation==5 | observation==6] = 0
  observation[observation==1 | observation==2] = 1
  
  fit <- survfit(Surv(stime, etime, observation) ~ 1, data=data1)
  Xtime = matrix(fit$time)
  Ytime = matrix(1 - fit$surv)

  out <- cbind(Xtime, Ytime)
  df <- data.frame(out)
  file_out = paste(Path_out, window, '_W1.csv', sep="")
  write.table(df, file_out, sep = ",", row.names = FALSE, col.names =FALSE)
  
  # ================= KM for W2
  term1 = matrix(data$tag)
  term2 = matrix(c(3,4,5,6,7,8))
  data2 = data[is.element(term1, term2),]
  
  y1 = data2$y1
  y2 = data2$y2
  y3 = data2$y3
  tags = data2$tag
  
  stime = y1+y2
  etime = y1+y2+y3
  observation = tags
  
  ID_delete = which(stime == etime)
  if(length(which(stime == etime)) > 0){
    stime = stime[-ID_delete]
    etime = etime[-ID_delete]
    observation = observation[-ID_delete]    
  } 
  observation[observation==5 | observation==6 | observation==8] = 0
  observation[observation==3 | observation==4 | observation==7] = 1
  
  fit <- survfit(Surv(stime, etime, observation) ~ 1)
  
  Xtime = matrix(fit$time)
  Ytime = matrix(1 - fit$surv)
  
  out <- cbind(Xtime, Ytime)
  df <- data.frame(out)
  file_out = paste(Path_out, window, '_W2.csv', sep="")
  write.table(df, file_out, sep = ",", row.names = FALSE, col.names =FALSE)
}
