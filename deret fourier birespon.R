library(ggplot2)
library(patchwork)

library(pracma)
library(readxl)
Data <- read_excel("D:/PUBLIKASI/DRAFT WORKSHOP Q1/data penelitian ketahanan pangan 2024.xlsx")
y1=as.matrix(Data[,2]) 
y2=as.matrix(Data[,3]) 
x1=as.matrix(Data[,4]) 
x2=as.matrix(Data[,5]) 
x3=as.matrix(Data[,6]) 
x4=as.matrix(Data[,7]) 

par(mfrow=c(2,2))

plot(x1,y1,
     pch=19,
     main="y1 vs x1",
     xlab="Agricultural Extension Workers",
     ylab="Rice Productivity")

plot(x2,y1,
     pch=19,
     main="y1 vs x2",
     xlab="Rainfall",
     ylab="Rice Productivity")

plot(x3,y1,
     pch=19,
     main="y1 vs x3",
     xlab="Rice Harvested Area",
     ylab="Rice Productivity")

plot(x4,y1,
     pch=19,
     main="y1 vs x4",
     xlab="Agricultural Enterprises",
     ylab="Rice Productivity")

par(mfrow=c(2,2))

plot(x1,y2,
     pch=19,
     main="y2 vs x1",
     xlab="Agricultural Extension Workers",
     ylab="Food Security Index")

plot(x2,y2,
     pch=19,
     main="y2 vs x2",
     xlab="Rainfall",
     ylab="Food Security Index")

plot(x3,y2,
     pch=19,
     main="y2 vs x3",
     xlab="Rice Harvested Area",
     ylab="Food Security Index")

plot(x4,y2,
     pch=19,
     main="y2 vs x4",
     xlab="Agricultural Enterprises",
     ylab="Food Security Index")
#==========================================
n<-length(y1)
korelasi<-function(Data)
{
  cat("\n Uji Korelasi Pearson \n")
  cat("\n=================================\n")
  alfa<-as.numeric(readline("\n\n Input nilai alfa :"))
  sy1y2<-cov(y1,y2) #covariansi
  sy1<-var(y1) #variansi y1
  sy2<-var(y2) #variansi y2
  korelasi<-sy1y2/sqrt(sy1*sy2) #kok nilainya beda dengan yang excel
  korelasi
  cat("\n koefisien korelasinya yaitu:", korelasi, "\n")
  t<-(korelasi*sqrt(n-2))/sqrt(1-(korelasi^2))
  v<-n-2
  ttabel<-qt(1-(alfa/2),v)
  cat("Hipotesis: \n")
  cat("Ho :rho = 0\n")
  cat("H1 :rho !=0\n")
  p_value<-round(2*pt(abs(t),v,lower.tail=FALSE),4)
  cat("\n=========================","\n")
  cat("nilai p value= ",p_value,"\n")
  cat("\n=========================","\n")
  cat("\n kesimpulan :\n")
  if(p_value<alfa)
  {
    cat("Tolak Ho\nada korelasi antara kedua variabel respon\n\n")
  }
  else
  {
    cat("Gagal Tolak Ho\n Tidak ada korelasi antara kedua variabel respon\n\n")
  }
}
korelasi(Data)

#Uji Bartlett Sphericity
#================================================================================
Rcorr=matrix(c(cor(y1,y1),cor(y1,y2),cor(y2,y1),cor(y2,y2)),2,2)
chihitung=-1*(n-1-((2*2+5)/6))*log(det(Rcorr))
{
  if (chihitung<=3.841)
  {
    cat("Ho gagal tolak sehingga antar variabel respon independen")}
  else
  {
    cat("Ho ditolak sehingga antar variabel respon dependen ")
  }
}

#============================================
N=length(y1)#Banyaknya Observasi
N=nrow(Data)
#matriks pembobot (matriks varians kovarians)
s1=matrix(0,N,N)
s2=matrix(0,N,N)
s3=matrix(cov(y1,y2),N,N)
s4=matrix(cov(y1,y2),N,N)

for (p in 1:N)
{
  s1[p,p]=var(y1)
  s2[p,p]=var(y2)
}
s11=cbind(s1,s3)
s22=cbind(s4,s2)

library(MASS)
W=rbind(s11,s22)

deretfourier<-function(y1, y2,x1, x2, x3, x4, K)
{
  y=rbind(y1,y2)
  N=length(y1)
  K=1 #osilasi
  a<-(4*(K+1))+1 #hanya ada 1 konstanta
  C<-matrix(0,N,a)
  hasil<-matrix(0,K,2)
  for (k in 1:K)
  {
    for(i in 1:N)
    {
      for(j in 1:k)
      {
        C[i,1]<-x1[i]
        C[i,2]=0.5
        C[i,2+j]=cos(j*x1[i])
        C[i,3+k]<-x2[i]
        C[i,3+k+j]<-cos(j*x2[i])
        C[i,4+(2*k)]<-x3[i]
        C[i,4+(2*k)+j]<-cos(j*x3[i])
        C[i,5+(3*k)]<-x4[i]
        C[i,5+(3*k)+j]<-cos(j*x4[i])
      }
    }
    D=matrix(0,N,a)
    D1=rbind(C,D)
    D2=rbind(D,C)
    Dfou=cbind(D1,D2)#matriks deret fourier
    Dfou
    I<-diag(2*N)
    beta<-Dfou%*%inv(t(Dfou)%*%W%*%Dfou)%*%t(Dfou)%*%W
    #parameter estimasi deret fourier
    beta1<-inv(t(Dfou)%*%W%*%Dfou)%*%t(Dfou)%*%W%*%y
    ytopi<-beta%*%y
    #ytopi<-Dfou%*%beta
    error<-y-ytopi
    df<-(((2*N)^-1)*sum(diag(I-beta)))^2
    MSE<-((2*N)^-1)*(t(error)%*%error)
    GCV<-MSE/df
    SSE=sum((y-ytopi)^2)
    SSR=sum((ytopi-mean(y))^2)
    SST=SSE+SSR
    #SST=sum((y-mean(y))^2)
    R2=((1-(SSE/SST))*100)
    
    #simpan hasil GCV
    hasil[k,1]<-k
    hasil[k,2]<-GCV
    
    # Tampilkan nilai parameter dan hasil regresi
    cat("Hasil untuk K =", k, "\n")
    cat("Parameter:\n", beta1, "\n")
    cat("GCV:", GCV, "\n")
    cat("R2:", R2, "%\n")
  }
  print(hasil)
  GCV2<-min(hasil[,2])
  s<-1
  repeat{
    if(hasil[s,2]==GCV2)
    {
      kOpt<-hasil[s,1]
      GCVOpt<-GCV2
      break
    }
    else s<-s+1
  }
  cat("nilai l optimal adalah \t",kOpt,"\n")
  cat("nilai GCV min \t",GCVOpt,"\n")
  
  # Kembalikan nilai parameter optimal
  return(list(beta1 = beta1, GCVOpt = GCVOpt, kOpt = kOpt, R2 = R2))
}
hasil=deretfourier(y1,y2,x1,x2,x3,x4,1)
hasil

#plot data
#PLOTTING
#plot y1
y1=y[1:24]
y1topi=ytopi[1:24]

# Menentukan rentang ylim berdasarkan nilai minimum dan maksimum dari y1 dan y1topi
ylim_range <- range(c(y1, y1topi))

# Buat plot dengan rentang ylim yang disesuaikan
plot(y1, main="y1 Actual vs. y1 Predictions", xlab="Observation", ylab="Rice productivity", 
     type="l", lwd=2, col="red", cex.main=1.5, cex.lab=1.2, cex.axis=1.1, ylim=ylim_range)

# Tambahkan prediksi dengan garis dan titik lebih besar
lines(y1topi, col="blue", type="o", lwd=2, pch=16, cex=1.2)

# Tambahkan grid untuk memperjelas visualisasi
grid()

# Tambahkan legenda dengan ukuran dan warna lebih tebal
legend(x=62,y=25, legend=c("Actual", "Predictions"), col=c("red", "blue"), 
       lwd=2, cex=0.75, pch=c(NA,16), pt.cex=1.2)

#plot y2
y2=y[25:48]
y2topi=ytopi[25:48]

# Menentukan rentang ylim berdasarkan nilai minimum dan maksimum dari y1 dan y1topi
ylim_range <- range(c(y2, y2topi))

# Buat plot dengan rentang ylim yang disesuaikan
plot(y2, main="y2 Actual vs. y2 Predictions", xlab="Observation", ylab="Food Security Index", 
     type="l", lwd=2, col="red", cex.main=1.5, cex.lab=1.2, cex.axis=1.1, ylim=ylim_range)

# Tambahkan prediksi dengan garis dan titik lebih besar
lines(y2topi, col="blue", type="o", lwd=2, pch=16, cex=1.2)

# Tambahkan grid untuk memperjelas visualisasi
grid()

# Tambahkan legenda dengan ukuran dan warna lebih tebal
legend(x=62,y=15, legend=c("Actual", "Predictions"), col=c("red", "blue"), 
       lwd=2, cex=0.75, pch=c(NA,16), pt.cex=1.2)


