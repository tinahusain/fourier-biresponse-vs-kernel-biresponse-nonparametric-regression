library(pracma)
Data<- read_excel("D:/ITH/PDP INTERNAL ITH/Pendanaan 2025/PDP 2025/Koding dan referensi/data/data gabungan.xlsx")
View(Data)
y1=as.matrix(Data[,11])
y2=as.matrix(Data[,3])
#hitung korelasi antar kedua variabel respon
cor(y1,y2)
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

y=rbind(y1,y2)
xk=as.matrix(Data[,c(4,5,10,11)]) #variabel kernel
n=length(y) #jumlah observasi (n=48)
nk=ncol(xk) #banyaknya variabel kernel
int.ker=50 #jumlah pembagi titik bandwith yang diinginkan

#~~penentuan bandwith~~
#untuk respon 1
bwy1=matrix(0,int.ker,nk)
for(i in 1:nk)
{
  bwy1[,i]=seq(0,(max(xk[,i])-min(xk[,i])),length.out=int.ker)
}
bwy1=as.matrix(bwy1[2:(int.ker-1),])
nbandy1=nrow(bwy1)

#untuk respon 2
bwy2=matrix(0,int.ker+2,nk)
for(i in 1:nk)
{
  bwy2[,i]=seq(0,(max(xk[,i])-min(xk[,i])),length.out=int.ker+2)
}
bwy2=as.matrix(bwy2[2:(int.ker-1),])
nbandy2=nrow(bwy2)
bw=cbind(bwy1,bwy2) #gabungkan bandwith respon 1 dan respon 2
nband=nrow(bw) #48
vk=blkdiag(xk,xk) #48x6
sk=ncol(vk)
#matrix
m1.nn=matrix(1, nrow=n, ncol=n) # matriks 1 nxn
m1.n1=matrix(1, nrow=n) # matriks 1 nx1
mi.nn=diag(1,n,n) # matriks identitas nxn
e=1
#Desain matriks V pada kernel

MSE=matrix(0,nband)
GCV=matrix(0,nband)
R2=matrix(0,nband) 
code=matrix(0,nband,sk)
for(j in 1:nband)
{
  sum.q.phi=0
for(k in 1:sk)
{
  v.diag=diag(vk[,k]) #akan berubah menjadi matriks diagonal ukuran 48x48
  V1=m1.nn%*%v.diag
  z1=(t(V1)-V1)/bw[j,k] #bw adalah bandwith pada kernel
  K=(1/sqrt(2*pi))*exp(-1/2*z1^2) #fungsi kernel gaussian
  K.z1=(1/bw[j,k])*K #fungsi kernel
  pembobot.penyebut=diag(c(1/n*K.z1%*%m1.n1))%*%m1.nn
  q.phi=1/n*K.z1/pembobot.penyebut
  
  #penimbang q(phi).1
  sum.q.phi=sum.q.phi+q.phi #nilai kernel pada setiap varaiabel
}
#penimbang kernel gabungan
q.phi=sum.q.phi/sk #nilaI kernel rata-rata
Vker=q.phi
MSE[e]=((2*n)^-1)*(t(y)%*%t(mi.nn-Vker)%*%(mi.nn-Vker)%*%y)
db=(((2*n)^-1)*sum(diag(mi.nn-Vker)))^2
GCV[e]=MSE[e]/db
code[e,]=c(bw[j,1],bw[j,2],bw[j,3],bw[j,4])
e=e+1
}
optimum=cbind(code,MSE,GCV)
GCVmin=optimum[order(optimum[,10]),] #mengurutkan nilai GCV minimum
GCVmin[1,]

band.opt=GCVmin[1,1:8]

#validasi nilai GCV terkecil
sum.q.phi=0
for (k in 1:sk)
{
  v.diag=diag(vk[,k]) #akan berubah menjadi matriks diagonal ukuran 48x48
  V1=m1.nn%*%v.diag
  z1=(t(V1)-V1)/band.opt[k] #bw adalah bandwith pada kernel
  K=(1/sqrt(2*pi))*exp(-1/2*z1^2) #fungsi kernel gaussian
  K.z1=(1/band.opt[k])*K #fungsi kernel
  pembobot.penyebut=diag(c(1/n*K.z1%*%m1.n1))%*%m1.nn
  q.phi=1/n*K.z1/pembobot.penyebut
  
  #penimbang q(phi).1
  sum.q.phi=sum.q.phi+q.phi #nilai kernel pada setiap varaiabel
}
#penimbang kernel gabungan
q.phi=sum.q.phi/sk #nilaI kernel rata-rata
Vker=q.phi

yhat=Vker%*%y
SSR=sum((yhat-mean(y))^2)
SSE=sum((y-yhat)^2)
R2=(SSR/(SSR+SSE))*100
R2

#regresi linear
# Buat data frame dari xk dan tambahkan y1 dan y2
df = as.data.frame(xk)
colnames(df) = c("X1","X2","X3","X4")  # Ganti sesuai nama variabel bebas kamu
df$y1 = y1
df$y2 = y2
mlm_model = lm(cbind(y1, y2) ~ X1 + X2 + X3 + X4, data = df)
summary(mlm_model)


# Pisahkan kembali yhat menjadi yhat1 dan yhat2
yhat1 = yhat[1:24]
yhat2 = yhat[25:48]

# Respon asli
y1 = y[1:24]
y2 = y[25:48]

png("respon1.png", width=800, height=600)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(y1, type="l", col="blue", lwd=2, ylab="Nilai", main="Respon 1: y1 vs yhat1")
lines(yhat1, col="red", lwd=2, lty=2)
legend("topleft", legend=c("y1 (aktual)", "yhat1 (prediksi)"), col=c("blue", "red"), lty=c(1,2), lwd=2)
dev.off()



# Plot untuk respon 1
plot(y1, type="l", col="blue", lwd=2, ylab="Nilai", main="Respon 1: y1 vs yhat1")
lines(yhat1, col="red", lwd=2, lty=2)
legend("topleft", legend=c("y1 (aktual)", "yhat1 (prediksi)"), col=c("blue", "red"), lty=c(1,2), lwd=2)

# Plot untuk respon 2
plot(y2, type="l", col="blue", lwd=2, ylab="Nilai", main="Respon 2: y2 vs yhat2")
lines(yhat2, col="red", lwd=2, lty=2)
legend("topleft", legend=c("y2 (aktual)", "yhat2 (prediksi)"), col=c("blue", "red"), lty=c(1,2), lwd=2)


