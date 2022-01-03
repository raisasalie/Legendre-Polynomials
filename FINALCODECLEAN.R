rm(list = ls(all = TRUE))
install.packages("wesanderson")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)

# QUESTION 1 i)
# Legendre polynomials
Legendre=function(x,n){
  val=0
  for(i in 0:n){
    val=val+((x^i)*choose(n,i)*choose((n+i-1)/2,n))
  }
  return((2^n)*val)
}

# plot Legendre polynomials of different orders
x<-seq(-1,1,0.01)
tans <- c("blue", "tan3", "darkgreen", "orange","darkblue")
plot(x,Legendre(x,1),type="l",
     xlab=expression(x),
     ylab=expression(L[n](x)),
     col=tans[1],lwd=2,
     main="Legendre polynomials")
lines(x,Legendre(x,0),type="l",col="black",lwd=2)
lines(x,Legendre(x,2),type="l",col=tans[2],lwd=2)
lines(x,Legendre(x,3),type="l",col=tans[3],lwd=2)
lines(x,Legendre(x,4),type="l",col=tans[4],lwd=2)
lines(x,Legendre(x,5),type="l",col=tans[5],lwd=2)
legend("bottomright", title = expression(q),
       legend = c(0:5), 
       lty = 1, lwd = 2, 
       col=c("black", tans), 
       horiz = T)

# QUESTION 1 ii)
# creates random polynomial target functions
rLegfunc=function(x,n){
  val=0
  vec=runif(n,-1,1)
  for(i in 1:n){
    val=val+(Legendre(x,i)*vec[i])
  }
  return(val)
}

rPolyfunc=function(x,n){
  val=0
  vec=runif(n,-1,1)
  for(i in 1:n){
    val=val+(x^i*vec[i])
  }
  return(val)
}

# plot random target functions 
x <- seq(-1,1,by = 0.01)
set.seed(2020)
par(mfrow=c(1,2))
# polynomial targets
tans <- c("blue", "darkgreen", "tan3")
plot(x,rPolyfunc(x,4),type="l",
     ylab=expression(f(x)),
     xlab=expression(x),
     col=tans[1], 
     ylim=c(-2,2), 
     lwd=2,
     main="Random polynomial target functions")
lines(x,rPolyfunc(x,5),type="l",col=tans[2], lwd=2)
lines(x,rPolyfunc(x,6),type="l",col=tans[3], lwd=2)
legend("bottomleft", title = expression(Q[f]),
       legend = c(5:7), 
       lty = 1, lwd = 2, 
       col=tans, horiz = T)
# Legendre polynomial targets
plot(x,rLegfunc(x,4),type="l",
     ylab=expression(f(x)),
     xlab=expression(x),
     col=tans[1], lwd=2,
     ylim=c(-2,2),
     main="Random Legendre target functions")
lines(x,rLegfunc(x,5),type="l",col=tans[2], lwd=2)
lines(x,rLegfunc(x,6),type="l",col=tans[3], lwd=2)
par(mfrow=c(1,1))

#---------------------------------------------
# QUESTION 2i)
# generate random 10th order target function over [-1,1]
set.seed(2020)
func <- rLegfunc(x, 10)
plot(x, func, type='l')

#simulates a dataset with target "func" of size n with noise sig
generator=function(n,x,func,sig){
  l<-length(func)
  dat<-matrix(rep(NA,2*n),ncol=n)
  xdat<-floor(runif(n)*l)+1
  ydat<-func[xdat]+rnorm(n,0,sig)
  xdat<-x[xdat]
  Data<-data.frame(xdat,ydat)
  return(Data)
}

#fitted model of the data
fit=function(x,model){
  v=0
  for(i in 1:length(model$coefficient)){
    v=v+(model$coefficient[i]*(x^(i-1)))
  }
  return(v)
}

# gives the bias for a given fitted model (DOES NOT ACCOUNT FOR VARIABILITY IN THE DATA)
fdiff=function(x,target,model){
  f=fit(x,model)
  return((t(f-target)%*%(f-target))*(x[2]-x[1]))
}

# -------------------

Nseq <- seq(20,110,1)
sigseq <- seq(0.2,1.1,0.01)

# results
res1long <- data.frame(N=0, sigma=0, i=0,overfit=0)
colnames(res) <- paste(Nseq)

set.seed(2020)
for (i in 1:100){
  for (n in 1:length(Nseq)){
    N <- Nseq[n]
    for (s in 1:length(sigseq)){
      sig <- sigseq[s] 
      
      #generate data
      data <- generator(n = N, x = x, func = func, sig = sig)
      # fit second/ten order polynomial to data 
      two <- lm(I(data$ydat) ~ data$xdat + I(data$xdat^2))
      ten <- lm(I(data$ydat) ~ data$xdat + I(data$xdat^2) + 
                  I(data$xdat^3) + I(data$xdat^4)+ 
                  I(data$xdat^5)+ I(data$xdat^6)+
                  I(data$xdat^7)+ I(data$xdat^8)+ 
                  I(data$xdat^9)+ I(data$xdat^10))
      err.out2<-fdiff(x,func,two)+(sig^2) # Out of sample error
      err.out10<-fdiff(x,func,ten)+(sig^2) # Out of sample error
      # difference between E_out for each model
      overfit <- err.out10 - err.out2
      # store 
      res1long <- rbind(res1long, c(N=N, sigma = sig, i=i ,overfit = overfit)) 
    }
  }
}

beep()
res1long <- res1long[-1,]

# find means across different i
res1df <- data.frame(sigma =0, N=0, E_overfit = 0)
for (s in sigseq){
  sig.res <- res1long[which(res1long$sigma==s),]
  for (this.N in Nseq){
    res1df <- rbind(res1df, 
                    cbind(sigma = s, 
                          N =this.N, 
                          E_overfit = mean(sig.res[which(sig.res$N==this.N),]$overfit)))
  }
}
res1df <- res1df[-1,]

# scale for better colour mapping
res1df$E_overfit <- log(res1df$E_overfit - min(res1df$E_overfit))

# make a palette 
pal <- wes_palette("Darjeeling2", 100000, type = "continuous")

# colour map
plot1 <- ggplot(res1df, aes(N, sigma)) +
  geom_raster(aes(fill = E_overfit)) +
  theme_minimal() + 
  ylab(expression(sigma))+
  scale_fill_gradientn(colors=pal, name = "colour map")+
  xlim(20,110) +
  ylim(0.2,1.1) +
  ggtitle(expression(log(E[out](g[10]) - E[out](g[2]) - 
                           min(group("(",E[out](g[10]) - E[out](g[2]),")")) + epsilon)))

plot1

#--------------------
# QUESTION 2 ii.)
Qfseq <- seq(1,40,by=1)
Nseq <- seq(20,60,by=1)
sig <- 0.2
Nreps <- 20 # number of repeats per value of Qf
Ndatasets <- 10 # number of datasets per Qf, N combination 

# storage for results
res2long <- c()

set.seed(2020)
for (q in 1:length(Qfseq)){
  Qf <- Qfseq[q]
  func <- rLegfunc(x, Qf)
  
  for (i in 1:Nreps){
    for (n in 1:length(Nseq)){
      N <- Nseq[n]
      
      for (j in 1:Ndatasets){
        # generate data from target function
        data <- generator(n = N, x = x, func = func, sig = sig)
        # fit second/tenth order polynomial to data 
        two <- lm(I(data$ydat) ~ data$xdat + I(data$xdat^2))
        ten <- lm(I(data$ydat) ~ data$xdat + I(data$xdat^2) + I(data$xdat^3) + 
                    I(data$xdat^4)+ I(data$xdat^5)+ 
                    I(data$xdat^6)+ I(data$xdat^7)+ 
                    I(data$xdat^8)+ I(data$xdat^9)+ 
                    I(data$xdat^10))
        err.out2<-fdiff(x,func,two)+(sig^2) # Out of sample error
        err.out10<-fdiff(x,func,ten)+(sig^2) # Out of sample error
        # difference between E_out
        overfit <- err.out10 - err.out2
        res2long <- rbind(res2long, c(Qf, N, i, j, overfit))
        
      }
    }
  }
}

colnames(res2long) <- c("Qf", "N", "i", "j", "overfit")
res2long <- as.data.frame(res2long)

# find means across different i,j
res2df <- data.frame(Qf =0, N=0, E_overfit = 0)
for (this.q in Qfseq){
  Qf.res <- res2long[which(res2long$Qf==this.q),]
  for (this.N in Nseq){
    res2df <- rbind(res2df, 
                    cbind(Qf = this.q, 
                          N =this.N, 
                          E_overfit = mean(Qf.res[which(Qf.res$N==this.N),]$overfit)))
  }
}
res2df <- res2df[-1,]

# scale for better colour mapping 
res2df$E_overfit <-log((res2df$E_overfit - min(res2df$E_overfit)) / (-min(res2df$E_overfit)))

# create a palette
pal2 <- wes_palette("Darjeeling2", n = 100000, type = "continuous")

# colour map
plot2<- ggplot(res2df, aes(N, Qf)) +
  geom_raster(aes(fill = E_overfit)) +
  theme_minimal() +
  ylab(expression(Q[f]))+
  scale_fill_gradientn(colors=pal2, name="colour map")+
  xlim(20,60) +
  ylim(0,40) +
  ggtitle(expression(log(frac(E[out](g[10]) - E[out](g[2]) - min(group("(",E[out](g[10]) - 
                                                                         E[out](g[2]),")")) + epsilon, 
                              - min(group("(",E[out](g[10]) - E[out](g[2]),")"))))))

plot2
