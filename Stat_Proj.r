cat("\f")
rm(list=ls())
if(dev.cur()>1) dev.off()
library(ggplot2)
library(colorRamps)
library(colorspace)
library(plotly)
library(RColorBrewer)
library(gridExtra)

#Beale function is in domain [-4.5,4.5] for both x and y axes. Global min is at (3, 0.5)
beale <- function(x){
  x1 = x[1]
  x2 = x[2]
  y <- (1.5 - x1 + x1*x2)^2 +(2.25 - x1 + x1*x2^2)^2 +(2.625 - x1 + x1*x2^3)^2
  return(y)
}

#Visualizing Beale function surface in 3D
x1 = seq(-4.5,4.5,length.out=100)
x2 = seq(-4.5,4.5,length.out=100)
z = as.data.frame(matrix(0,nrow=length(x1)*length(x2),ncol=3))
count=1
for(i in 1:length(x1)){
  for(j in 1:length(x2)){
    z[count,] = c(x1[i],x2[j],beale(c(x1[i],x2[j])))
    count=count+1
  }
}
colnames(z)=c("x","y","z")

# plot_ly(z=~z) %>% add_contour(x=x1,y=x2)
#Try to add 3D grid lines later with scatter3d + add_surface=T
# plot_ly(type = 'scatter3d' , x = , y =, z = ,mode = 'lines')

#Stochastic Gradient Descent algorithm
stochastic_gradient_descent <- function(theta_init,step=.0001,eps=1e-4,tol=1e-3){
  thetas = matrix(NA,nrow=1000,ncol=2)
  fs = rep(NA,1000)
  dfs = matrix(NA,1000,2)
  fs = rep(NA,1000)
  x = as.vector(theta_init)
  f.x = beale(x)
  fs[1] = f.x
  f.x.plus.eps = c(beale(x+c(eps,0)), beale(x+c(0,eps)))
  df.x =  (f.x.plus.eps - f.x)/eps
  thetas[1,] = theta_init
  dfs[1,] = df.x
  count=1
  # any(abs(dfs[count]) > tol)
  while(1){
    x = thetas[count,]
    f.x = fs[count]
    df.x = dfs[count,] 
    count= count+1
    # print(count)
    x.new = x - step*df.x
    thetas[count,] = x.new
    f.x = beale(x.new)
    fs[count] = f.x
    f.x.plus.eps = c(beale(x.new+c(eps,0)), beale(x.new+c(0,eps)))
    df.x =  (f.x.plus.eps - f.x)/eps
    if(any(!is.finite(df.x))) break
    dfs[count,] = df.x
    if(count%%1000==0){
      fs = c(fs,rep(NA,1000))
      dfs = rbind(dfs,matrix(NA,1000,2))
      thetas = rbind(thetas,matrix(NA,nrow=1000,ncol=2))
    }
    if(abs(fs[count]-fs[count-1])<tol) break
    if(count==10000) break
  }
  thetas = na.omit(thetas)
  attributes(thetas)$na.action = NULL
  fs = as.vector(na.omit(fs))
  dfs = na.omit(dfs)
  attributes(dfs)$na.action = NULL
  return (data.frame("thetas"=thetas,"fs"=fs,"dfs"= dfs))
}

x.in = c(-3,3)
t1 = Sys.time()
Final.res = stochastic_gradient_descent(x.in)
t2 = Sys.time()
print(t2-t1)
sgd.res = Final.res
print("SGD res:")
print(tail(sgd.res,1))
# df=sgd.res
# colnames(df) = c("Beta1","Beta2","Fx","dfx1","dfx2")
# plots=list()
# plots[[1]] = ggplot(df,aes(1:nrow(df),Fx))+geom_point()+geom_line()+xlab("i")+ylab("f(x)")+theme_classic()
# plots[[2]] = ggplot(df,aes(1:nrow(df),dfx1))+geom_point()+geom_line()+xlab("i")+ylab("f'(x1)")+theme_classic()
# plots[[3]] = ggplot(df,aes(1:nrow(df),dfx2))+geom_point()+geom_line()+xlab("i")+ylab("f'(x2)")+theme_classic()
# grid.arrange(grobs=plots,nrow=1,ncol=3)

#Stochastic Gradient Descent algorithm with momentum
momentum <- function(theta_init,step=.000001,eps=1e-4,tol=1e-3,gamma=0.9){
  thetas = matrix(NA,nrow=1000,ncol=2)
  fs = rep(NA,1000)
  dfs = matrix(NA,1000,2)
  fs = rep(NA,1000)
  x = as.vector(theta_init)
  f.x = beale(x)
  fs[1] = f.x
  f.x.plus.eps = c(beale(x+c(eps,0)), beale(x+c(0,eps)))
  df.x =  (f.x.plus.eps - f.x)/eps
  thetas[1,] = theta_init
  dfs[1,] = df.x
  count=1
  prev.update = 0
  while(1){
    x = thetas[count,]
    f.x = fs[count]
    df.x = dfs[count,] 
    count= count+1
    # print(count)
    x.new = x - gamma*prev.update- step*df.x
    prev.update = x-x.new
    thetas[count,] = x.new
    f.x = beale(x.new)
    fs[count] = f.x
    f.x.plus.eps = c(beale(x.new+c(eps,0)), beale(x.new+c(0,eps)))
    df.x =  (f.x.plus.eps - f.x)/eps
    if(any(!is.finite(df.x))) break
    dfs[count,] = df.x
    if(count%%1000==0){
      fs = c(fs,rep(NA,1000))
      dfs = rbind(dfs,matrix(NA,1000,2))
      thetas = rbind(thetas,matrix(NA,nrow=1000,ncol=2))
    }
    if(abs(fs[count]-fs[count-1])<tol) break
    if(count==10000) break
  }
  thetas = na.omit(thetas)
  attributes(thetas)$na.action = NULL
  fs = as.vector(na.omit(fs))
  dfs = na.omit(dfs)
  attributes(dfs)$na.action = NULL
  return (data.frame("thetas"=thetas,"fs"=fs,"dfs"= dfs))
}


t1 = Sys.time()
Final.res = momentum(x.in)
mom.res = Final.res
print("Momentum res:")
print(tail(mom.res,1))
t2 = Sys.time()
print(t2-t1)
# df = mom.res
# colnames(df) = c("Beta1","Beta2","Fx","dfx1","dfx2")
# plots=list()
# plots[[1]] = ggplot(df,aes(1:nrow(df),Fx))+geom_point()+geom_line()+xlab("i")+ylab("f(x)")+theme_classic()
# plots[[2]] = ggplot(df,aes(1:nrow(df),dfx1))+geom_point()+geom_line()+xlab("i")+ylab("f'(x1)")+theme_classic()
# plots[[3]] = ggplot(df,aes(1:nrow(df),dfx2))+geom_point()+geom_line()+xlab("i")+ylab("f'(x2)")+theme_classic()
# grid.arrange(grobs=plots,nrow=1,ncol=3)

#Nested Accelerated Gradient
nag <- function(theta_init,step=.000002,eps=1e-4,tol=1e-3,gamma=0.9){
  thetas = matrix(NA,nrow=1000,ncol=2)
  fs = rep(NA,1000)
  dfs = matrix(NA,1000,2)
  fs = rep(NA,1000)
  x = as.vector(theta_init)
  f.x = beale(x)
  fs[1] = f.x
  f.x.plus.eps = c(beale(x+c(eps,0)), beale(x+c(0,eps)))
  df.x =  (f.x.plus.eps - f.x)/eps
  thetas[1,] = theta_init
  dfs[1,] = df.x
  count=1
  prev.update = 0
  while(1){
    x = thetas[count,]
    f.x = fs[count]
    df.x = dfs[count,] 
    count= count+1
    # print(count)
    x.new = x - gamma*prev.update- step*df.x
    prev.update = x-x.new
    thetas[count,] = x.new
    #Changing current theta to approx future theta for next update (NAG step)
    x.new = x.new - gamma*prev.update
    f.x = beale(x.new)
    fs[count] = f.x
    f.x.plus.eps = c(beale(x.new+c(eps,0)), beale(x.new+c(0,eps)))
    df.x =  (f.x.plus.eps - f.x)/eps
    if(any(!is.finite(df.x))) break
    dfs[count,] = df.x
    if(count%%1000==0){
      fs = c(fs,rep(NA,1000))
      dfs = rbind(dfs,matrix(NA,1000,2))
      thetas = rbind(thetas,matrix(NA,nrow=1000,ncol=2))
    }
    if(abs(fs[count]-fs[count-1])<tol) break
    if(count==10000) break
  }
  thetas = na.omit(thetas)
  attributes(thetas)$na.action = NULL
  fs = as.vector(na.omit(fs))
  dfs = na.omit(dfs)
  attributes(dfs)$na.action = NULL
  return (data.frame("thetas"=thetas,"fs"=fs,"dfs"= dfs))
}


t1 = Sys.time()
Final.res = nag(x.in)
t2 = Sys.time()
print(t2-t1)
nag.res = Final.res
print("Nag res:")
print(tail(nag.res,1))
# colnames(df) = c("Beta1","Beta2","Fx","dfx1","dfx2")
# plots=list()
# plots[[1]] = ggplot(df,aes(1:nrow(df),Fx))+geom_point()+geom_line()+xlab("i")+ylab("f(x)")+theme_classic()
# plots[[2]] = ggplot(df,aes(1:nrow(df),dfx1))+geom_point()+geom_line()+xlab("i")+ylab("f'(x1)")+theme_classic()
# plots[[3]] = ggplot(df,aes(1:nrow(df),dfx2))+geom_point()+geom_line()+xlab("i")+ylab("f'(x2)")+theme_classic()
# grid.arrange(grobs=plots,nrow=1,ncol=3)

#adadelta
adadelta <- function(theta_init,eps=1e-4,tol=1e-3,err=1e-6,gamma=0.9){
  thetas = matrix(NA,nrow=1000,ncol=2)
  G =c(0,0)
  E = c(0,0)
  delta.x=c(0,0)
  fs = rep(NA,1000)
  dfs = matrix(NA,1000,2)
  fs = rep(NA,1000)
  x = as.vector(theta_init)
  f.x = beale(x)
  fs[1] = f.x
  f.x.plus.eps = c(beale(x+c(eps,0)), beale(x+c(0,eps)))
  df.x =  (f.x.plus.eps - f.x)/eps
  thetas[1,] = theta_init
  dfs[1,] = df.x
  count=1
  # any(abs(dfs[count]) > tol)
  while(1){
    x = thetas[count,]
    f.x = fs[count]
    df.x = dfs[count,] 
    G = gamma*G + (1-gamma)*(df.x)**2
    delta.x = - sqrt((E+err)/(G+err)) * df.x
    count= count+1
    x.new = x +delta.x
    # #print(count)
    E = gamma*E + (1-gamma)*(delta.x)**2
    thetas[count,] = x.new
    f.x = beale(x.new)
    fs[count] = f.x
    f.x.plus.eps = c(beale(x.new+c(eps,0)), beale(x.new+c(0,eps)))
    df.x =  (f.x.plus.eps - f.x)/eps
    if(any(!is.finite(df.x))) break
    dfs[count,] = df.x
    if(count%%1000==0){
      fs = c(fs,rep(NA,1000))
      dfs = rbind(dfs,matrix(NA,1000,2))
      thetas = rbind(thetas,matrix(NA,nrow=1000,ncol=2))
    }
    if(abs(fs[count]-fs[count-1])<tol) break
    # if(fs[count]>fs[count-1]) print("failed descent")
    if(count==5*10**4) break
  }
  thetas = na.omit(thetas)
  attributes(thetas)$na.action = NULL
  fs = as.vector(na.omit(fs))
  dfs = na.omit(dfs)
  attributes(dfs)$na.action = NULL
  return (data.frame("thetas"=thetas,"fs"=fs,"dfs"= dfs))
}


t1 = Sys.time()
Final.res = adadelta(x.in)
t2 = Sys.time()
print(t2-t1)
adadelta.res = Final.res
print("adadelta res:")
print(tail(adadelta.res,1))
# colnames(df) = c("Beta1","Beta2","Fx","dfx1","dfx2")
# plots=list()
# plots[[1]] = ggplot(df,aes(1:nrow(df),Fx))+geom_point()+geom_line()+xlab("i")+ylab("f(x)")+theme_classic()
# plots[[2]] = ggplot(df,aes(1:nrow(df),dfx1))+geom_point()+geom_line()+xlab("i")+ylab("f'(x1)")+theme_classic()
# plots[[3]] = ggplot(df,aes(1:nrow(df),dfx2))+geom_point()+geom_line()+xlab("i")+ylab("f'(x2)")+theme_classic()
# grid.arrange(grobs=plots,nrow=1,ncol=3)

adam <- function(theta_init,eps=1e-4,tol=1e-3,err=1e-8,B1=0.9,B2=0.999,step=.001){
  thetas = matrix(NA,nrow=1000,ncol=2)
  m = c(0,0)
  v = c(0,0)
  fs = rep(NA,1000)
  dfs = matrix(NA,1000,2)
  fs = rep(NA,1000)
  x = as.vector(theta_init)
  f.x = beale(x)
  fs[1] = f.x
  f.x.plus.eps = c(beale(x+c(eps,0)), beale(x+c(0,eps)))
  df.x =  (f.x.plus.eps - f.x)/eps
  thetas[1,] = theta_init
  dfs[1,] = df.x
  count=1
  # any(abs(dfs[count]) > tol)
  while(1){
    x = thetas[count,]
    f.x = fs[count]
    df.x = dfs[count,] 
    m = B1*m + (1-B1)*df.x
    v = B2*v + (1-B2)*df.x**2
    mb = m/(1-B1**count)
    vb = v/(1-B2**count)
    count= count+1
    x.new = x - step * ((vb**0.5 + err)**-1 )*mb
    # print(count)
    thetas[count,] = x.new
    f.x = beale(x.new)
    fs[count] = f.x
    f.x.plus.eps = c(beale(x.new+c(eps,0)), beale(x.new+c(0,eps)))
    df.x =  (f.x.plus.eps - f.x)/eps
    if(any(!is.finite(df.x))) break
    dfs[count,] = df.x
    if(count%%1000==0){
      fs = c(fs,rep(NA,1000))
      dfs = rbind(dfs,matrix(NA,1000,2))
      thetas = rbind(thetas,matrix(NA,nrow=1000,ncol=2))
    }
    if(fs[count]>fs[count-1]) print("failed descent")
    if(abs(fs[count]-fs[count-1])<tol) break
    
    if(count==10000) break
  }
  thetas = na.omit(thetas)
  attributes(thetas)$na.action = NULL
  fs = as.vector(na.omit(fs))
  dfs = na.omit(dfs)
  attributes(dfs)$na.action = NULL
  return (data.frame("thetas"=thetas,"fs"=fs,"dfs"= dfs))
}


t1 = Sys.time()
Final.res = adam(x.in)
t2 = Sys.time()
print(t2-t1)
adam.res = Final.res
print("adam res:")
print(tail(adam.res,1))
# colnames(df) = c("Beta1","Beta2","Fx","dfx1","dfx2")
# plots=list()
# plots[[1]] = ggplot(df,aes(1:nrow(df),Fx))+geom_point()+geom_line()+xlab("i")+ylab("f(x)")+theme_classic()
# plots[[2]] = ggplot(df,aes(1:nrow(df),dfx1))+geom_point()+geom_line()+xlab("i")+ylab("f'(x1)")+theme_classic()
# plots[[3]] = ggplot(df,aes(1:nrow(df),dfx2))+geom_point()+geom_line()+xlab("i")+ylab("f'(x2)")+theme_classic()
# grid.arrange(grobs=plots,nrow=1,ncol=3)


  # stat_contour(z,aes(x=x,y=y,z=log10(z),color=..level..)) +guides(color = guide_colorbar("Value")) +
ggplot() + theme_classic()+ xlab("x1") +  ylab('x2') +
  stat_contour(data=z,aes(x=x,y=y,z=log10(z),color=..level..)) +guides(color = guide_colorbar("Value"))+
  geom_line(data=sgd.res,aes(x=thetas.1,y=thetas.2),color="black") +
  geom_line(data=mom.res,aes(x=thetas.1,y=thetas.2),color="blue") +
  geom_line(data=nag.res,aes(x=thetas.1,y=thetas.2),color="red") +
  geom_line(data=adadelta.res,aes(x=thetas.1,y=thetas.2),color="green") +
  geom_line(data=adam.res,aes(x=thetas.1,y=thetas.2),color="yellow")+
labs(subtitle="Black:SGD Blue:Mom Red:NAG Green:adadelta Yellow:Adam")
   
preds <- function(x)  {return(as.numeric(x))}

predictions = rbind("SGD"=
                  preds(tail(sgd.res[,1:3],1)),
                "Mom"=
                  preds(tail(mom.res[,1:3],1)),
                "nag"=
                  preds(tail(nag.res[,1:3],1)),
                "adadelta"=
                  preds(tail(adadelta.res[,1:3],1)),
                "adam"=
                  preds(tail(adam.res[,1:3],1))
)
print((predictions))
