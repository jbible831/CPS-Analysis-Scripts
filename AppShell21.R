

#setwd("C:\\Users\\biblejd\\Desktop\\CPS reboot\\Analysis")



load('UseThisOne')



library("matrixStats")
library(geepack)
library(glmmML)
library(lme4)
library(glmm)




head(lag[,1:3], n=100)


#lose 99 to non reported alcohol or druguse
#good bit of overlap 
#77 and 82 respectively from each

gmat<-as.matrix(cbind(
lag$maternalid,
lag$pregnum,
lag$ptb,
lag$pre_ptb,
lag$ipiM,
lag$ipiMa,
lag$ipiMb,
lag$momage,
lag$dsmoke,
lag$dalcohol,
lag$ddruguse,
lag$dpredm,
lag$dchronhtn,
lag$dhxheartdis,
lag$dhxrenaldis,
lag$dhxdepression,
lag$dhxseizure,
lag$dhxthyroid,
lag$dhxasthma,
lag$dantiuti,
lag$dantistd,
lag$dantibleed,
lag$dprevia,
lag$daccreta, 
lag$lagage
))
colnames(gmat)<-
c('maternalid','pregnum','ptb','pre_ptb','ipiM','ipiMa','ipiMb',
	'momage','dsmoke','dalcohol','ddruguse',
	'dpredm','dchronhtn','dhxheartdis','dhxrenaldis',
	'dhxdepression','dhxseizure','dhxthyroid','dhxasthma',
	'dantiuti','dantistd','dantibleed','dprevia','daccreta', 
	'lagage')


uid<-unique(lag$maternalid)
indices<-lapply(uid, FUN= function(i) l<- which(lag$maternalid==i))
drops<-c()
for (i in 1:length(uid)){
if(head(lag[unlist(indices[i]),2])[1] != 2) drops<- c(drops,unlist(indices[i]))
}





#
##
###
#####
#######
#########
############
##############
#################
##################
###################
####################
nullis<-lag[-drops,]
lag<-nullis
#######
#####
###
##
#


#cc<-glmmML(ptb~	pregnum+pre_ptb + momage+ipiMa+ipiMb+
#	prebmi_new + momraceother +dnotpinsurance , data=nullis, 
#	, family=binomial, cluster= maternalid, method='ghq', n.points=15)



#parameters for unadusted model (excluding) gaptime 
#taken from glmmML (most stable) 
# pre_ptb_1, momage_1, ipia,  ipib


uid<-unique(lag$maternalid)



lag$bday<-as.Date(as.character(lag$birthday), "%m/%d/%Y")
mbday<-max(as.Date(as.character(lag$birthday), "%m/%d/%Y"))
lag$t2end<- as.numeric(((mbday+.01)-lag$bday))/(365/12)

mat<-as.matrix(cbind(
	lag$maternalid,
	lag$pregnum,
	lag$ptb,
	lag$pre_ptb,
	lag$momage,
	lag$ipiMa,
	lag$ipiMb,
	lag$prebmi_new,
	lag$momraceother,
	lag$dnotpinsurance,
	lag$t2end, 
	lag$lagage
	))
colnames(mat)<-c("maternalid","pregnum","ptb",
	"pre_ptb","momage","ipiMa","ipiMb","prebmi_new","momraceother",
	"dnotpinsurance", "t2end", "lagage"
	)



rest<-lapply(uid, FUN= function(i) l<- which(lag$maternalid==i))
last<-unlist(lapply(uid, FUN= function(i) l<- tail(which(lag$maternalid==i), n=1)))
gapmat<-matrix(0, ncol= ncol(mat)+1, nrow= length(c(last, unlist(rest))))
colnames(gapmat)<-c("maternalid","pregnum","ptb",
	"pre_ptb","momage","ipiMa","ipiMb","prebmi_new","momraceother",
	"dnotpinsurance", "t2end", "lagage","cens"
	)


k<-1
i<-1
for (i in 1:length(uid)){
	newage<-round(sum(mat[unlist(last[i]),11]/12))
	lvec<-rbind(c(mat[unlist(last[i]),],1))+c(0,1,0,0,0,rep(0,6),mat[unlist(last[i]),5]-mat[unlist(last[i]),12] , 0)
 	fvec<-matrix(c((mat[unlist(rest[i]),]),rep(0, length(unlist(rest[i])))), nrow=length(unlist(rest[i])))
	ttei<-rbind(fvec, lvec)
	gapmat[k:(k+nrow(ttei)-1),]<-as.matrix(ttei)
	k<-k+nrow(ttei)
}
head(gapmat,n=100)
##
#
#
#
#
#
#
#

#
##
#

#
#
##
##
#
#
#
#
#
#

#
hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}


gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}

mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}






#Number of quad points to use

idinds<-lapply(unique(gapmat[,1]), FUN= l<-function(i) which(gapmat[,1]==i))


#Likelihood



#######################################################
###the original version of this can be found 		#
###Here\/ 								#
###C:\Users\biblejd\Desktop\NEW CPS\older\newcont	#
#######################################################

points<-2
lff<-function(b){
#Extracts parameters (just saves you the trouble of having to index everything)
par<-b
	beta<-par[2:9]
	alpha1<- par[10]
	tau<-par[11:14]
	gam<-par[15:23]
	sigma1<-exp(par[1])**2
	sigma2<-exp(par[24])**2
	alpha2<-par[26]	

#	xs<-x*(exp((par[1])))
	#variance components (std-dev pars log trans, rho-par logit trans)
	dud<-2*(plogis(par[25]/10)-.5)
#	dud<-((-rhotol<dud)&(dud<rhotol))*dud+(1-(dud<rhotol)) * rhotol + (1-(dud>-rhotol))* -rhotol
	rho<-dud
		sigma<- matrix(c(sigma1,rho*sqrt(sigma1 *sigma2),rho*sqrt(sigma1*sigma2),  sigma2) , ncol=2, nrow=2)
		mghp<-mgauss.hermite(points, c(0,0), sigma, prune=NULL)
		b1<-mghp$points[,1]
		b2<-mghp$points[,2]
		w<-mghp$weights	
#
	gapcon<-function(a1, a2){
#cont. portion
		p<-(plogis((tau[1]+gapmat[,c(2,12, 4)]%*%tau[-1])    +(a1)*alpha1+(a2) *alpha2))
		p[which(gapmat[,2]==2)]<-1
		p[which(gapmat[,2]==6)]<-0
#Trans portion 
		(
		(
		(plogis((gam[1]+gapmat[,-c(1,3,11,12:13)]%*%gam[-1])+a2)*gapmat[,3])+
		((1-plogis((gam[1]+gapmat[,-c(1,3,11,12:13)]%*%gam[-1])+a2)))*
		(1-gapmat[,3])
		)
		)*
#gap and continuation portion 
		(
		(p*dgamma(rowSums(gapmat[,6:7]),shape=exp(beta[1]+gapmat[,c(2,4,12,8,9,10)]%*%beta[-c(1,8)] +a1), scale= exp(beta[8]))*(1-gapmat[,13])
		)
		)+ 
		((p*(1-pgamma(gapmat[,11],shape=exp(beta[1]+gapmat[,c(2,4,12,8,9,10)]%*%beta[-c(1,8)]+a1), scale= exp(beta[8])))+(1-p))*gapmat[,13])
		
	}


qpnts<-cbind(b1, b2)


stuff<-matrix(unlist(lapply(1:nrow(qpnts), FUN= raw<- function(i) (gapcon(qpnts[i,1], qpnts[i,2]))   )), ncol=nrow(gapmat), byrow=T)
stuff2<-matrix(0, ncol=length(idinds), nrow=points**2)
for (i in 1:length(idinds)) { 
	stuff2[,i]<- rowProds(stuff[, unlist(idinds[i])])
}

	#likelihood sort of for the gh contributions
			ll<-sum(-log(rowSums(t(stuff2)%*% diag(w))))

#neg2loglikebi<- -log(rowSums(t(stuff)%*% diag(w)))


	return(ll)
}




par1<-c(c(
#[1] log(sigma1)
-5,
#[2:9] gaptime model pars
0.943279886,  0.120101050, -0.058825879,  0.002050170,-0.001321664, -0.093218646, -0.149117802,  1.817231205,
#[10] alpha1 the shared parameter with gap and cont. models
0,
#[11:14] cont. model parameters 
0.932616052, -0.048161928, -0.047874830, -0.237927623 ),
#[15:23] trans model parameters
c(-1.876926864, 0.138553283,  1.681749831, -0.022364707, -0.070807026,  0.008916330,-0.000306610, -0.043459195,  0.065590602), 
#[24] log(sigma2)
-0.536636180 ,
#[25] transformed(rho) i.e. the correlation between gap and trans re
0, 
#[26] alpha2 
-0.218195032 
)


par2<-c(c(
#[1] log(sigma1)
-1.7801,
#[2:9] gaptime model pars
0.900137,0.198217,-0.04793,
0.001007,-0.00088,-0.08606,-0.15517,1.70224,
#[10] alpha1 the shared parameter with gap and cont. models
-4.28543,
#[11:14] cont. model parameters 
2.228983,-0.56691,-0.04042,-0.29551),
#[15:23] trans model parameters
c(-1.73793,0.122712,1.684995,-0.02226,
-0.06785,0.008652,-0.0006,-0.04382,
0.065834), 
#[24] log(sigma2)
-5,
#[25] transformed(rho) i.e. the correlation between gap and trans re
0, 
#[26] alpha2 relation between trans and cont
0
)

par3<-c(c(
#[1] log(sigma1)
-1.7801,
#[2:9] gaptime model pars
0.900137,0.198217,-0.04793,
0.001007,-0.00088,-0.08606,-0.15517,1.70224,
#[10] alpha1 the shared parameter with gap and cont. models
-4.28543,
#[11:14] cont. model parameters 
2.228983,-0.56691,-0.04042,-0.29551),
#[15:23] trans model parameters
c(-1.73793,0.122712,1.684995,-0.02226,
-0.06785,0.008652,-0.0006,-0.04382,
0.065834), 
#[24] log(sigma2)
-0.536636180,
#[25] transformed(rho) i.e. the correlation between gap and trans re
0, 
#[26] alpha2 relation between trans and cont
-0.218195032
)



par4<-c(c(
#[1] log(sigma1)
0,
#[2:9] gaptime model pars
0.900137,0.198217,-0.04793,
0.001007,-0.00088,-0.08606,-0.15517,1.70224,
#[10] alpha1 the shared parameter with gap and cont. models
0,
#[11:14] cont. model parameters 
2.228983,-0.56691,-0.04042,-0.29551),
#[15:23] trans model parameters
c(-1.73793,0.122712,1.684995,-0.02226,
-0.06785,0.008652,-0.0006,-0.04382,
0.065834), 
#[24] log(sigma2)
0,
#[25] transformed(rho) i.e. the correlation between gap and trans re
0, 
#[26] alpha2 relation between trans and cont
0
)


points<- 21

lff(par4)
lff(par3)
lff(par2)
lff(par1)


#DUMZ<-1

dudz<-paste('par',DUMZ, sep='')
lff(eval(parse(text=dudz)))




start<-proc.time()
dumb2<-optim(eval(parse(text=dudz)), lff, method="BFGS", control=c(maxit=500, trace=T)  , hessian=T)
start-proc.time()
save(dumb2, file=paste("res","APP", sep="") )



