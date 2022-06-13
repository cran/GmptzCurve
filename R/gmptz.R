#' Gompetrz Curve Fitting
#' @export
#' @param y a numeric vector
#' @author ARNAB ROY, DEBARGHYA BAUL
#' @importFrom graphics legend matplot
#' @importFrom stats aggregate
#' @description This function  fits the Gompertz Curve in Time Series Data along with estimates of the parameters and predicted value.
#' @details The Gompertz curve is a type of mathematical model for a Time Series. It is a sigmoid function which describes growth as being slowest at the start and end of a given time period. The equation of Gompertz curve is given by , Y=k*a^(b^x)
#' @details Taking Logarithm on both sides we get, LogY=Logk+(b^x)*Loga ; OR, Z=A+B(C^x) , where Z=LogY , A=logk , B=loga ,C=b; It is a form of Modified Exponential.So ,now we can apply Method Of Group Average.
#' @return a.hat , b.hat , k.hat   :    the estimated values of the parameters a,b and k.
#' @return  predicted.value  :   the predicted values of y
#' @examples  p=c(12,15,16,18,16,21,25,27,29,30,35,36)
#' @examples  gmptz(p)


gmptz=function(y){
  if(length(y)%%3==0){
    y=y
  }

  else if(length(y)%%3==1){
    y=y[-1]
  }
  else {
    y=y[c(-1,-2)]}
  Z=log(y)


  q=length(y)/3
  u=rep(1:3,each=q,times=1)
  S=aggregate(Z,by=list(u),sum)$x

  D1=S[2]-S[1]

  D2=S[3]-S[2]

  D=D2/D1

  C=D^(1/q)

  B=(D1*(C-1))/(C*((D-1)^2))

  A=(S[1]-(D1/(D-1)))/q

  b=C
  a=exp(B)
  k=exp(A)


  t=1:length(y)
  Y2.cal=k*(a^(b^t))




  matplot(t,cbind(y,Y2.cal),xlab = "YEAR" , ylab="y",main="GOMPERTZ CURVE",
          type="l",lwd=c(2,3),lty=c(1,2),col=c(1,2))
  legend("topleft",legend=c("Time Series Plot", "Gompetrz Curve"),
         lty=c(1,2),
         col=c(1,2))
  list(a.hat=a,b.hat=b,k.hat=k,predicted.value=Y2.cal)
}
