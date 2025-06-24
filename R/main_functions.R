F_stat<-function(d_mat,predictor,confonding_stat=FALSE,confonding_var=NA,r){
  n<-nrow(d_mat)
  sim_mat<-(-1/2)*d_mat^2
  ide_mat<-diag(1,n,n)
  one_vec<-matrix(1, nrow = n, ncol = 1)
  pre_mat<-as.matrix(scale(predictor, center=T, scale = FALSE))
  G_mat<-(ide_mat-(one_vec%*%t(one_vec)/n))%*%sim_mat%*%(ide_mat-(one_vec%*%t(one_vec)/n))
  G_mat_eigen<-eigen(G_mat)$values
  G_mat_eigen<-abs(G_mat_eigen)
  G_mat_vector<-matrix(eigen(G_mat)$vectors,n,n)
  G_mat_r_eigen<-G_mat_eigen^r
  G_mat_r<-(ide_mat-(one_vec%*%t(one_vec)/n))%*%G_mat_vector%*%diag(G_mat_r_eigen)%*%t(G_mat_vector)%*%(ide_mat-(one_vec%*%t(one_vec)/n))
  if(confonding_stat==T){
    confonding_mat<-as.matrix(scale(confonding_var ,center=TRUE, scale = FALSE))
    X_mat<-cbind(pre_mat,confonding_mat)
    H_mat_confonding<-confonding_mat%*%ginv(t(confonding_mat)%*%confonding_mat)%*%t(confonding_mat)
    H_mat_X<-X_mat%*%ginv(t(X_mat)%*%X_mat)%*%t(X_mat)
    H_mat<-H_mat_X-H_mat_confonding
  }
  else{
    H_mat<-pre_mat%*%ginv(t(pre_mat)%*%pre_mat)%*%t(pre_mat)
  }
  F_stat<-tr(H_mat%*%G_mat_r)
  return(list(F_stat=F_stat,H_mat=H_mat, G_mat_r=G_mat_r))
}

permutation_null<-function(A,W){
  n=nrow(A)
  Fstar=tr(A%*%W)
  mean.null=tr(A)*tr(W)/(n-1)	## mean

  T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
  Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
  temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
  temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
  temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
  temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
  temp2=temp21*temp22/temp23
  variance.null=temp1+temp2		## variance

  T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
  T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
  t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
  t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
  t3=24*(n^2-n-4)*(U*Bs+B*Us)
  t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
  t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
  t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
  t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
  t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
  t8=24*(t81+t82)
  t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
  t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
  t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
  t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
  t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
  t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
  t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
  t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
  t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
  t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
  t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
  t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
  t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
  t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
  t20=-(n-2)*(t201+t202+t203)
  temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
  temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
  mom3=temp31/temp32
  skewness.null= (mom3-3*mean.null*variance.null-mean.null^3)/variance.null^1.5

  return(list(expection_null=mean.null, variance_null=variance.null, ske_null=skewness.null))
}

#' Main function
#' @description Main function for calculating EMANOVA p-value. Detail description of method can be found in
#' paper "Ensemble test for microbiome data“(currently under review).
#' @param multi_d_mat Input three dimensional distance matrix calculate based on microbiome OTU data.
#' The weighted UniFrac distance and unweighted UniFrac distance distance can be calculated by GUniFrac function in package GUniFrac.
#' The Bray-Curtis distance can be calculated by vegdist function in package vegan. See examples for an illustration.
#' The dimension of the multi_d_mat should be in n times n times number of distance used,
#' where n represents number of samples.
#' @param predictor Predictor variable could be binary(0-1) or continous.
#' @param confonding_stat Variables indicate whether confonding variables should be included
#' @param confonding_var Confonding variables, used only when confonding_stat=TRUE
#' @param r_vec Parameters vector of r values, suggested using value 0.125,0.25,0.5,1,2.
#' @return P-value calculated by EMANOVA
#' @import GUniFrac
#' @import PearsonDS
#' @importFrom vegan vegdist
#' @importFrom abind abind
#' @import fBasics
#' @import MASS
#' @example inst/examples/test_examples.R
#' @export

EMANOVA<-function(multi_d_mat,predictor,confonding_stat=FALSE,confonding_var=NA,r_vec){
  number_of_distances<-length(multi_d_mat[1,1,])
  number_of_r<-length(r_vec)
  p_value_mat<-matrix(NA,nrow = number_of_distances, ncol = number_of_r)
  gamma_hat_mat<-matrix(NA,nrow = number_of_distances, ncol = number_of_r)
  for (i in 1:number_of_distances) {
    for (j in 1:number_of_r) {
      result<-F_stat(multi_d_mat[,,i],predictor,confonding_stat=confonding_stat,confonding_var=confonding_var,r_vec[j])
      moments_result<-permutation_null(result$H_mat,result$G_mat_r)
      gamma_hat<-as.numeric(moments_result$ske_null)
      k_hat<-4/(gamma_hat^2)
      scaled_sample<-(result$F_stat-moments_result$expection_null)/sqrt(moments_result$variance_null)
      if(gamma_hat>0){
        p_value_mat[i,j]<-ppearsonIII(q=scaled_sample,shape=k_hat,location=-sqrt(k_hat),scale=1/sqrt(k_hat),lower.tail = FALSE)
      }
      else{
        p_value_mat[i,j]<-ppearsonIII(q=scaled_sample,shape=k_hat,location=sqrt(k_hat),scale=-1/sqrt(k_hat),lower.tail = FALSE)
      }
      gamma_hat_mat[i,j]<-gamma_hat
    }
  }
  cauchy_statistic<-(1/(number_of_distances*number_of_r))*sum(tan((0.5-p_value_mat)*pi))
  pvalue<-0.5-(atan(cauchy_statistic)/pi)
  return(pvalue)
}

