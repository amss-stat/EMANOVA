library(GUniFrac)
data("throat.otu.tab")
data("throat.tree")
data("throat.meta")

####transfer the outcome variable into binary.
smoke_stat<-rep(0,60)
for (i in 1:60) {
  if(throat.meta$SmokingStatus[i]==throat.meta$SmokingStatus[1]){
    smoke_stat[i]<-0
  }
  else{
    smoke_stat[i]<-1
  }
}

####Calculate relative abundance.
row.sum<-apply(throat.otu.tab,1,sum)
std_out<-throat.otu.tab/row.sum

####calculated the distance, D.05 is weighted UniFrac distance,
####Du unweighted UniFrac distance, D.BC is Bray-Curtis distance.
unifracs = GUniFrac::GUniFrac(std_out, throat.tree, alpha = c(0,0.5,1))$unifracs
Du  = unifracs[,,"d_UW"]
D.05= unifracs[,,"d_0.5"]
D.BC= as.matrix(vegan::vegdist(std_out, method="bray"))

combine_distance<-abind::abind(D.05, Du, D.BC, along = 3)
EMANOVA(combine_distance, smoke_stat, confonding_stat = FALSE, r_vec = c(0.125,0.25,0.5,1,2))
