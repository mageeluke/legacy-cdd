neighbours_2<-function(dt,r,f,alpha,beta){ ## r = radius 
  
  # matrixes are a lot faster than data.frames so we create one to hold our results from loop
  temp<- NULL
  no.col <- 4   #  is the number of index values
  temp.mat <- matrix( nrow=nrow(f) , ncol=no.col)
  
  for(i in 1:nrow(f)){
    
    # data.row <- dt[i,]
    # spp <- data.row$sp
    # focal.tag <- data.row$tag
    # x.row <- data.row$x
    # y.row <- data.row$y
    
    data.row <- f[i,]
    spp <- data.row$sp
    focal.tag <- data.row$tag
    x.row <- data.row$x
    y.row <- data.row$y
    status <-data.row$status
    ba <- data.row$ba
    
    
    square <- dt[x <= x.row+r & x >= x.row-r & y <= y.row+r & y >= y.row-r,]
    
    square$dist <- sqrt(( ( square$x - x.row )^2 + ( square$y - y.row )^2))   # calculate distance between neighbors and focal tree  
    # this part uses data.table syntax notice the generic and familial references
    neighb <- square[dist<= r & square$tag!=focal.tag,]   # remove focal tree and only keep individuals that are with "radius" meters of the focal tree 
    
    neighb[,ch:= ifelse(sp==spp, "c","h")]  # this uses data.table notation
    neighb[,st:= ifelse(status=="A", "A","D")]
    
    # cBA <- neighb[ch=="c" & st=="A", sum(ba, na.rm = T)]  #calculate basal area of conspecific neighbors
    # 
    # hBA <- neighb[ch=="h" & st=="A", sum(ba, na.rm=T)]    # calculate basal area of heterospecific neighbors
    # 
    # lcBA <- neighb[ch=="c" & st=="D", sum(ba, na.rm=T)]
    # 
    # lhBA <- neighb[ch=="h" & st=="D", sum(ba, na.rm=T)]
    #tBA <- sum(neighb$ba, na.rm=T)                                 # total basal area
    
    # nC <- neighb[ch == "c", .N]  # data.table for count rows grouped by condition
    # 
    # nH <- neighb[ch == "h", .N]  
    # 
    # nT <- nrow(neighb)
    
    idw_cons <- sum( ifelse(neighb$ch=="c" & neighb$st=="A", neighb$ba^alpha/neighb$dist^beta, 0), na.rm = T) # distance weighted
    
    idw_hets <- sum( ifelse(neighb$ch=="h" & neighb$st=="A", neighb$ba^alpha/neighb$dist^beta, 0), na.rm = T)
    
    idw_llc <- sum( ifelse(neighb$ch=="c" & neighb$st=="D", neighb$ba^alpha/neighb$dist^beta, 0), na.rm = T)
    
    idw_lhc <- sum( ifelse(neighb$ch=="h" & neighb$st=="D", neighb$ba^alpha/neighb$dist^beta, 0), na.rm = T)
    # distance weighted
    
    #idw_tots <- sum(neighb$ba/neighb$dist, na.rm = T)                             # distance weighted
    
    temp.mat[i, ] <- c(idw_cons, idw_hets, idw_llc,idw_lhc) # put them in the matrix
    
  }
  
  return(temp.mat) # return the matrix
  
}
