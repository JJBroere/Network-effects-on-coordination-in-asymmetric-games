# **************************************
# SIMULATION
# **************************************

# load the libraries
library(igraph)
library(permute)
library(ggplot2)

#load script with functions
setwd("")
source("Functions iterated maps.R")

for(m in 1:1000){
  iter=60
 
  g <- barabasi.game(20,directed = F)
  
  #Utility matrix for different groups
  score.matrix <- get.score.matrix(0,1,0.7,0)
  score.matrix2 <- get.score.matrix(0,0.7,1,0)
  
  #choose update rule Function in file: "Functions iterated maps.R"
  #updaterule <- "myopic"
  #updaterule <- "Unconditinalimitation"
  #updaterule <- "myopicprop"
  #updaterule <- "myopicBoS"
  updaterule<- "myopicpropBoS2"
  
  Network <- get.adjacency(g, type= "both", edges=F)
  BB<-get.edgelist(g)
  CC <- BB[,c(2,1)]
  Edgelist <- rbind(BB,CC)
  
  #Clustering algorithm
  fc <- walktrap.community(as.undirected(g))
  
  # n loop for different starting positions
  for(n in 1:100){
    
    #Save network characteristics
    netid <- c(rep(m,nrow(Network)))
    nodeid<- c(1:20)
    ModG <- c(rep(modularity(fc),nrow(Network)))
    ModL <- fc$modularity
    
    EV <- evcent(g, directed=F)$vector
    Bet <- betweenness(g)
    Clos <- closeness(g)
    Deg <-degree(g)
    TTL <-transitivity(g, type= "local")
    TTG <-transitivity(g, type= "global")
    AVpath <- c(rep(average.path.length(g),nrow(Network)))
    DensN <- rep(graph.density(g), nrow(Network))
    CenN <-rep(centralization.degree(g)$centralization, nrow(Network))
    SegN <-rep(sum(shortest.paths(g) >3)/380, nrow(Network))
    
    #Genrate and randomize different initial configurations
    Initial <-rep(c("cooperate","defect"), nrow(Network)/2)
    SufK <- shuffle(nrow(Network))
    Permutations <- Initial[c(SufK)]
    ALPHA<-which(Permutations =="cooperate")
    BETA<- which(Permutations =="defect")
    ALPHA1<- which(Edgelist[,1] %in% ALPHA)
    BETA1<- which(Edgelist[,1] %in% BETA)
    total <- matrix(c(rep(NA,nrow(Network))), ncol=nrow(Network))
    prop <- Permutations
    prop[prop == "cooperate"] <- 1
    prop[prop == "defect"] <- 0
    prop <- as.numeric(prop)
    prop <- matrix(c(rep(prop, each = iter+1 )), ncol=nrow(Network))
    
    #Objects for storage
    action <- matrix(c(rep(0, nrow(Network)*(iter+1 ))), ncol=nrow(Network))
    action[1,] <- Permutations
    total.score.matrix <- matrix(c(rep(NA, nrow(Network)*iter)), nrow=nrow(Network))
    score.array <- array(c(rep(NA, 2*nrow(Edgelist)*iter)), dim= c( nrow(Edgelist),2,iter))
    
    #Computing scores for every iteration
    for(k in 1:iter){
      
      #computing scores for every edge tie
      #Functions in file "Functions iterated maps", should be in directery
      for(j in ALPHA1){
        score.array[j,,k] <- compute.scores(action[k,Edgelist[j,1]],action[k,Edgelist[j,2]],score.matrix)
      }
      for(j in BETA1){
        score.array[j,,k] <- compute.scores(action[k,Edgelist[j,1]],action[k,Edgelist[j,2]],score.matrix2)
      }
      
      #Computing average tie payoff for every node
      for(i in 1:nrow(Network)){
        In <- which(Edgelist[,1]==i)
        SIn <- matrix(c(rep(0,length(In))))
        for(l in 1: length(In)){
          SIn[l] <-  sum(score.array[In[l],1,k])
          
        }
        total.score.matrix[i,k] <- sum(SIn)/length(SIn)
      }
      
      #Different update rules to evaluatie utility stratergy
      #Functions in file "Functions iterated maps", should be in directery
      if(updaterule == "myopic"){
        action[k+1,] <- myopic(action[k,], t(total.score.matrix[,k]))
      }else if(updaterule == "myopicBoS"){
        action[k+1,] <- myopicBoS(action[k,], t(total.score.matrix[,k]))
      } else if(updaterule == "Unconditinalimitation"){
        action[k+1,] <- Unconditinalimitation(action[k,], t(total.score.matrix[,k]))
      } else if(updaterule == "myopicprop"){
        LAP <- myopicprop(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      } else if(updaterule == "myopicpropBoS"){
        LAP <- myopicpropBoS(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      } else if(updaterule == "myopicpropBoS2"){
        LAP <- myopicpropBoS2(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      }
      else{print("no update rule specified")}
    }
    total[1,]<- rowMeans(total.score.matrix)
    
    #Assesing convergence
    Conver <- rep(as.numeric(identical(action[iter,], action[iter+1,])),nrow(Network))
    
    Colour <- action
    Colour[Colour == "cooperate"] <- 1
    Colour[Colour == "defect"] <- 0
    
    dismat <-distances(g)
    maxdis<- c(rep(max(distances(g)), nrow(Network)))
    #Assasing for each node if it ended in the perferd equlibrium
    Change <- abs(as.numeric(Colour[1,])- as.numeric(Colour[iter+1,]))
    PerChange <- rep(sum(Change/length(Change)), nrow(Network))
    Perco <- sum(as.numeric(Colour[iter+1,]))/length(Colour[iter+1,])
    Percoop <- c(rep(Perco, nrow(Network)))
    #Geometric node characteristics
    Clusterp <- c(rep(1, nrow(Network)))
    inclusdeg <- c(rep(1, nrow(Network)))
    outclusdeg <- c(rep(1, nrow(Network)))
    PerSam <- c(rep(1, nrow(Network)))
    Leng1 <- c(rep(NA, nrow(Network)))
    Leng2 <- c(rep(NA, nrow(Network)))
    Leng3 <- c(rep(NA, nrow(Network)))
    Leng4 <- c(rep(NA, nrow(Network)))
    Leng5 <- c(rep(NA, nrow(Network)))
    Leng6 <- c(rep(NA, nrow(Network)))
    Leng7 <- c(rep(NA, nrow(Network)))
    Leng8 <- c(rep(NA, nrow(Network)))
    Leng9 <- c(rep(NA, nrow(Network)))
    Leng10 <- c(rep(NA, nrow(Network)))
    Leng11 <- c(rep(NA, nrow(Network)))
    Leng12 <- c(rep(NA, nrow(Network)))
    Lengc1 <- c(rep(NA, nrow(Network)))
    Lengc2 <- c(rep(NA, nrow(Network)))
    Lengc3 <- c(rep(NA, nrow(Network)))
    Lengc4 <- c(rep(NA, nrow(Network)))
    Lengc5 <- c(rep(NA, nrow(Network)))
    Lengc6 <- c(rep(NA, nrow(Network)))
    Lengc7 <- c(rep(NA, nrow(Network)))
    Percoopclus <- c(rep(NA, nrow(Network)))
    for(i in 1:nrow(Network)){
      IND <- which(fc$membership[i] == fc$membership,arr.ind=TRUE)
      SAME <- ifelse(Colour[1,IND]==Colour[1,i], 1, 0)
      Percoopclus[i]<-sum(as.numeric(Colour[iter+1,IND]))/length(Colour[iter+1,IND])
      if(sum(as.numeric(SAME))==1){
        Clusterp[i]<-0
      } else{
        Clusterp[i] <- (sum(as.numeric(SAME))-1)/(length(SAME)-1)}
      
      INDE <- which(Edgelist[,1] == i,arr.ind=TRUE)
      INDER<-Edgelist[INDE,2]
      Sam <- ifelse(Colour[1,INDER]==Colour[1,i], 1, 0)
      PerSam[i] <- sum(Sam)/ length(Sam)
      
      inclusdeg[i] <- sum(as.numeric(Edgelist[INDE,2]%in%IND))/ length(as.numeric(Edgelist[INDE,2]%in%IND))
      outclusdeg[i] <- 1-inclusdeg[i]
      
      #calculating length correlation
      
      for(j in 1:max(dismat[i,])){
        Len <-  which(as.numeric(dismat[i,]) == j,arr.ind=TRUE)
        Samlen <- ifelse(as.numeric(Colour[(iter+1),Len])==as.numeric(Colour[(iter+1),i]), 1, -1)
        if(j==1){
          Leng1[i]<- sum(Samlen)/length(Samlen)
        }else
          if(j==2){
            Leng2[i]<- sum(Samlen)/length(Samlen)
          }else
            if(j==3){
              Leng3[i]<- sum(Samlen)/length(Samlen)
            }else
              if(j==4){
                Leng4[i]<- sum(Samlen)/length(Samlen)
              }else
                if(j==5){
                  Leng5[i]<- sum(Samlen)/length(Samlen)
                }else
                  if(j==6){
                    Leng6[i]<- sum(Samlen)/length(Samlen)
                  }else
                    if(j==7){
                      Leng7[i]<- sum(Samlen)/length(Samlen)
                    }else
                      if(j==8){
                        Leng8[i]<- sum(Samlen)/length(Samlen)
                      }else
                        if(j==9){
                          Leng9[i]<- sum(Samlen)/length(Samlen)
                        }else
                          if(j==10){
                            Leng10[i]<- sum(Samlen)/length(Samlen)
                          }else
                            if(j==11){
                              Leng11[i]<- sum(Samlen)/length(Samlen)
                            }else
                              if(j==12){
                                Leng12[i]<- sum(Samlen)/length(Samlen)
                              }else{break}
        
        
      }
      for(j in 1:max(dismat[i,])){
        Len <-  which(as.numeric(dismat[i,]) == j,arr.ind=TRUE)
        Blu<- Len%in%IND
        Samlen <- ifelse(as.numeric(Colour[(iter+1),Len[Blu]])==as.numeric(Colour[(iter+1),i]), 1, -1)
        if(j==1){
          Lengc1[i]<- sum(Samlen)/length(Samlen)
        }else
          if(j==2){
            Lengc2[i]<- sum(Samlen)/length(Samlen)
          }else
            if(j==3){
              Lengc3[i]<- sum(Samlen)/length(Samlen)
            }else
              if(j==4){
                Lengc4[i]<- sum(Samlen)/length(Samlen)
              }else
                if(j==5){
                  Lengc5[i]<- sum(Samlen)/length(Samlen)
                }else
                  if(j==6){
                    Lengc6[i]<- sum(Samlen)/length(Samlen)
                  }else
                    if(j==7){
                      Lengc7[i]<- sum(Samlen)/length(Samlen)
                    }else{break}
        
        
      }
    }
    
    #Combining data    
    dataf <- cbind(netid,nodeid, Change, EV, Bet, Deg, Clusterp,PerSam, inclusdeg, outclusdeg,  
                   TTL, TTG, AVpath,Clos, DensN, CenN, SegN, Conver, Percoop,Percoopclus, Leng1, Leng2, Leng3,
                   Leng4, Leng5, Leng6,Leng7, Leng8, Leng9,Leng10, Leng11, Leng12, maxdis,
                   Lengc1, Lengc2, Lengc3, Lengc4, Lengc5, Lengc6,Lengc7, ModG, ModL)
    if(n ==1){dat <- dataf}
    else{dat <- rbind(dat, dataf)}
  }
  if(m ==1){data <- dat}
  else{data <- rbind(data, dat)}
  print(m)
}

save(data, file="dataScaleFns7.Rda")

for(m in 1:1000){
  iter=60
 
  g <- barabasi.game(20,directed = F)
  
  
  
  #Utility matrix for different groups
  score.matrix <- get.score.matrix(0,1,0.5,0)
  score.matrix2 <- get.score.matrix(0,0.5,1,0)
  
  #choose update rule Function in file: "Functions iterated maps.R"
  #updaterule <- "myopic"
  #updaterule <- "Unconditinalimitation"
  #updaterule <- "myopicprop"
  #updaterule <- "myopicBoS"
  updaterule<- "myopicpropBoS2"
  
  Network <- get.adjacency(g, type= "both", edges=F)
  BB<-get.edgelist(g)
  CC <- BB[,c(2,1)]
  Edgelist <- rbind(BB,CC)
  
  #Clustering algorithm
  fc <- walktrap.community(as.undirected(g))
  
  # n loop for different starting positions
  for(n in 1:100){
    
    #Save network characteristics
    netid <- c(rep(m,nrow(Network)))
    nodeid<- c(1:20)
    ModG <- c(rep(modularity(fc),nrow(Network)))
    ModL <- fc$modularity
    
    EV <- evcent(g, directed=F)$vector
    Bet <- betweenness(g)
    Clos <- closeness(g)
    Deg <-degree(g)
    TTL <-transitivity(g, type= "local")
    TTG <-transitivity(g, type= "global")
    AVpath <- c(rep(average.path.length(g),nrow(Network)))
    DensN <- rep(graph.density(g), nrow(Network))
    CenN <-rep(centralization.degree(g)$centralization, nrow(Network))
    SegN <-rep(sum(shortest.paths(g) >3)/380, nrow(Network))
    
    #Genrate and randomize different initial configurations
    Initial <-rep(c("cooperate","defect"), nrow(Network)/2)
    SufK <- shuffle(nrow(Network))
    Permutations <- Initial[c(SufK)]
    ALPHA<-which(Permutations =="cooperate")
    BETA<- which(Permutations =="defect")
    ALPHA1<- which(Edgelist[,1] %in% ALPHA)
    BETA1<- which(Edgelist[,1] %in% BETA)
    total <- matrix(c(rep(NA,nrow(Network))), ncol=nrow(Network))
    prop <- Permutations
    prop[prop == "cooperate"] <- 1
    prop[prop == "defect"] <- 0
    prop <- as.numeric(prop)
    prop <- matrix(c(rep(prop, each = iter+1 )), ncol=nrow(Network))
    #Initialize starting probabilities
    #if(updaterule == "myopicpropBoS"){ 
    #  prop <- matrix(c(rep(0.5, nrow(Network)*(iter+1 ))), ncol=nrow(Network))
    #}
    #if(updaterule == "myopicpropBoS2"){ 
    # prop <- matrix(c(rep(0.5, nrow(Network)*(iter+1 ))), ncol=nrow(Network))
    #}
    
    #Objects for storage
    action <- matrix(c(rep(0, nrow(Network)*(iter+1 ))), ncol=nrow(Network))
    action[1,] <- Permutations
    total.score.matrix <- matrix(c(rep(NA, nrow(Network)*iter)), nrow=nrow(Network))
    score.array <- array(c(rep(NA, 2*nrow(Edgelist)*iter)), dim= c( nrow(Edgelist),2,iter))
    
    #Computing scores for every iteration
    for(k in 1:iter){
      
      #computing scores for every edge tie
      #Functions in file "Functions iterated maps", should be in directery
      for(j in ALPHA1){
        score.array[j,,k] <- compute.scores(action[k,Edgelist[j,1]],action[k,Edgelist[j,2]],score.matrix)
      }
      for(j in BETA1){
        score.array[j,,k] <- compute.scores(action[k,Edgelist[j,1]],action[k,Edgelist[j,2]],score.matrix2)
      }
      
      #Computing average tie payoff for every node
      for(i in 1:nrow(Network)){
        In <- which(Edgelist[,1]==i)
        SIn <- matrix(c(rep(0,length(In))))
        for(l in 1: length(In)){
          SIn[l] <-  sum(score.array[In[l],1,k])
          
        }
        total.score.matrix[i,k] <- sum(SIn)/length(SIn)
      }
      
      #Different update rules to evaluatie utility stratergy
      #Functions in file "Functions iterated maps", should be in directery
      if(updaterule == "myopic"){
        action[k+1,] <- myopic(action[k,], t(total.score.matrix[,k]))
      }else if(updaterule == "myopicBoS"){
        action[k+1,] <- myopicBoS(action[k,], t(total.score.matrix[,k]))
      } else if(updaterule == "Unconditinalimitation"){
        action[k+1,] <- Unconditinalimitation(action[k,], t(total.score.matrix[,k]))
      } else if(updaterule == "myopicprop"){
        LAP <- myopicprop(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      } else if(updaterule == "myopicpropBoS"){
        LAP <- myopicpropBoS(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      } else if(updaterule == "myopicpropBoS2"){
        LAP <- myopicpropBoS2(action[k,], t(total.score.matrix[,k]), prop[k,])
        action[k+1,] <- LAP[[1]]
        prop[k+1,] <- LAP[[2]]
      }
      else{print("no update rule specified")}
    }
    total[1,]<- rowMeans(total.score.matrix)
    
    #Assesing convergence
    Conver <- rep(as.numeric(identical(action[iter,], action[iter+1,])),nrow(Network))
    
    Colour <- action
    Colour[Colour == "cooperate"] <- 1
    Colour[Colour == "defect"] <- 0
    
    dismat <-distances(g)
    maxdis<- c(rep(max(distances(g)), nrow(Network)))
    #Assasing for each node if it ended in the perferd equlibrium
    Change <- abs(as.numeric(Colour[1,])- as.numeric(Colour[iter+1,]))
    PerChange <- rep(sum(Change/length(Change)), nrow(Network))
    Perco <- sum(as.numeric(Colour[iter+1,]))/length(Colour[iter+1,])
    Percoop <- c(rep(Perco, nrow(Network)))
    #Geometric node characteristics
    Clusterp <- c(rep(1, nrow(Network)))
    inclusdeg <- c(rep(1, nrow(Network)))
    outclusdeg <- c(rep(1, nrow(Network)))
    PerSam <- c(rep(1, nrow(Network)))
    Leng1 <- c(rep(NA, nrow(Network)))
    Leng2 <- c(rep(NA, nrow(Network)))
    Leng3 <- c(rep(NA, nrow(Network)))
    Leng4 <- c(rep(NA, nrow(Network)))
    Leng5 <- c(rep(NA, nrow(Network)))
    Leng6 <- c(rep(NA, nrow(Network)))
    Leng7 <- c(rep(NA, nrow(Network)))
    Leng8 <- c(rep(NA, nrow(Network)))
    Leng9 <- c(rep(NA, nrow(Network)))
    Leng10 <- c(rep(NA, nrow(Network)))
    Leng11 <- c(rep(NA, nrow(Network)))
    Leng12 <- c(rep(NA, nrow(Network)))
    Lengc1 <- c(rep(NA, nrow(Network)))
    Lengc2 <- c(rep(NA, nrow(Network)))
    Lengc3 <- c(rep(NA, nrow(Network)))
    Lengc4 <- c(rep(NA, nrow(Network)))
    Lengc5 <- c(rep(NA, nrow(Network)))
    Lengc6 <- c(rep(NA, nrow(Network)))
    Lengc7 <- c(rep(NA, nrow(Network)))
    Percoopclus <- c(rep(NA, nrow(Network)))
    for(i in 1:nrow(Network)){
      IND <- which(fc$membership[i] == fc$membership,arr.ind=TRUE)
      SAME <- ifelse(Colour[1,IND]==Colour[1,i], 1, 0)
      Percoopclus[i]<-sum(as.numeric(Colour[iter+1,IND]))/length(Colour[iter+1,IND])
      if(sum(as.numeric(SAME))==1){
        Clusterp[i]<-0
      } else{
        Clusterp[i] <- (sum(as.numeric(SAME))-1)/(length(SAME)-1)}
      
      INDE <- which(Edgelist[,1] == i,arr.ind=TRUE)
      INDER<-Edgelist[INDE,2]
      Sam <- ifelse(Colour[1,INDER]==Colour[1,i], 1, 0)
      PerSam[i] <- sum(Sam)/ length(Sam)
      
      inclusdeg[i] <- sum(as.numeric(Edgelist[INDE,2]%in%IND))/ length(as.numeric(Edgelist[INDE,2]%in%IND))
      outclusdeg[i] <- 1-inclusdeg[i]
      
      #calculating length correlation
      
      for(j in 1:max(dismat[i,])){
        Len <-  which(as.numeric(dismat[i,]) == j,arr.ind=TRUE)
        Samlen <- ifelse(as.numeric(Colour[(iter+1),Len])==as.numeric(Colour[(iter+1),i]), 1, -1)
        if(j==1){
          Leng1[i]<- sum(Samlen)/length(Samlen)
        }else
          if(j==2){
            Leng2[i]<- sum(Samlen)/length(Samlen)
          }else
            if(j==3){
              Leng3[i]<- sum(Samlen)/length(Samlen)
            }else
              if(j==4){
                Leng4[i]<- sum(Samlen)/length(Samlen)
              }else
                if(j==5){
                  Leng5[i]<- sum(Samlen)/length(Samlen)
                }else
                  if(j==6){
                    Leng6[i]<- sum(Samlen)/length(Samlen)
                  }else
                    if(j==7){
                      Leng7[i]<- sum(Samlen)/length(Samlen)
                    }else
                      if(j==8){
                        Leng8[i]<- sum(Samlen)/length(Samlen)
                      }else
                        if(j==9){
                          Leng9[i]<- sum(Samlen)/length(Samlen)
                        }else
                          if(j==10){
                            Leng10[i]<- sum(Samlen)/length(Samlen)
                          }else
                            if(j==11){
                              Leng11[i]<- sum(Samlen)/length(Samlen)
                            }else
                              if(j==12){
                                Leng12[i]<- sum(Samlen)/length(Samlen)
                              }else{break}
        
        
      }
      for(j in 1:max(dismat[i,])){
        Len <-  which(as.numeric(dismat[i,]) == j,arr.ind=TRUE)
        Blu<- Len%in%IND
        Samlen <- ifelse(as.numeric(Colour[(iter+1),Len[Blu]])==as.numeric(Colour[(iter+1),i]), 1, -1)
        if(j==1){
          Lengc1[i]<- sum(Samlen)/length(Samlen)
        }else
          if(j==2){
            Lengc2[i]<- sum(Samlen)/length(Samlen)
          }else
            if(j==3){
              Lengc3[i]<- sum(Samlen)/length(Samlen)
            }else
              if(j==4){
                Lengc4[i]<- sum(Samlen)/length(Samlen)
              }else
                if(j==5){
                  Lengc5[i]<- sum(Samlen)/length(Samlen)
                }else
                  if(j==6){
                    Lengc6[i]<- sum(Samlen)/length(Samlen)
                  }else
                    if(j==7){
                      Lengc7[i]<- sum(Samlen)/length(Samlen)
                    }else{break}
        
        
      }
    }
    
    #Combining data    
    dataf <- cbind(netid,nodeid, Change, EV, Bet, Deg, Clusterp,PerSam, inclusdeg, outclusdeg,  
                   TTL, TTG, AVpath,Clos, DensN, CenN, SegN, Conver, Percoop,Percoopclus, Leng1, Leng2, Leng3,
                   Leng4, Leng5, Leng6,Leng7, Leng8, Leng9,Leng10, Leng11, Leng12, maxdis,
                   Lengc1, Lengc2, Lengc3, Lengc4, Lengc5, Lengc6,Lengc7, ModG, ModL)
    if(n ==1){dat <- dataf}
    else{dat <- rbind(dat, dataf)}
  }
  if(m ==1){data <- dat}
  else{data <- rbind(data, dat)}
  print(m)
}

save(data, file="dataScaleFns5.Rda")