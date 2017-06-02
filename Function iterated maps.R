# **************************************
# FUNCTIONS FOR GAMES PLAYED ON GRAPHS
# **************************************

#function the payoff matrix
get.score.matrix <- function(T=9,R=5,P=1,S=0){
  # create score matrix
  score.matrix <- matrix(0,nr=2,nc=2)
  rownames(score.matrix) <- colnames(score.matrix) <- c("cooperate","defect")
  score.matrix[1,1] <- R
  score.matrix[2,2] <- P
  score.matrix[1,2] <- S
  score.matrix[2,1] <- T
  return(score.matrix)
}

get.score.matrix2 <- function(T=9,R=5,P=1,S=0){
  # create score matrix
  score.matrix <- matrix(0,nr=2,nc=2)
  rownames(score.matrix) <- colnames(score.matrix) <- c("cooperate","defect")
  score.matrix[1,1] <- R
  score.matrix[2,2] <- P
  score.matrix[1,2] <- S
  score.matrix[2,1] <- T
  return(score.matrix)
}

#Compute score after interaction
compute.scores <- function(action.1,action.2,score.matrix){
  # if both bots decide to perform the same action
  score.1 <- score.matrix[action.1,action.2]
  score.2 <- score.matrix[action.2,action.1]
  return(c(score.1,score.2))
}

#Function of myopic response updating. Checks for each node if the payoff would be higher
#if they would have choosen the other stratergy. If yes, that becomes the new stratergy 
#for the next round, otherwise the node keeps the stratergy.
myopic <- function(action, total.score.matrix){
  #repeating actions with different startergy in order to calculate myopic best resonce dynamics
  action3 <-matrix(rep(NA,nrow(Network)))
  for(i in 1:nrow(Network)){
    action2 <- action
    #Change every iteration one of the stratergies
    action2[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    
    #Compute score wirh different stratergy
    score.myopic2 <- matrix(c(rep(NA, 2*nrow(Edgelist))), ncol=2)
    for(j in 1: nrow(Edgelist)){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix)
    }
    
    #Calculate mean edge score for the other stratergy
    In <- which(Edgelist[,1]==i)
    SIn <- matrix(c(rep(0,length(In))))
    for(j in 1: length(In)){
      SIn[j] <-  sum(score.myopic2[In[j],1])
    }
    scoremyopic <- sum(SIn)/length(SIn)
    
    #Compaire score.matrix.myopic to total.score.matrix to evaluate 
    #behavior and save the best response
    if(scoremyopic > total.score.matrix[i]){
      action3[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    }else{
      action3[i] <- action[i]
    }
  }
  action <- action3[]
  return(action)
}


#Function of unconditional imitation updating. Checks for each neighboring node the payoff, 
#if the payoff of a neigboring node is larger than the node itself, thanthe stratergy of the 
#neighbor node becomes the new stratergy for the next round, otherwise the node keeps the stratergy.
Unconditinalimitation <- function(action, total.score.matrix){
  action3 <-matrix(rep(NA,nrow(Network)))
  #checking neighboring nodes
  for(i in 1:nrow(Network)){
    RW <- which(Edgelist[,1]==i)
    MR <- matrix(c(rep(NA, length(RW))))
    for(j in 1:length(RW)){
      MR[j] <- Edgelist[RW[j],2]
    }
    #Determine highest payoff of neigbor node
    PR <- which.max(total.score.matrix[MR])
    #Accept neighbor stratergy if payoff of the neighbor is higher than the own payoff
    if(total.score.matrix[MR[PR]] > total.score.matrix[i]){
      action3[i] <- action[MR[PR]]
      #Otherwise keep current stratergy
    }else{
      action3[i] <- action[i]
    }
  }
  action <- action3[]
  return(action)
}


#Updating propensities with my opic best replay stratergies.
myopicprop <- function(action, total.score.matrix, prop){
  
  action3 <-matrix(rep(NA,nrow(Network)))
  prop2 <-matrix(rep(NA,nrow(Network)))
  for(i in 1:nrow(Network)){
    action2 <- action
    
    action2[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    
    #Compute score wirh different stratergy
    score.myopic2 <- matrix(c(rep(NA, 2*nrow(Edgelist))), ncol=2)
    for(j in 1: nrow(Edgelist)){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix)
    }
    
    #Calculate mean edge score for the other stratergy
    In <- which(Edgelist[,1]==i)
    SIn <- matrix(c(rep(0,length(In))))
    for(j in 1: length(In)){
      SIn[j] <-  sum(score.myopic2[In[j],1])
    }
    scoremyopic <- sum(SIn)/length(SIn)
    
    #Compaire score.matrix.myopic to total.score.matrix to evaluate 
    #behavior and update propensity score
    if(scoremyopic > total.score.matrix[i]){
      prop2[i] <- ifelse(action[i]=="cooperate", -0.1, 0.1)
      prop[i] <- prop[i] + prop2[i]
      if(prop[i] > 1){prop[i] <-1} #bounds
      if(prop[i] < 0){prop[i] <-0} #bounds
      action3[i] <- sample(c("cooperate","defect"),1, prob= c(prop[i],(1-prop[i])))
    }else{
      prop2[i] <- ifelse(action[i]=="cooperate", 0.1, -0.1)
      prop[i] <- prop[i] + prop2[i]
      if(prop[i] > 1){prop[i] <-1} #bounds
      if(prop[i] < 0){prop[i] <-0} #bounds
      action3[i] <- sample(c("cooperate","defect"),1, prob= c(prop[i],(1-prop[i])))
    }
  }
  action <- action3[]
  
  return(list(action, prop))
}

#Function of myopic response updating. Checks for each node if the payoff would be higher
#if they would have choosen the other stratergy. If yes, that becomes the new stratergy 
#for the next round, otherwise the node keeps the stratergy.
myopicT <- function(action, total.score.matrix){
  #repeating actions with different startergy in order to calculate myopic best resonce dynamics
  action3 <-matrix(rep(NA,nrow(Network)))
  for(i in 1:nrow(Network)){
    action2 <- action
    #Change every iteration one of the stratergies
    action2[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    
    #Compute score wirh different stratergy
    score.myopic2 <- matrix(c(rep(NA, 2*nrow(Edgelist))), ncol=2)
    for(j in 1: nrow(Edgelist)){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix)
    }
    
    #Calculate mean edge score for the other stratergy
    In <- which(Edgelist[,1]==i)
    if(length(In) ==0){
      action3[i] <- action[i]
    } else{
      SIn <- matrix(c(rep(0,length(In))))
      for(j in 1: length(In)){
        SIn[j] <-  sum(score.myopic2[In[j],1])
      }
      scoremyopic <- sum(SIn)/length(SIn)
      
      #Compaire score.matrix.myopic to total.score.matrix to evaluate 
      #behavior and save the best response
      if(is.na(total.score.matrix[i])){
        action3[i] <- action[i]
      }else
        if(scoremyopic > total.score.matrix[i]){
          action3[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
        }else{
          action3[i] <- action[i]
        }
    }}
  action <- action3[]
  return(action)
}

#Updating propensities with my opic best replay stratergies.
myopicpropBoS <- function(action, total.score.matrix, prop){
  
  action3 <-matrix(rep(NA,nrow(Network)))
  prop2 <-matrix(rep(NA,nrow(Network)))
  for(i in 1:nrow(Network)){
    action2 <- action
    
    action2[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    
    #Compute score wirh different stratergy
    score.myopic2 <- matrix(c(rep(NA, 2*nrow(Edgelist))), ncol=2)
    
    for(j in ALPHA1){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix)
    }
    for(j in BETA1){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix2)
    }
    #Calculate mean edge score for the other stratergy
    In <- which(Edgelist[,1]==i)
    SIn <- matrix(c(rep(0,length(In))))
    for(j in 1: length(In)){
      SIn[j] <-  sum(score.myopic2[In[j],1])
    }
    scoremyopic <- sum(SIn)/length(SIn)
    
    #Compaire score.matrix.myopic to total.score.matrix to evaluate 
    #behavior and update propensity score
    if(scoremyopic > total.score.matrix[i]){
      prop2[i] <- ifelse(action[i]=="cooperate", -0.1, 0.1)
      prop[i] <- prop[i] + prop2[i]
      if(prop[i] > 1){prop[i] <-1} #bounds
      if(prop[i] < 0){prop[i] <-0} #bounds
      action3[i] <- sample(c("cooperate","defect"),1, prob= c(prop[i],(1-prop[i])))
    }else{
      prop2[i] <- ifelse(action[i]=="cooperate", 0.1, -0.1)
      prop[i] <- prop[i] + prop2[i]
      if(prop[i] > 1){prop[i] <-1} #bounds
      if(prop[i] < 0){prop[i] <-0} #bounds
      action3[i] <- sample(c("cooperate","defect"),1, prob= c(prop[i],(1-prop[i])))
    }
  }
  action <- action3[]
  
  return(list(action, prop))
}


#Function of myopic response updating. Checks for each node if the payoff would be higher
#if they would have choosen the other stratergy. If yes, that becomes the new stratergy 
#for the next round, otherwise the node keeps the stratergy.
myopicBoS <- function(action, total.score.matrix){
  #repeating actions with different startergy in order to calculate myopic best resonce dynamics
  action3 <-matrix(rep(NA,nrow(Network)))
  for(i in 1:nrow(Network)){
    action2 <- action
    #Change every iteration one of the stratergies
    action2[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    
    #Compute score wirh different stratergy
    score.myopic2 <- matrix(c(rep(NA, 2*nrow(Edgelist))), ncol=2)
    for(j in ALPHA1){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix)
    }
    for(j in BETA1){
      score.myopic2[j,] <- compute.scores(action2[Edgelist[j,1]],action2[Edgelist[j,2]],score.matrix2)
    }
    
    #Calculate mean edge score for the other stratergy
    In <- which(Edgelist[,1]==i)
    SIn <- matrix(c(rep(0,length(In))))
    for(j in 1: length(In)){
      SIn[j] <-  sum(score.myopic2[In[j],1])
    }
    scoremyopic <- sum(SIn)/length(SIn)
    
    #Compaire score.matrix.myopic to total.score.matrix to evaluate 
    #behavior and save the best response
    if(scoremyopic > total.score.matrix[i]){
      action3[i] <- ifelse(action[i]=="cooperate", "defect", "cooperate")
    }else{
      action3[i] <- action[i]
    }
  }
  action <- action3[]
  return(action)
}