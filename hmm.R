promoters=read.csv("promoters.data", header = F, dec = ",",strip.white = TRUE, stringsAsFactors = FALSE)
promoters[1,]



head(promoters)
#dividing positive and negative set
positive_observations=subset(promoters, V1 == '+', 3)
negative_observations=subset(promoters, V1 == '-', 3)

#adding s and x at front and end
positive_observations=sapply(positive_observations,function(x) paste("S", x, "X", sep=""))
negative_observations=sapply(negative_observations,function(x) paste("S", x, "X", sep=""))

#splitting strings to characters
positive_observations=strsplit(positive_observations, "")
negative_observations=strsplit(negative_observations, "")
head(positive_observations[[1]], n = 15)

states=c("S", "X", "a", "c", "g", "t")
symbols=c("S", "X", "a", "c", "g", "t")
startingProbabilities= c(1,0,0,0,0,0)
emissionProbabilities=diag(6)
colnames(emissionProbabilities)=states
rownames(emissionProbabilities)=symbols
emissionProbabilities
#transmission probabilities
calculateTransitionProbabilities=function(data, states) {
    transitionProbabilities=matrix(0, length(states), length(states))
    colnames(transitionProbabilities)=states
    rownames(transitionProbabilities)=states
    for(index in 1:(length(data) - 1)) {
      current_state=data[index]
      next_state=data[index + 1]
      transitionProbabilities[current_state, next_state]=transitionProbabilities[current_state, next_state] + 1
    }
    transitionProbabilities=sweep(transitionProbabilities,1,rowSums(transitionProbabilities), FUN = "/")
    return(transitionProbabilities)
}

negative_observation=Reduce(function(x, y) c(x, y),negative_observations,c())
(transitionProbabilitiesNeg=calculateTransitionProbabilities(negative_observation, states))

library("HMM")
negative_hmm=initHMM(states,symbols,startProbs=startingProbabilities,transProbs=transitionProbabilitiesNeg,
                          emissionProbs=emissionProbabilities)

incorrect=0
for(obs in 1:length(positive_observations)) {
  positive_observation=Reduce(function(x, y) c(x, y),positive_observations[-obs],c())
  transitionProbabilitiesPos=calculateTransitionProbabilities(positive_observation, states)
  positive_hmm =initHMM(states,symbols,startProbs = startingProbabilities,transProbs = transitionProbabilitiesPos, 
                          emissionProbs = emissionProbabilities)
  
  test_observation=positive_observations[[obs]]
  final_index=length(test_observation)
  
  pos_probs=exp(forward(positive_hmm, test_observation))
  neg_probs=exp(forward(negative_hmm, test_observation))
  pos_seq_prob=sum(pos_probs[, final_index])
  neg_seq_prob=sum(neg_probs[, final_index])
  
  if(pos_seq_prob<neg_seq_prob)incorrect=incorrect + 1
  
}
incorrect


#train with all positives
positive_observation=Reduce(function(x, y) c(x, y),positive_observations,c())
transitionProbabilitiesPos=calculateTransitionProbabilities(positive_observation, states)
positive_hmm=initHMM(states,symbols,startProbs=startingProbabilities, transProbs = transitionProbabilitiesPos, 
                         emissionProbs = emissionProbabilities)

for (obs in 1:length(negative_observations)) {
  
  negative_observation=Reduce(function(x, y) c(x, y),negative_observations[-obs], c())
  transitionProbabilitiesNeg=calculateTransitionProbabilities(negative_observation, states)
  negative_hmm=initHMM(states, symbols, startProbs = startingProbabilities, 
                          transProbs = transitionProbabilitiesNeg, 
                          emissionProbs = emissionProbabilities)
  
  test_observation=negative_observations[[obs]]
  final_index=length(test_observation)
  
  pos_probs=exp(forward(positive_hmm,test_observation))
  neg_probs=exp(forward(negative_hmm,test_observation))
  pos_seq_prob=sum(pos_probs[, final_index])
  neg_seq_prob=sum(neg_probs[, final_index])
  
  if (pos_seq_prob > neg_seq_prob) incorrect=incorrect+1
  
}
incorrect
(cross_validation_accuracy= 1 - (incorrect/nrow(promoters)))