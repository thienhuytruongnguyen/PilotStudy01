### Markov model first order

#observed rainfall
obs <- RainDat

#fit MC model
MCparam <- fitMCModel(obs,threshold = 0)
