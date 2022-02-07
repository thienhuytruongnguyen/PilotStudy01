#------------------Runoff Model----------------------------------
par(mfrow = c(2,2))
for (i in 1:12){
  ##Set GR4J model parameter
  inputGR4J <-
    makeInputGR4J(
      P = RainDat$Value[indRainDate$i.mm[[i]]],
      Q = FlowDat$Value[indRainDate$i.mm[[i]]],
      E = EvapDat$Value[indRainDate$i.mm[[i]]],
      Date = FlowDat$Date_Time[indRainDate$i.mm[[i]]],
      A = whichSite$area
    )
  
  #Set start and end of period
  start <- FlowDat[indRainDate$i.mm[[i]],1][1] ; end <- FlowDat[indRainDate$i.mm[[i]],1][length(indRainDate$i.mm[[i]])]
  
  #months, set warm-up period, airGR prefer more than 12 months
  warmup <- 24
  
  #make parameter list
  paramGR4J <-
    getParamGR4J_SepMonth(
      inputGR4J = inputGR4J,
      start = start,
      end = end,
      warmup = warmup,
      parameter = "known",
      knownParam = catchmentInfo[1, c(9, 10, 11, 12)]
    )
  
  ##Run GR4J model
  outputGR4J <- runGR4J(paramGR4J)
  virObsFlow <- outputGR4J$Qsim
  
  
  simFlowRep <-
    getSimFlowRep(simRainRep = simRainRep[[1]][indRainDate$i.mm[[i]],], paramGR4J = paramGR4J)
  
  #plot(inputGR4J$Q[32:length(indRainDate$i.mm[[1]])],virObsFlow)
  plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
  
}


SSE_FDC_SepMonth<-
  function(theta,
           indRainDate,
           paramGR4J,
           inputGR4J,
           runOptionGR4J,
           virObsFDC) {
    
    #declare simrain dataframe
    simRainRep <- matrix(NA, nrow = indRainDate$nDy, ncol = rep)
    #Loop for each month
    for (i in 1:12){
      #Loop for each replicate
      for (j in 1:rep){
        #Create the occurence binary series
        #set.seed(68)
        U_t <- runif(length(indRainDate$i.mm[[i]]),0,1)
        bin <- MCmodel_C(length(U_t), occurParam[i,1], occurParam[i,2], U_t)
        #make rain ts from gamma distribution
        randRain <- rgamma(length(bin[bin==1]), amountParam[i,1], amountParam[i,2])
        
        bin[bin==1] <- randRain
        #matching
        simRainRep[indRainDate$i.mm[[i]],j] <- bin
      }
    }
    return(simRainRep)
    
    #Passing element in theta to WGEN parameter
    #Occurence model parameters
    paramMC <- matrix(NA,12,2)
    paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
    #Amount model parameters
    paramAmount <- matrix(NA,12,2)
    paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
    
    #Generate sim rain with given parameters above (a vector)
    simRainRep <- WGEN_V4.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
    #Generate sim flow with sim rain
    #add simRain to paramGR4J options
    inputGR4J[[2]] <- simRainRep[,1]
    #RunGR4J model with updated sim rain
    outputGR4J <- airGR::RunModel_GR4J(InputsModel = inputGR4J, RunOptions = runOptionGR4J, Param = paramGR4J)
    #Get sim flow from output GR4J
    simFlowRep <- outputGR4J$Qsim
    
    #Calculate Exceedance Probability for sim flow and virobs flow
    simFDC <- getExceedProb_V2.0(simFlowRep)
    
    #Calculate the Sum of square Error
    err <- simFDC$flow - virObsFDC
    SSE <- sum(err^2)
    #SSE <- SSE*100
    return(SSE)
  }