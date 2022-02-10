#------------------Runoff Model----------------------------------
optimParamWGEN <- paramWGEN

for (i in 1:12){
  ##Set GR4J model parameter
  inputGR4J_SepMonth <-
    makeInputGR4J(
      P = RainDat$Value[indRainDate$i.mm[[i]]],
      Q = FlowDat$Value[indRainDate$i.mm[[i]]],
      E = EvapDat$Value[indRainDate$i.mm[[i]]],
      Date = FlowDat$Date_Time[indRainDate$i.mm[[i]]],
      A = whichSite$area
    )
  
  #Set start and end of period
  start <- FlowDat[indRainDate$i.mm[[i]],1][1] ; end <- FlowDat[indRainDate$i.mm[[i]],1][length(indRainDate$i.mm[[i]])]

  
  #make parameter list
  paramGR4J_SepMonth <-
    getParamGR4J_SepMonth(
      inputGR4J = inputGR4J_SepMonth,
      start = start,
      end = end,
      parameter = "known",
      knownParam = catchmentInfo[1, c(9, 10, 11, 12)]
    )
  
  ##Run GR4J model
  outputGR4J_SepMonth <- runGR4J(paramGR4J_SepMonth)
  virObsFlowSM <- outputGR4J_SepMonth$Qsim
  virObsFDCSM <- getExceedProb_V2.0(virObsFlowSM)
  
  #Initiate parameter to optimise
  iniTheta <- rep(0,4) 
  iniTheta[1] <- paramWGEN[i,1]; iniTheta[2] <-paramWGEN[i,2];
  iniTheta[3] <- paramWGEN[i,3]; iniTheta[4] <- paramWGEN[i,4]
  lowerTheta <- rep(0,4)
  lowerTheta[1:2] = 0; lowerTheta[3:4] = 1e-10
  upperTheta <- rep(0,4)
  upperTheta[1:2] = 1; upperTheta[3:4] = Inf
  
  #Run SCE optim
  optResult <- hydromad::SCEoptim(FUN = SSE_WeightedFDC_SingleMonth,
                                  par = iniTheta,
                                  obsRain = RainDat[indRainDate$i.mm[[i]],2],
                                  paramGR4J = paramGR4J_SepMonth[[1]],
                                  inputGR4J = paramGR4J_SepMonth[[3]],
                                  runOptionGR4J = paramGR4J_SepMonth[[4]],
                                  virObsFDC = virObsFDCSM,
                                  lower = lowerTheta,
                                  upper = upperTheta)
  
  optimParamWGEN[i,] <- optResult$par
  remove(paramGR4J_SepMonth, inputGR4J_SepMonth)
}

thetaTrialSepmonth <- rep(0,48)
thetaTrialSepmonth[1:12]<-optimParamWGEN[1:12,1]; thetaTrialSepmonth[13:24]<-optimParamWGEN[1:12,2]
thetaTrialSepmonth[25:36]<-optimParamWGEN[1:12,3]; thetaTrialSepmonth[37:48]<-optimParamWGEN[1:12,4]


#Testing----
iniTheta <- rep(0,4); iniTheta[1] <- paramWGEN[8,1]; iniTheta[2] <-paramWGEN[8,2];
iniTheta[3] <- paramWGEN[8,3]; iniTheta[4] <- paramWGEN[8,4]
virObsFDC <- getExceedProb_V2.0(virObsFlow)
SSEsingleMonth <- SSE_FDC_SingleMonth(theta = theta,
                                      obsRain = RainDat[indRainDate$i.mm[[8]],2],
                                      paramGR4J = paramGR4J_SepMonth[[1]],
                                      inputGR4J = paramGR4J_SepMonth[[3]],
                                      runOptionGR4J = paramGR4J_SepMonth[[4]],
                                      virObsFDC = virObsFDC$flow)
#Run optim----

lowerTheta <- rep(0,4)
lowerTheta[1:2] = 0; lowerTheta[3:4] = 1e-10
upperTheta <- rep(0,4)
upperTheta[1:2] = 1; upperTheta[3:4] = Inf

#Run SCE optim
optResult <- hydromad::SCEoptim(FUN = SSE_FDC_SingleMonth,
                                par = iniTheta,
                                obsRain = RainDat[indRainDate$i.mm[[8]],2],
                                paramGR4J = paramGR4J[[1]],
                                inputGR4J = paramGR4J[[3]],
                                runOptionGR4J = paramGR4J[[4]],
                                virObsFDC = virObsFDC$flow,
                                lower = lowerTheta,
                                upper = upperTheta)

#Testing----
theta <- optResult$par
#Passing element in theta to WGEN parameter
occurParam <- vector(length = 2)
occurParam[1] <- theta[1]; occurParam[2] <- theta[2]

amountParam <- vector(length = 2)
amountParam[1] <- theta[3]; amountParam[2] <- theta[4]





simRainRep <-
  WGEN_SingleMonth(
    occurParam = occurParam,
    amountParam = amountParam,
    obsRain = RainDat[indRainDate$i.mm[[8]], 2],
    rep = 100
  )

simFlowRep <-
  getSimFlowRep(simRainRep = simRainRep, paramGR4J = paramGR4J)

#plot(inputGR4J$Q[32:length(indRainDate$i.mm[[1]])],virObsFlow)
plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
