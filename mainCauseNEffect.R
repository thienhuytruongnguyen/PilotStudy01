#Get annual maxima of virtual observed flow
annualMaxFlow <- getAnnualMaxima(indObsDate = indFlowDate,
                                 value = virObsFlow)
#locate the index of the annual maxima events within the virtual observed flow
indAnnualMaxFlow <- match(c(annualMaxFlow), virObsFlow)
#Get the date that the virtual observed event happened
dateAnnualMaxFlow <-
  RainDat$Date_Time[paramGR4J[[2]]][indAnnualMaxFlow]
#Get annual maxima of the observed rainfall
annualMaxRain <- getAnnualMaxima(indObsDate = indRainDate,
                                 value = RainDat$Value[paramGR4J[[2]]])
#Locate the index of the annual maxima events within the rainfall timeseries
indAnnualMaxRain <-
  match(c(annualMaxRain), RainDat$Value[paramGR4J[[2]]])
#Get the date that the virtual observed event happened
dateAnnualMaxRain <-
  RainDat$Date_Time[paramGR4J[[2]]][indAnnualMaxRain]

#Get annualmaxima of simRainRep
annualMaxRainRep <- data.frame(matrix(NA,nrow = 55,ncol = ncol(simRainRep)))

for (i in 1:100){
  annualMaxRainRep[,i] <- getAnnualMaxima(indObsDate = indRainDate,
                                      value = simRainRep[paramGR4J[[2]],i])
}
medianSimAnnualMaxima <- apply(annualMaxRainRep,1,percentile50)
for (i in 1:100){
  indSimAnnualMaxRain <- match(c(medianSimAnnualMaxima), simRainRep[,i])
}
indSimAnnualMaxRain <- match(c(medianSimAnnualMaxima), simRainRep[,1])

simAnnualMax <- getAnnualMaxima(indObsDate = indRainDate,
                                value = simRainRep[paramGR4J[[2]],1])
indSimAnnualMaxRain <- match(c(simAnnualMax), simRainRep[paramGR4J[[2]],1])
simDateAnnualMaxRain <-
  RainDat$Date_Time[paramGR4J[[2]]][indSimAnnualMaxRain]
#Plotting

plot(
  dateAnnualMaxFlow[1:15],
  annualMaxFlow[1:15],
  type = "b",
  ylim = c(0, 50),
  xaxt = "n",
  xlab = "",
  ylab = ""
)
axis(1, at = dateAnnualMaxFlow[1:15], labels = FALSE)
text(
  x = dateAnnualMaxFlow[1:15],
  y = par("usr")[3] - 5,
  labels = format(dateAnnualMaxFlow[1:15], "%d %b %y"),
  srt = 70,
  xpd = NA,
  adj = 0.965,
  cex = .7
)
par(new = TRUE)
plot(
  dateAnnualMaxRain[1:15],
  annualMaxRain[1:15],
  type = "b",
  col = "blue",
  ylim = rev(c(0,200)),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)
axis(3, at = dateAnnualMaxRain[1:15], labels = FALSE, col = "blue")
axis(4, seq(0, max(annualMaxRain[1:15]), 20))
text(
  x = dateAnnualMaxRain[1:15],
  y = par("usr")[4] - 35,
  labels = format(dateAnnualMaxRain[1:15], "%d %b %y"),
  srt = 70,
  xpd = NA,
  adj = 0.965,
  cex = .7,
  col = "blue"
)
par(new = TRUE)
plot(
  simDateAnnualMaxRain[1:15],
  simAnnualMax[1:15],
  type = "b",
  col = "red",
  ylim = rev(range(simAnnualMax[1:15])),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)

indVirObsFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]])
virObsFlow74Jul <- virObsFlow[indVirObsFlowDate$i.ym[[10]][[7]]]
virObsFlow74Aug <- virObsFlow[indVirObsFlowDate$i.ym[[10]][[8]]]
virObsFlow74Sep <- virObsFlow[indVirObsFlowDate$i.ym[[10]][[9]]]
inteVOFlow74 <- c(virObsFlow74Jul,virObsFlow74Aug,virObsFlow74Sep)

dateFlow74Jul <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[7]]]
dateFlow74Aug <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[8]]]
dateFlow74Sep <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[9]]]
inteDateFlow74 <- c(dateFlow74Jul,dateFlow74Aug,dateFlow74Sep)

rain74Jul <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[7]]]
rain74Aug <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[8]]]
rain74Sep <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[10]][[9]]]
inteRain74 <- c(rain74Jul,rain74Aug,rain74Sep)

plot(
  inteDateFlow74,
  inteVOFlow74,
  type = "b",
  ylim = c(0, 100),
  ylab = "",
  xlab = "",
  yaxt = "n"
)
axis(2, seq(0, 50, 5))
par(new = TRUE)
plot(
  inteDateFlow74,
  inteRain74,
  type = "b",
  ylim = rev(c(0, 200)),
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  col = "blue"
)
axis(4, seq(0, max(inteRain74), 10))




indVirObsFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]])
virObsFlow75Jul <- virObsFlow[indVirObsFlowDate$i.ym[[11]][[7]]]
virObsFlow75Aug <- virObsFlow[indVirObsFlowDate$i.ym[[11]][[8]]]
virObsFlow75Sep <- virObsFlow[indVirObsFlowDate$i.ym[[11]][[9]]]
virObsFlow75Oct <- virObsFlow[indVirObsFlowDate$i.ym[[11]][[10]]]
virObsFlow75Nov <- virObsFlow[indVirObsFlowDate$i.ym[[11]][[11]]]
inteVOFlow75 <- c(virObsFlow75Jul,virObsFlow75Aug,virObsFlow75Sep,virObsFlow75Oct,virObsFlow75Nov)

dateFlow75Jul <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[7]]]
dateFlow75Aug <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[8]]]
dateFlow75Sep <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[9]]]
dateFlow75Oct <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[10]]]
dateFlow75Nov <- RainDat$Date_Time[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[11]]]
inteDateFlow75 <- c(dateFlow75Jul,dateFlow75Aug,dateFlow75Sep,dateFlow75Oct,dateFlow75Nov)

rain75Jul <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[7]]]
rain75Aug <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[8]]]
rain75Sep <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[9]]]
rain75Oct <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[10]]]
rain75Nov <- RainDat$Value[paramGR4J[[2]]][indVirObsFlowDate$i.ym[[11]][[11]]]
inteRain75 <- c(rain75Jul,rain75Aug,rain75Sep,rain75Oct,rain75Nov)

plot(
  inteDateFlow75,
  inteVOFlow75,
  type = "l",
  ylim = c(0, 100),
  ylab = "",
  xlab = "",
  yaxt = "n"
)
axis(2, seq(0, 50, 5))
par(new = TRUE)
plot(
  inteDateFlow75,
  inteRain75,
  type = "l",
  ylim = rev(c(0, 200)),
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  col = "blue"
)
axis(4, seq(0, max(inteRain75), 10))

virObsFlow78 <- virObsFlow[4623:4867]
dateFlow78 <- RainDat$Date_Time[paramGR4J[[2]]][4623:4867]
rain78 <- RainDat$Value[paramGR4J[[2]]][4623:4867]
plot(
  dateFlow78,
  virObsFlow78,
  type = "b",
  ylim = c(0, 100),
  ylab = "",
  xlab = "",
  yaxt = "n"
)
axis(2, seq(0, 50, 5))
par(new = TRUE)
plot(
  dateFlow78,
  rain78,
  type = "b",
  ylim = rev(c(0, 200)),
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  col = "blue"
)
axis(4, seq(0, max(inteRain75), 10))