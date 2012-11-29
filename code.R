##################################################################
## Source code for the paper: "Time dependent correlation between AC
## power fluctuations time series from a network of DC/AC inverters
## in a large PV plant"

## Copyright (C) 2012 Oscar Perpiñán Lamigueiro

## This program is free software you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation; either version 2 of the License,
## or (at your option) any later version.
   
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
####################################################################

library(lattice)
library(latticeExtra)
library(parallel)
library(solaR)
library(wmtsa)
library(cluster)
library(car)
library(zoo)

## Change to the folder where the local copy is located. The data
## folder of the GitHub repository contains three files of
## intermediate results: meteo.RData (clustering), wavCorDF.RData and
## corDistLong.RData (wavelet correlation). With these three files
## most of the figures can be reproduced. Due to large size and
## confidential information requirements, the full collection of data
## can not be published. Sample data is available from the authors
## upon request.
setwd('~/Investigacion/DocsPropios/Moura/R')

################################################################################
## FUNCTIONS
################################################################################

## Clean a time series
no0 <- function(x, threshold=0.1){
  idx <- seq_along(x)
  len=length(x)
  mx=max(abs(x), na.rm=TRUE)      ## maximum value
  thr=mx*threshold        ## noise threshold
  xrle <- rle(abs(x)>thr) ## parts of the signal above the threshold
  xrle.len <- xrle$lengths  ##lengths of these parts
  xrle.value <- xrle$values 
  n=length(xrle.value)      ## number of parts

  hd <- 1:xrle.len[1]         ##head of the signal
  tl <- (len-xrle.len[n]):len ##tail of the signal
  
  if(!xrle.value[1] & !xrle.value[n]) {
    idx[-c(hd, tl)] ## head and tail are below the threshold, suppressed.
  } else if(!xrle.value[1] & xrle.value[n]) {
    return(idx[-hd]) ## head below threshold, suppressed.
  } else if(xrle.value[1] & !xrle.value[n]) {
    return(idx[-tl]) ## tail below threshold, suppressed.
  } else return(idx) ## unmodified
}

## wavelet variance
wavVar.uni <- function(x){
    delta <- 1/frequency(x)
    aux <- wavVar(x)
    unbiased.data <- aux$block$unbiased
    scales <- as.numeric(aux$scales)*delta
    conf.low <- aux$confidence$n3$low
    conf.high <- aux$confidence$n3$high
    df <- data.frame(var=unbiased.data, scale=scales, low=conf.low, high=conf.high)
    rownames(df) <- NULL
    df
  }

##Wavelet variance of a multivariate time series (using wavVar.uni)
wavVar.multi <- function(x){
  delta <- 1/frequency(x)
  result0 <- apply(x, 2, wavVar.uni)
  result <- do.call('rbind', result0)
  rownames(result) <- NULL
  result$ID <- rep(names(result0), each=dim(result)[1]/length(result0))
  result
}

## wavelet coefficients
modwt <- function(x, ...){
  y <- wavMODWT(x, ...)
  dat <- do.call(cbind, y$data)
  row.names(dat) <- NULL
  z <- cbind(x, dat)
  z
  }

correlogram <- function(x, scales=list(cex=0.5, x=list(rot=90)), par.settings=BTCTheme, ...){
  ord <- order.dendrogram(as.dendrogram(hclust(dist(x))))
  levelplot(x[ord, ord], scales=scales, xlab='', ylab='',
            par.settings=par.settings, ...)
}

################################################################################
## WAVELET CORRELATION
################################################################################

## Change to the folder where data is available. Sample data is
## available from the authors upon request
old <- setwd('data/')

inverters <- read.csv2('coord_reticulas.csv')

## Distances between inverters
distances <- as.vector(dist(inverters[,2:3]))

## Labels of combinations between inverters
labels <- inverters$reticula
ids <- outer(labels, labels, paste, sep='_')
ids <- ids[lower.tri(ids)] ## drop repetitions

## Angle between inverters (positive radians)
x <- inverters[,2]
xx <- outer(x, x, function(x, y)abs(x-y))
y <- inverters[,3]
yy <- outer(y, y, function(x, y)abs(x-y))
az <- atan2(yy, xx)
az <- as.vector(az[lower.tri(az)])

##change to the sequence of available dates
fechas <- seq(as.Date('2010/06/15'), as.Date('2011/12/15'), by='1 day')

for (d in seq_along(fechas)){
  ## Open the correspondent file
  fecha <- fechas[d]
  print(fecha)
  ips <- paste('192.168.26', 16:18, sep='.')
  fechaFormat <- format(fecha, '%d%m%Y')
  fechaIP <- paste(fechaFormat, ips, sep='_')
  fich <- paste(fechaIP, '.csv', sep='')
  folder <- paste(format(fecha, '%Y'),
                  paste('CSV', format(fecha, '%m'), sep='_'),
                  fechaIP,
                  sep='/')
  URLs<- paste(folder, fich, sep='/')

  try({
    z <- lapply(URLs, function(x) read.zoo(x,
                                           tz='UTC', format='%d/%m/%Y %H:%M:%S',
                                           header=TRUE, sep=','))


    power <- do.call(cbind, lapply(z, function(x){ ## Inverter power
      idCols <- grep('Potencia', names(x))
      variadores <- x[ ,idCols] ##Each inverter contains 4 parts
      labels <- substr(names(variadores), 9, 10)
      groups <- split(1:ncol(variadores), labels)
      invPower <- sapply(groups, ##The total power is the sum of the parts
                         function(cols)rowSums(variadores[,cols]))
      res <- zoo(invPower, index(x))
    }))

    idx0 <- no0(rowMeans(power))
    powerClean <- power[idx0,]
    outlier <- which(powerClean>600) ## | power<=0 )
    powerClean[outlier] <- NA

    ## MODWT for each inverter modwtPower is a list: each component is
    ## a zoo object with the signal and the wavelet coefficients
    modwtPower <- lapply(powerClean, modwt)
    ## Reorder inverters according to inverters$labels to use safely distances and az
    idxSort <- order(names(modwtPower))
    modwtPower <- modwtPower[idxSort]
    delta <- 1/frequency(powerClean)


    nLevels <- 9
    scales <- paste('d', seq_len(nLevels), sep='')

    ## wavDetail is a list: each component is a zoo with the wavelet
    ## coefficients of each scale for each inverter
    wavDetail <- list() 
    for (s in scales){
      wavDetail[[s]] <- do.call(cbind, lapply(modwtPower, function(x)x[,s]))
    }
                    
    wavCor <- lapply(scales, function(s)cor(wavDetail[[s]], use='pairwise.complete.obs')) 
    names(wavCor) <- scales


    ##corDist is a data.frame with the wavelet correlation for each scale
    correlation <- sapply(wavCor, function(x)x[lower.tri(x)])
    corDist <- as.data.frame(cbind(correlation, distances, az))
    corDist$ids <- ids
    corDist$class <- cut(distances, breaks=seq(min(distances), max(distances), length=50))

    corDistLong <- reshape(corDist,
                           varying=list(1:9), v.names='correlation',
                           timevar='scale', times=names(corDist)[1:9],
                           direction='long')

    save(powerClean, modwtPower, wavDetail, wavCor, corDist, corDistLong,
         file=paste(fecha, 'RData', sep='.'))
  })
}

setwd(old)

## Confidence intervals
s <- 1:9 ##scales
at <- seq(0, 1, .1) ##correlation values
N <- dim(powerClean)[1] ##number of samples

g <- expand.grid(at=at, scale=s)

p=0.975

g <- with(g,{
          Nhat <- trunc(N/(2^s))
          q <- qnorm(p)/sqrt(Nhat-3)
          lower <- tanh(atanh(at)-q[scale])
          higher <- tanh(atanh(at)+q[scale])
          int <- abs(higher - lower)
          data.frame(at, scale, Nhat, q, lower, higher, int)
          })

## Figure 2
contourplot(int~scale*at, data=g,
            cut=35, lwd=0.4, labels=list(cex=0.7),
            xlab='Wavelet scale', ylab='Wavelet correlation')

################################################################################
## CLUSTERING
################################################################################
## Latitude of the plant
lat=38.2

## Available dates
fechas <- seq(as.Date('2010/06/15'), as.Date('2011/12/15'), by='1 day')

##The result of this section is available with
## load('data/meteo.RData')

## change to the folder where data is available.
old <- setwd('data/')

meteo <- mclapply(fechas, function(fecha){
  print(fecha)
  ips <- paste('192.168.26', 16:18, sep='.')
  fechaFormat <- format(fecha, '%d%m%Y')
  fechaIP <- paste(fechaFormat, ips, sep='_')
  fich <- paste(fechaIP, '.csv', sep='')
  folder <- paste(format(fecha, '%Y'),
                  paste('CSV', format(fecha, '%m'), sep='_'),
                  fechaIP,
                  sep='/')
  URLs<- paste(folder, fich, sep='/')
  

  res <- NULL
  try({
    z <- lapply(URLs, function(x) read.zoo(x,
                                           tz='UTC', format='%d/%m/%Y %H:%M:%S',
                                           header=TRUE, sep=','))
    z <- do.call(cbind, z)
    rad <- z[,which(names(z)=="Estaciones.EstPiranometros.Rglobal")]
    rad[rad<0] <- NA

    G0d <- sum(rad, na.rm=1)/(3600/5)
    ktd <- G0d/as.numeric(fSolD(lat=lat, index(rad)[1])$Bo0d)

    radWavVar <- wavVar.uni(rad)$var

    wind <- z[,which(names(z)=="Estaciones.EstMetereologica.VelocidadViento")]
    wind[wind<0] <- NA
    windRange <- range(wind, na.rm=1)
    windMean <- mean(wind, na.rm=1)

    res <- c(G0d, ktd, radWavVar, windRange, windMean)
    })
  res
  }, mc.cores=detectCores())
setwd(old)

fechasOK <- fechas[!sapply(meteo, is.null)]
meteo <- do.call(rbind, meteo)

meteoZ <- zoo(meteo, as.Date(fechasOK))
names(meteoZ) <- c('G0d', 'ktd', paste('wavVar', 1:10, sep=''), 'minWind', 'maxWind', 'meanWind')

## Clustering procedure

foo <- function(x){
  lambda <- powerTransform(x~1)
  res <- bcPower(x, coef(lambda))
}

## Drop Minimum Wind Speed 
trans <- lapply(subset(meteoZ, select=-minWind), foo)
trans <- as.data.frame(trans)

densityplot(as.formula(paste('~', paste(names(trans), collapse='+'), sep='')), data=trans,
            scales=list(y=list(relation='free'), x=list(relation='free')),
            allow.multiple=TRUE,outer=TRUE)

splom(trans,
      panel=panel.hexbinplot,
      colramp=BTC,
      diag.panel = function(x, ...){
        yrng <- current.panel.limits()$ylim
        d <- density(x, na.rm=TRUE)
        d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
        panel.lines(d)
        diag.panel.splom(x, ...)
      },
      pscale=0, varname.cex=0.7)

## Three clusters, PAM method 
nCl=3
## 1=Low, 2=High, 3=Med

## With maxWind and meanWind
pamMeteo0 <- pam(subset(trans, select=-c(G0d)), nCl, stand=TRUE)
silh0 <- silhouette(pamMeteo0)
summary(silh0)$clus.avg.widths

plot(silh0, main='', col=brewer.pal(n=3, name='Dark2'))

## without G0d, maxWind and meanWind
pamMeteo1 <- pam(subset(trans, select=-c(G0d, maxWind, meanWind)), nCl, stand=TRUE)

silh1 <- silhouette(pamMeteo1)
summary(silh1)$clus.avg.widths

plot(silh1, main='', col=brewer.pal(n=3, name='Dark2'))

## We choose cluster without wind
pamMeteo <- pamMeteo1

meteoZ$cluster <- pamMeteo$clustering
trans$cluster <- pamMeteo$clustering

wavVars <- paste('wavVar', 1:10, sep='', collapse='+')
form <- as.formula(paste('ktd', wavVars, sep='~'))
xyplot(form, data=trans, groups=cluster, auto.key=list(space='right'),
       cex=0.5, alpha=0.5, xlab='')

##inspired from https://stat.ethz.ch/pipermail/r-help/2006-April/104068.html
strip.variance <- function(which.given, which.panel, var.name, 
                       factor.levels, ...) {
  expr <- paste("nu[G(0)](lambda[", substr(factor.levels, 7, 8), "])", sep='', collapse=',')
  expr <- paste("expression(", expr, ")", sep = "")
  fl <- eval(parse(text = expr))
  strip.default(which.given, which.panel, var.name, fl, ...)
}

## Figure 5
form <- as.formula(paste('~', wavVars, sep=''))
densityplot(form, data=trans, groups=cluster,
            auto.key=list(corner=c(x=0.6, y=0), text=c('Low', 'High', 'Medium')),
            scales=list(y=list(relation='free'),
              x=list(relation='free')),
            strip=strip.variance,
            xlab='')

## Figure 4
nc <- ncol(trans)
splom(trans[,-nc], groups=trans$cluster,
      auto.key=list(space='right', text=c('Low', 'High', 'Medium')),
      cex=0.2, alpha=0.3, xlab='',
      varname.cex=0.6, pscale=0,
      diag.panel = function(x, varname, ...){
        yrng <- current.panel.limits()$ylim
        d <- density(x, na.rm=TRUE)
        d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
        panel.lines(d, col='lightgray')
        if (grepl('wavVar', varname)) {
          s <- substr(varname, 7, 8)
          varname <- paste("expression(nu(lambda[", s, "]))", sep='')
          varname <- eval(parse(text = varname))
          }
        diag.panel.splom(x, varname,...)
      }
      )

##################################################################
## WAVELET CORRELATION
##################################################################

## The results of this section is available with
## load('data/wavCorDF.RData')
## load('data/corDistLong.RData')

## Days from each cluster
leeZ <- function(d){
  ips <- paste('192.168.26', 16:18, sep='.')
  fechaFormat <- format(d, '%d%m%Y')
  fechaIP <- paste(fechaFormat, ips, sep='_')
  fich <- paste(fechaIP, '.csv', sep='')
  folder <- paste(format(d, '%Y'),
                  paste('CSV', format(d, '%m'), sep='_'),
                  fechaIP,
                  sep='/')
  URLs<- paste(folder, fich, sep='/')
  
  z <- lapply(URLs, function(x) read.zoo(x,
                                         tz='UTC', format='%d/%m/%Y %H:%M:%S',
                                         header=TRUE, sep=','))
  z <- do.call(cbind, z)
  rad <- z[,which(names(z)=="Estaciones.EstPiranometros.Rglobal")]
  rad[rad<0] <- NA
  rad
  }


meteoZ[,c('ktd', 'cluster')]

d1 <- as.Date('2011-11-19')## Low
d2 <- as.Date('2011-06-04')## High
d3 <- as.Date('2011-04-11')## Medium

## Change to the folder where data is available. Sample data is
## available from the authors upon request.
old <- setwd('data/')
rad1 <- leeZ(d1)
rad2 <- leeZ(d2)
rad3 <- leeZ(d3)
setwd(old)

## Maximum Fluctuation
wav1 <- modwt(rad1)
maxFluc1 <- colMaxs(coredata(wav1[,2:10]), na.rm=1)
wav2 <- modwt(rad2)
maxFluc2 <- colMaxs(coredata(wav2[,2:10]), na.rm=1)
wav3 <- modwt(rad3)
maxFluc3 <- colMaxs(coredata(wav3[,2:10]), na.rm=1)

maxFluc <- data.frame(Low=maxFluc1, High=maxFluc2, Med=maxFluc3)/1000*100
maxFluc$scale <- paste('d', 1:9, sep='')

mathLabs <- function(x){
  paste("expression(",
                   paste("lambda[", substr(x, 2, 3), "]", sep='', collapse=',')
                   , ")", sep='')
  }
scaleLabs <- mathLabs(maxFluc$scale)

## Figure 7
dotplot(scale ~ Low + High + Med, data=maxFluc,
        type='o', xlab='Maximum fluctuation (%)',
        scales=list(y=list(labels=eval(parse(text=scaleLabs))))) +
    glayer(panel.text(x[1], y[1], group.value, pos=1, cex=0.7))


index(rad2) <- index(rad3) <- index(rad1) ##todos en el mismo día para verlos bien
rad <- cbind(rad1, rad2, rad3)
names(rad) <- c('Low', 'High', 'Med') ##cluster 1 es baja fluctuación, etc.

xyplot(rad, superpose=TRUE, auto.key=FALSE) +
  glayer(panel.text(x[2200], y[2200], group.value, pos=2, cex=1))


load(paste(d1,'.RData', sep=''))
wavCorLow <- wavCor
corDistLongLow <- corDistLong

load(paste(d2,'.RData', sep=''))
wavCorHigh <- wavCor
corDistLongHigh <- corDistLong

load(paste(d3,'.RData', sep=''))
wavCorMed <- wavCor
corDistLongMed <- corDistLong

toDF <- function(i, lista){
  x <- lista[[i]]
  ord <- order.dendrogram(as.dendrogram(hclust(dist(x))))
  xx <- x[ord, ord]
  data.frame(corr=c(xx), row=c(row(xx)), col=c(col(xx)), scale=paste('d', i, sep=''))
}

wavCorLowDF <- lapply(1:9, toDF, wavCorLow)
wavCorLowDF <- do.call(rbind, wavCorLowDF)
wavCorLowDF$level <- 'Low'

wavCorMedDF <- lapply(1:9, toDF, wavCorMed)
wavCorMedDF <- do.call(rbind, wavCorMedDF)
wavCorMedDF$level <- 'Med'

wavCorHighDF <- lapply(1:9, toDF, wavCorHigh)
wavCorHighDF <- do.call(rbind, wavCorHighDF)
wavCorHighDF$level <- 'High'

wavCorDF <- rbind(wavCorLowDF, wavCorMedDF, wavCorHighDF)
wavCorDF$level <- factor(wavCorDF$level, levels=c('Low', 'Med', 'High'))

## Next figures can be produced with:
## load('data/wavCorDF.RData')
## load('data/corDistLong.RData')

##inspired from https://stat.ethz.ch/pipermail/r-help/2006-April/104068.html
strip.math <- function(which.given, which.panel, var.name, 
                       factor.levels, ...) {
  expr <- paste("lambda[", substr(factor.levels, 2, 3), "]", sep='', collapse=',')
  expr <- paste("expression(", expr, ")", sep = "")
  fl <- eval(parse(text = expr))
  strip.default(which.given, which.panel, var.name, fl, ...)
}


## Figure 8
useOuterStrips(
  levelplot(corr~row*col|scale + level, data=wavCorDF,
            xlab='', ylab='', scales=list(draw=FALSE),
            aspect='xy', par.settings=BTCTheme()),
  strip=strip.math)


## Matrix of correlation-distance
corDistLong <- rbind(
                 cbind(corDistLongLow[,c('distances', 'scale', 'correlation')], level='Low'),
                 cbind(corDistLongMed[,c('distances', 'scale', 'correlation')], level='Med'),
                 cbind(corDistLongHigh[,c('distances', 'scale', 'correlation')], level='High')
                 )

## Figure 9 (without exponential model)
useOuterStrips(
  xyplot(correlation~distances|scale*level, data=corDistLong,
         type=c('p', 'smooth'), auto.key=list(space='right'), ylab='',
         par.settings=rasterTheme(
           col.points=brewer.pal(n=9, name='Blues')[7],
           col.line='black',
           alpha.points=0.5, cex=0.5, lwd=3)),
  strip=strip.math)


################################################################################
## EXPONENTIAL MODEL
################################################################################

### pp. 32 y ss. de Libro "Nonlinear regression with R"
expModel <- function(predictor, a, b, c) a + b * exp(-predictor/c)

expModelInit <- function(mCall,LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  a <- mean(tail(xy$y), na.rm=1) ##asimptota
  xy$y <- xy$y - a
  lmFit <- lm(log(xy[, "y"]) ~ xy[, "x"])
  coefs <- coef(lmFit)
  b <- exp(coefs[1])
  c <- -1/coefs[2]
  value <- c(a, b, c)
  names(value) <- mCall[c("a", "b", "c")]
  value
}

SSexp <- selfStart(expModel, expModelInit, c("a", "b", "c"))

## Nonlinear regression function with multiple models
fooNLS <- function(data){
  control <- nls.control(warnOnly=TRUE)
  modelAsymp <- try(nls(correlation ~ SSasymp(distances, Asym, R0, lrc), data=data, control=control))
  modelExp <- try(nls(correlation~SSexp(distances, a, b, c), data=data, control=control))
  models <- list(modelAsymp, modelExp)
  classModel <- lapply(models, class)
  idx <- which(classModel!='try-error')[1]
  models[[idx]]
}

corList <- split(corDistLong, corDistLong[,c('scale', 'level')])
nlsFull <- lapply(corList, fooNLS)

## Coefficients for each model
coefs <- lapply(nlsFull, function(x){
  if (is.null(x)) {
    coefs <- rep(NA, 4)
    } else {
      coefs <- coef(x)
      model <- all.names(formula(x))[3]
      coefs <- switch(model,
                      SSexp={
                        c(coefs, 1)}, ##pwr=1
                      SSweibull={
                        coefs[3] <- 1/exp(coefs[3]) ##exp(lrc) 
                        coefs[2] <- -coefs[2] ## -Drop
                        coefs},
                      SSasymp={
                        coefs[2] <- coefs[2] - coefs[1] ##R0 - Asymp
                        coefs[3] <- 1/exp(coefs[3]) ##exp(lrc)
                        c(coefs, 1)} ##pwr=1
                      )
            }
  coefs
  })

nms <- data.frame(do.call(rbind, strsplit(names(coefs), '\\.')))
names(nms) <- c('scale', 'level')

coefs <- data.frame(do.call(rbind, coefs))
row.names(coefs) <- NULL
names(coefs) <- c('a', 'b', 'c', 'pwr')
coefs <- cbind(coefs, nms)

coefs$model <- lapply(nlsFull, function(x)all.names(formula(x))[3])
coefs$ABequal <- with(coefs, abs((a+b)/a)<1)
coefs$AB <- with(coefs, a+b)

## inspection
coefs$valid <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
                 FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
                 FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)

scaleLabs <- mathLabs(unique(coefs$scale[coefs$valid]))
coefs$level <- factor(coefs$level, levels=c('Low', 'High', 'Med'))

## Figure 10
trellis.device(pdf, file='range.pdf')
xyplot(log(c)~scale, groups=level, data=coefs, subset=valid,
       ylab='Logarithm of the range factor',
       xlab='Wavelet scale',
       scales=list(x=list(labels=eval(parse(text=scaleLabs)))),
       type=c('p', 'r', 'g'),
       auto.key = list(x = .8, y = .1))
dev.off()


## Figure 9
predFull <- lapply(nlsFull, function(x)if(is.null(x)) rep(NA, 2415) else predict(x))
corDistLong$pred <- do.call(c, predFull)
useOuterStrips(xyplot(correlation+pred~distances|scale*level,
                      data=corDistLong[order(corDistLong$distances),],
                      type=c('p', 'l'), distribute.type=TRUE, ylab='',
                      par.settings=rasterTheme(
                        col.points=brewer.pal(n=9, name='Blues')[7],
                        col.line='black',
                        alpha.points=0.5, cex=0.5, lwd=3)),
               strip=strip.math)
