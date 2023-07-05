
buildTetrachoricPPMC = function(predictiveSamples, rawData, ...) {
  nObs = nrow(rawData)
  nItems = ncol(rawData)
  itemNames = names(rawData)
  
  # build posterior predictive samples matrix for tetrachoric correlations
  ttcPPMC = matrix(data = NA, nrow = nrow(predictiveSamples), ncol = nItems * (nItems - 1) / 2)
  tempNames = NULL
  for (item1 in 1:(nItems - 1)){
    for (item2 in (item1 + 1):nItems){
      tempNames = c(tempNames, paste0("TTC_", itemNames[item1], "_", itemNames[item2]))
    }
  }
  colnames(ttcPPMC) = tempNames
  simData = matrix(data = NA, nrow = nObs, ncol = nItems)
  
  pb = txtProgressBar()
  for (draw in 1:nrow(predictiveSamples)){
    for (item in 1:nItems){
      simData[, item] = predictiveSamples[draw, paste0("simY[", item, ",", 1:nObs, "]")]
    }
    
    col = 1
    for (item1 in 1:(nItems-1)){
      for (item2 in (item1+1):nItems){
        ttcPPMC[draw, col] = blatent:::findTetrachoric(data = simData, var1 = item1, var2 = item2, correct = .5)
        col = col + 1
      }
    }
    setTxtProgressBar(pb, draw/nrow(predictiveSamples))
  }
  
  # build data frame from draw
  col = 1
  obsTTC = matrix(data = NA, nrow = 1, ncol = nItems * (nItems - 1) / 2)
  for (item1 in 1:(nItems-1)){
    for (item2 in (item1+1):nItems){
      obsTTC[1, col] = blatent:::findTetrachoric(data = rawData, var1 = item1, var2 = item2, correct = .5)
      col = col + 1
    }
  }

  
  PPMCstats = NULL
  residTTC = matrix(data = 0, nrow = nItems, ncol = nItems)
  colnames(residTTC) = rownames(residTTC) = itemNames
  
  
  for (column in 1:length(obsTTC)){
    # build empirical distribution
    ttcdist = ecdf(ttcPPMC[, column])
    
    if (ttcdist(obsTTC[column]) > .95 || ttcdist(obsTTC[column]) < .05) {
      ttcFlag = TRUE
    } else {
      ttcFlag = FALSE
    }
    
    item1 = strsplit(x = colnames(ttcPPMC)[column], split = "_")[[1]][2]
    item2 = strsplit(x = colnames(ttcPPMC)[column], split = "_")[[1]][3]
    residTTCval = mean(ttcPPMC[, column]) - obsTTC[column]
    residTTC[item1, item2] = residTTC[item2, item1] = residTTCval
    
    # tabulate item correlation model fit statistics
    PPMCstats = rbind(PPMCstats, 
                      data.frame(
                        item1 = item1,
                        item2 = item2,
                        obsTTC = obsTTC[column],
                        estTTCmean = mean(ttcPPMC[, column]),
                        estTTCsd = sd(ttcPPMC[, column]),
                        estTTC05 = quantile(ttcPPMC[, column], probs = 0.05),
                        estTTC95 = quantile(ttcPPMC[, column], probs = 0.95),
                        residTTC = residTTCval,
                        absResidTTC = abs(mean(ttcPPMC[, column]) - obsTTC[column]),
                        PPPvalue = ttcdist(obsTTC[column]),
                        ttcFlag = ttcFlag
                      )
    )
    rownames(PPMCstats) = NULL
  }
  
  close(pb)
  
  # order by discrepancy 
  PPMCstats = PPMCstats[order(PPMCstats$absResidTTC, decreasing = TRUE), ]
  return(list(PPMCstats = PPMCstats, residTTC = residTTC))
}