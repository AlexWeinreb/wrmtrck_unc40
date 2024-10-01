
####### The main function that will be used to process each individual frame. #########

processFrame = function(f, plot=FALSE) {
  
  returnVal = list( frame=f, error=0, area = NA, 
                    dL=NA, midPoint=c(NA,NA), endPoints = matrix(NA,2,2), 
                    meanAngle=NA, bbSplineWeights=rep(NA, ncol(bbBasis)),
                    intensity=rep(NA, nrow(intFit)), 
                    pharIntensity=rep(NA, length(pharynxPointsAlongBackbone)), 
                    widthSplineWeights=matrix(NA, nrow=ncol(widthBasis), ncol=2)
  )
  
  cat("\nFrame " %+% f %+% " |")
  rawFull = readJPEG(f)                     # 3-4 ms/frame
  rawFull = rawFull[nrow(rawFull):1,]       # reorient image so image axes are consistent with stage axes
  raw = EBImage::resize(rawFull, dim(rawFull)[1]/imgDs)
  thresholded = zeroEdges(threshold(1-raw, imgThreshold))   # 4 ms/frame (at full res)
  cat("| loaded and thresholded. ")
  
  dilated = woRmTools::boxDilate(thresholded, round(dilrodeIters/imgDs));
  eroded  = woRmTools::boxErode(dilated, round(dilrodeIters/imgDs));     # 80 ms/frame for both dilation and erosion (10 steps, full res)
  cat("| eroded+dilated")
  
  labeled = 1.0*bwlabel(zeroEdges(eroded));
  
  if (max(labeled)>1) {
    labeled[labeled != which.max(computeFeatures.shape(labeled)[,'s.area']) ] = 0 # set to zero any stuff that doesn't belong to the largest component 
    labeled=labeled/max(labeled)
  }
  cat("| labeled. ")
  
  perimeterImage = woRmTools::thinImage(woRmTools::boxDilate(labeled, 1) - labeled, 5)
  #returnVal$perimeterLength = imgDs * sum(perimeterImage)
  returnVal$area = imgDs^2 * sum(labeled)
  returnVal$medBg = median(raw[!labeled])
  
  thinned = woRmTools::thinImage(labeled, max(dim(labeled)));
  cat("| thinned. ")
  
  if(all(thinned == 0)){
    stop("Problem with frame ", f)
  }
  graphRaw = getGraphFromThinnedImage(thinned) 
  nodes = which(thinned > 0, arr.ind=TRUE) # 3ms
  #returnVal$nodes = nodes
  #returnVal$graph = graphRaw
  
  # check if the graph is a tree. if its not, its got a cycle, so we check all possible trees we could get by breaking cycles.
  if (!isTree(graphRaw)){
    returnVal$error = 1
    if (plot){
      image(x = 1:dim(raw)[1], y=1:dim(raw)[2], z=raw, col=gray(0:100/100), 
            main='frame'%+%f, xlim=c(0,max(dim(raw))), ylim=c(0,max(dim(raw))))
      points(which(perimeterImage > 0, arr.ind=TRUE), pch=16, cex=.3, col='red')
      points(nodes, pch=16, cex=.1, col='green')
    }
    return(returnVal)
  }
  
  
  # get the centerline as the longest shortest path in the tree, smooth and interpolate it.
  centerLine = nodes[getTreeDiameterPath(graphRaw), ]
  if(nrow(centerLine) < centerLineFilterLength/imgDs - 1){
    stop("Problem with frame: ", f)
  }
  
  centerLineFiltered = runmean2(centerLine, k=centerLineFilterLength/imgDs, endrule='keep')
  centerLineResampled = woRmTools::interpArc(centerLineFiltered, nOut=nBbAngles+1)
  colnames(centerLineResampled) = c("y", "x")
  
  angles = signal::unwrap(atan2(diff(centerLineResampled[,1]), diff(centerLineResampled[,2])))
  meanAngle = mean(angles)
  angles = angles - meanAngle
  
  bbSplineWeights = bbFit %*% angles
  dL = sum(sqrt(diff(centerLineResampled[,1])^2 + diff(centerLineResampled[,2])^2)) / (nrow(centerLineResampled)-1)
  anglesPredicted = bbBasis %*% bbSplineWeights + meanAngle
  
  centerLinePredicted = cbind(y = centerLineResampled[1,1] + c(0, dL* cumsum(sin(anglesPredicted))), 
                              x = centerLineResampled[1,2] + c(0, dL* cumsum(cos(anglesPredicted))) )
  
  if (plot){
    image(x = 1:dim(raw)[1], y=1:dim(raw)[2], z=raw, col=gray(0:100/100), main='frame'%+%f, xlim=c(0,max(dim(raw))), ylim=c(0,max(dim(raw))))
    points(which(perimeterImage > 0, arr.ind=TRUE), pch=16,cex=.3, col='red')
    points(nodes, pch=16, cex=.1, col='green')
    points(centerLineResampled, col='blue',type='l')
    points(centerLinePredicted, col='red',type='l')
    points(centerLineResampled[1,1], centerLineResampled[1,2], pch=8, col='red', cex=2)
  }
  
  # Set up return values related to centerline. 
  # note that we rescale some of these so they are returned at the native resolution (and not in downsampled coordinates)
  returnVal$dL                 = dL * imgDs;                                                                 #48 bytes    
  returnVal$midPoint           = as.integer(centerLineResampled[nrow(centerLineResampled)/2, ] * imgDs);     #48 bytes
  returnVal$endPoints          = as.integer(centerLineResampled[c(1,nrow(centerLineResampled)), ] * imgDs);  #56 bytes
  returnVal$meanAngle          = meanAngle;                                                                  #48 bytes
  returnVal$bbSplineWeights    = as.numeric(bbSplineWeights);                                                #168 bytes
  
  # make sure the centerline (compressed representation) remains in the frame
  withinFrame = function(x) x[1] %within% c(0, nrow(thinned)+1) & x[2] %within% c(0, ncol(thinned)+1)
  if(!all(apply(centerLinePredicted, 1, withinFrame))){
    returnVal$error = 2
    return(returnVal)
  }
  
  # here we will get the width of the worm along the centerline
  widths = woRmTools::getWormWidths(labeled, dL, centerLinePredicted[1,], anglesPredicted)
  widthSplineWeights = widthFit %*% widths[1:nBbAngles, ]
  
  # using the full image, get the avg intensity along a short path orthogonal to the backbone at each backbone point
  intensity = sapply(1:nrow(bbBasis), function(i) {
    a = anglesPredicted[i] +  pi/2
    pts = t( centerLinePredicted[i,]*imgDs + outer(c(sin(a), cos(a)), (-4:4)) )
    pts = pts[ pts[,1] >= 1 & pts[,1] <= nrow(rawFull) & pts[,2] >= 1 & pts[,2] <= ncol(rawFull), ] # make sure its within the confines of the image
    mean(rawFull[pts]) 
  })
  
  returnVal$widthSplineWeights = matrix(widthSplineWeights, ncol=2) * imgDs;                                 #712 bytes
  returnVal$intensity          = as.integer(2^16* (intFit %*% intensity))                                    #248 bytes
  returnVal$pharIntensity      = as.integer(2^16* intensity[pharynxPointsAlongBackbone])                     #1040 bytes
  
  return(returnVal);
}

