## convention throughout is that image coordinates (row, col) = (y,x), which are oriented the same way as stage coordinates
## convention is also that stage coordinates refer to the center of the frame


# Specify file locations
csvFile =       ""
jpegDirectory = "E:/images/06232018/06232018am/images/"
lawnBoundaryFile = NULL                             # set to NULL if not measured

# Imaging/instrument parameters
umPerPx = 1.44               # microns per pixel
umPerStageUnit = 1.25        # microns per stage unit
frameRate = 20               # Hz
boundingBoxPadding = 40      # pixels
fullImageWidth = 1280        # pixels
fullImageHeight = 1024       # pixels


# Image analysis parameters
imgDs = 2                    # downsampling factor, in pixels
#imgThreshold = 0.5           # brightness (0-1)
dilrodeIters = 10            # in pixels in NATIVE resolution (not downsampled)
centerLineFilterLength = 31  # in pixels in NATIVE resolution (not downsampled)

nBbAngles = 1000             # number of backbone angles. Typically 500 or 1000
bbKnots = c(0, 0.05, 1:9/10, 0.95, 1) # locations at which 2nd derivative is allowed to be discontinuous. empirically 13 bbKnots seems sufficiently flexible
widthKnots = 0:30/30         # location of knots for recording body width along centerline. more granular than the angles because we want to detect eggs (which are ~5% of a worm's length)
intensityKnots = seq(0,1, length.out=50)

####################################################
############ Find image files ######################

setwd(jpegDirectory)                                           # set the jpeg directory as current working directory
jpegFiles = list.files(jpegDirectory, pattern=".jpeg")         # get a list of all the jpegs there
jpegFrameIndices = sapply(strsplit(jpegFiles,"_"), function(x) as.numeric(x[1])) # extract the frame index
jpegFiles = jpegFiles[order(jpegFrameIndices)]                 # make sure files are ordered by increasing frame index
jpegFrameIndices = jpegFrameIndices[order(jpegFrameIndices)]   # make sure frame indices are sorted as well 

#Get pixel intensity information from video
pixelIntMat = matrix(nrow = 100, ncol=500)
framInd = seq(from=10, to=length(jpegFiles), length.out=100) #100 frames to sample

for(aa in 1:length(framInd)) {
  
  frameHere = framInd[aa]
  JPEGdata = readJPEG(jpegFiles[frameHere])
  JPEGdata = as.vector(round(JPEGdata*255)) #convert pixel intensities to vector on 0-255 scale
  pixelsToSample = round(seq(from=10, to=length(JPEGdata),length.out=500)) # 500 pixels from throughout frame
  pixelIntMat[aa,] = JPEGdata[pixelsToSample] #write to output matrix
}

ThirtyFifthPercHere = quantile(as.vector(pixelIntMat),0.35)

if(ThirtyFifthPercHere<=254) {
  imgThreshold = 0.5 # brightness (0-1)
} else  {
  imgThreshold = 0.3 # brightness (0-1)
}
######################################################################################
#### load in the stage data, keep only rows corresponding to frames we've recorded ###

stageData = read.csv(csvFile, header=FALSE) 
names(stageData) = c("frameIndex", "xstage","ystage", "zstage", "laserOn","autofocusing", "cx","cy", "area", "bbc1", "bbr1", "bbc2","bbr2","AR", "bodyAngle")
rownames(stageData) = stageData$frameIndex
options(scipen=5)

# truncate either stageData or jpegFiles or both to ensure that we only keep data that is present in both
frameIndices = intersect(jpegFrameIndices, stageData$frameIndex)
jpegFiles = jpegFiles[which(jpegFrameIndices %in% frameIndices)]
jpegFrameIndices = jpegFrameIndices[which(jpegFrameIndices %in% frameIndices)]
stageData = stageData[as.character(jpegFrameIndices), ]


#####################################################
############ Load necessary libraries ###############
loadPackage = function(x) {
  if(!require(x,character.only=TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

if(!require(EBImage)){ # for almost all of the image processing tools 
  if (!require(BiocManager)) install.packages('BiocManager')
  BiocManager::install("EBImage")
}
loadPackage("splines")   # for generating spline bases - function bs() for B-splines
loadPackage("jpeg")      # for reading in jpegs
loadPackage("caTools")   # for runmean()
loadPackage("signal")    # for filter() and butter()
loadPackage("igraph")    # for estimating centerline and finding shortest paths
loadPackage("woRmTools") # written by Nate, contains a handful of image processing functions in C. 
loadPackage("parallel")  # enables parallel image processing, provides significant speedup
loadPackage("Rcpp")
#loadPackage("robust")
loadPackage("zoo")
loadPackage("mgcv")      # for in.out()
#loadPackage("animation") # TODO: IS THIS NEEDED?
loadPackage("MASS")
loadPackage("mmand")     # for dilate/erode of 1D signal

#####################################################
#### Set up basic functions we'll use throughout ####

"%+%" <- function(x,y) paste(x,y,sep='')         # basic string concatenation operator

"%within%" = function(x,y) x>y[1] & x<y[2]       # basic operator to find if x is in range (a,b)

threshold = function(x,t) {0.5+0.5*sign(x-t)}    # threshold a vector or matrix, return binary vector/matrix

# ensure that the edge columns and rows are zero (woRmTools functions do not handle points on edge of image well)
zeroEdges = function(x){                         
  x[,1] = 0; 
  x[1,] = 0; 
  x[nrow(x),] = 0; 
  x[,ncol(x)] = 0; 
  return(x)
}

# check if an index is a local maximum (or minimum) in a vector
is.localmax = function(i, context, k, FUN=max, ...){
    context = c(rep(NA,k), context, rep(NA,k)) # pad with NAs
    i = i + k # account for the NA padding
    contextMatrix = sapply(i, function(x) {context[(x-k):(x+k)]})
    context[i] == apply(contextMatrix, 2, FUN, ...)
}

# check if a graph contains only one connected component
isConnected = function (g) return(max(components(g)$membership) == 1)

isTree = function(g) return(all(get.adjacency(g) == get.adjacency(mst(g))))

# difference in angles (can't exceed pi) 
angleDiff = function(x,y) pmin((x-y) %% (2*pi), (y-x) %% (2*pi))

# takes a matrix with a collection of bright pixels, turns each bright pixel into a node on a graph, connected to its neighbors.
# this function is slow and could probably be improved
getGraphFromThinnedImage = function(thinned){
  nodes = which(thinned>0, arr.ind=TRUE) # 3ms
  d = as.matrix(dist(nodes)) # 25-30ms
  d[lower.tri(d,diag=TRUE)] = 10
  edgeList = which(d<= 1.1, arr.ind=TRUE)
  k = nrow(edgeList)+1
  possibleEdges = which(d > 1.1 & d < 1.5, arr.ind=TRUE)
  edgeList = rbind(edgeList, matrix(0,nrow(possibleEdges),2))
  for (i in 1:nrow(possibleEdges)){
    pt1 = nodes[possibleEdges[i,1], ]
    pt2 = nodes[possibleEdges[i,2], ]
    if (!(thinned[pt1[1], pt2[2]]) & !(thinned[pt2[1], pt1[2]])){
      edgeList[k,] = possibleEdges[i,]
      k=k+1
    }
  }
  edgeList=edgeList[1:(k-1),]
  g = graph_from_edgelist(edgeList, directed=FALSE)
}

# returns node indices (in order) that constitute the longest shortest path in the tree.
# the centerline is going to be the path across the tree at its diameter
getTreeDiameterPath = function(g){
  end1 = which.max( dfs(g, root=1, dist=TRUE)$dist    )
  end2 = which.max( dfs(g, root=end1, dist=TRUE)$dist )
  return(shortest_paths(g, end1, end2)$vpath[[1]])
}

# this is mostly a wrapper for runmean, but fixes the fact that the runmean end rule has a varying group delay.
runmean2 = function(x, k, ...){
  x = as.matrix(x)
  xFilt = runmean(x, k=k, ...)
  for (j in 1:ncol(x)){
    for (i in 1:ceiling(k/2)){ 
      xFilt[i,         j] = mean(x[1:(2*i-1), j])
      xFilt[nrow(x)-i+1, j] = mean(x[nrow(x)-(0:(2*i-2)), j])
    }
  }
  return(xFilt)
}

reversePosture = function(wormPosture){
  wormPosture$bbSplineWeights = rev(wormPosture$bbSplineWeights)
  wormPosture$widthSplineWeights = wormPosture$widthSplineWeights[ncol(widthBasis):1, 2:1]
  wormPosture$intensity = rev(wormPosture$intensity)
  wormPosture$pharIntensity = rev(wormPosture$pharIntensity)
  wormPosture$endPoints = wormPosture$endPoints[c(2,1,4,3)]
  wormPosture$meanAngle = wormPosture$meanAngle + pi
  return(wormPosture)
}

paddedDiff = function(x, lag) {
  c( rep(NA, ceiling(lag/2)), diff(x, lag=lag), rep(NA, floor(lag/2)))
}

minDist = function(pt, ptSet){
  pt = as.numeric(pt)
  ptSet=as.matrix(ptSet)
  distSquared = (pt[1] - ptSet[,1])^2 + (pt[2] - ptSet[,2])^2
  return(sqrt(min(distSquared)))
}

eventsToTimeseries = function(t,n,val=1){
  a = rep(0,n);
  a[t] = val;
  return(a)
}

centroidFromPolygon = function(x,y){
  # x[1],y[1] must equal x[n], y[n] !!
  n = length(x)
  A = 1/2 * sum(x[-n]*y[-1] - x[-1]*y[-n])
  cx = 1/(6*A) * sum( (x[-n]+x[-1]) * (x[-n]*y[-1] - x[-1]*y[-n]) )
  cy = 1/(6*A) * sum( (y[-n]+y[-1]) * (x[-n]*y[-1] - x[-1]*y[-n]) )
  c(x=cx, y=cy)
}

####################################################################################
################ Precompute matrices for backbone and width fitting ################
bbBasis = bs(seq(0,1,length=nBbAngles), knots = bbKnots, degree=2); 
bbBasis = bbBasis[,-ncol(bbBasis)] # backbone basis
bbFit   = solve(t(bbBasis) %*% bbBasis) %*% t(bbBasis)

widthBasis = bs(seq(0,1,length=nBbAngles), knots = widthKnots, degree=2); 
widthBasis = widthBasis[,-ncol(widthBasis)]
widthFit = solve(t(widthBasis) %*% widthBasis) %*% t(widthBasis)

intBasis = bs(seq(0,1,length=nBbAngles), knots = intensityKnots, degree=2); 
intBasis = intBasis[,-ncol(intBasis)]
intFit = solve(t(intBasis) %*% intBasis) %*% t(intBasis)

pharynxPointsAlongBackbone = c(seq(1, 0.25*nBbAngles), seq(0.75*nBbAngles, nBbAngles))




######################################################################################
######################### Load in the lawn boundary ##################################

if (!is.null(lawnBoundaryFile)){
  lawnBoundaryData = read.csv(lawnBoundaryFile,header=FALSE)
  lawnBoundaryData = unique(lawnBoundaryData[,2:3]) *umPerStageUnit/umPerPx # convert to pixels, drop unnecessary columns and redundant rows
  lawnBoundaryData = as.matrix(lawnBoundaryData)                            # make sure its in matrix form for subsequent functions
  colnames(lawnBoundaryData) = c('x','y')                                   # name the columns appropriately
  plot(lawnBoundaryData, pch='.')                                           # plot it
  
  lawnBoundaryData = runmean(as.matrix(lawnBoundaryData), k=2*frameRate)    # smooth it
  points(lawnBoundaryData, pch='.', col='red')                              # plot the 2s-smoothed ata

  lawnBoundaryData = woRmTools::interpArc(rbind(lawnBoundaryData, lawnBoundaryData[1,]), nOut=1000) # downsample the boundary
  points(lawnBoundaryData, pch=16, col='blue', cex=0.2)                           # plot the downsampled data
  
  lawnCenter = centroidFromPolygon(lawnBoundaryData[,1], lawnBoundaryData[,2]) # find the center
}


#######################################################################################
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




#######################################################################################
################### Run processFrame() on each and every frame ########################

# when parallelized on 8-core machine, takes 20 ms/frame at half-resolution (imgDs=2)

#### Unparallelized version, for testing ####
#lapply(jpegFiles[round(seq(1,length(jpegFiles), length.out=20))], processFrame, plot=TRUE)
#for (i in 1:length(jpegFiles)) processFrame(jpegFiles[i]) # for debugging only


nFramesToProcess = length(jpegFiles) 
#framesToProcess = round(seq(1,length(jpegFiles),length.out=nFramesToProcess))     # use this to analyze evenly-spaced frames
framesToProcess = 1:nFramesToProcess                                               # use this to analyze the first n frames   


no_cores <- detectCores() - 1                                                      # Calculate the number of cores
cl <- makeCluster(no_cores)                                                        # Initiate cluster
clusterExport(cl, setdiff(ls(), c("jpegFiles","stageData", "wormPostures")))       # load data into cluster for parallel image analysis
a = clusterEvalQ(cl, library(signal))                                              # load all relevant libraries on cluster
a = clusterEvalQ(cl, library(caTools))
a = clusterEvalQ(cl, library(jpeg))
a = clusterEvalQ(cl, library(splines))
a = clusterEvalQ(cl, library(woRmTools))
a = clusterEvalQ(cl, library(EBImage))
a = clusterEvalQ(cl, library(igraph))

t=proc.time()                                                                      # time how long analysis takes
wormPostures = parLapply(cl, jpegFiles[framesToProcess], processFrame, plot=FALSE) # run analysis!
print(proc.time()-t)
stopCluster(cl)

saveRDS(wormPostures, jpegDirectory %+% "../wormPosturesRaw.rds")





















# First, reorient worms so that they are consistently oriented frame-to-frame (typically 170ms for 1000 frames) 
for (k in 2:length(wormPostures)){
  prevAngles = wormPostures[[k-1]]$meanAngle + wormPostures[[k-1]]$bbSplineWeights
  newAngles =  wormPostures[[k]]$meanAngle + wormPostures[[k]]$bbSplineWeights
  d1 = sum(angleDiff(prevAngles, newAngles))
  d2 = sum(angleDiff(prevAngles, pi + rev(newAngles)))
  if (is.na(d1) | is.na(d2)) next 
  if (d2 < d1){
    wormPostures[[k]] = reversePosture(wormPostures[[k]])
  }    
  if (k %% 1e5 == 0)
    cat(".")
}

# Second, make sure each "segment" of frames between self-intersecting frames is consistent 
# with all the segments before it. For each contiguous series of frames, flip the orientation if it
# makes the intensity values more consistent with the previous ones. 
errors = sapply(wormPostures, function(x) x$error)
contiguousSegmentStarts = which(c(errors,1) == 0 & c(1,errors) != 0)   # current is 0, prev (padded w 1) is error
contiguousSegmentEnds =   which(c(errors,1) != 0 & c(1,errors) == 0)-1 # current (padded) is error, prev is 0



prevAvgIntensity = rep(0, nrow(intFit))
n=0
for (j in 1:length(contiguousSegmentStarts)){
  a = contiguousSegmentStarts[j]
  b = contiguousSegmentEnds[j]
  i = seq(a, b, length.out = min(b-a+1, 1000))
  intensityMatrix = t(sapply(wormPostures[i], function(x) if(is.null(x$intensity)) rep(NA,nBbAngles) else x$intensity))
  avgIntensity = colMeans(intensityMatrix, na.rm=TRUE)
  
  if (j == 1) { # first segment - take our best to guess which is head/tail
    headIntensity = mean(avgIntensity[seq(1,length(avgIntensity)*0.1)])
    tailIntensity = mean(avgIntensity[seq(length(avgIntensity)*0.9, length(avgIntensity))])
    if (tailIntensity > headIntensity){    # if (the tail looks brighter than the head)
      avgIntensity = rev(avgIntensity)     #   flip average intensity vector by 180 degrees (reverse it)
      for (k in a:b)                       #   flip every frame within this segment so that the brighter part is the head
        wormPostures[[k]] = reversePosture(wormPostures[[k]])
      
    }
  } 
  else {  # later segment - match intensity profile to previous contiguousSegment (to ensure worm never flips 180)
    if ( sum((avgIntensity-prevAvgIntensity)^2) > sum((rev(avgIntensity)-prevAvgIntensity)^2)){
      avgIntensity = rev(avgIntensity)     #   flip average intensity vector by 180 degrees (reverse it)
      for (k in a:b)                       #   flip every frame within this segment so that the brighter part is the head
        wormPostures[[k]] = reversePosture(wormPostures[[k]])
    }
  }
  
  prevAvgIntensity = (b-a+1)/(n+b-a+1)*avgIntensity + n/(n+b-a+1)*prevAvgIntensity
  n = n + (b-a+1) # keep track of total number of frames averaged
}


# create indices for easy downsampling
i.3 = framesToProcess[seq(1, length(framesToProcess), length.out=1e3)]
i.4 = framesToProcess[seq(1, length(framesToProcess), length.out=1e4)]
i.5 = framesToProcess[seq(1, length(framesToProcess), length.out=1e5)]
i.6 = framesToProcess[seq(1, length(framesToProcess), length.out=1e6)]

i.3c = framesToProcess[(length(framesToProcess)/2)+ (-5e2:5e2)]
i.4c = framesToProcess[(length(framesToProcess)/2)+ (-5e3:5e3)]
#i.5c = framesToProcess[(length(framesToProcess)/2)+ (-5e4:5e4)]
#i.6c = framesToProcess[(length(framesToProcess)/2)+ (-5e5:5e5)]

# Verify that all frames have been appropriately oriented by making an intensity heatmap! 
# Then, have the user inspect  the results and see if we need to flip the overall result.
intensityMatrix = t(sapply(wormPostures[i.4], function(x) if(is.null(x$intensity)) rep(NA,ncol(intBasis)) else x$intensity))
image(x=i.4/frameRate, z=intensityMatrix, xlab='time (s)', ylab='body index (nose to tail)', col=gray(0:100/100))

while (readline("\nIs the worm correctly oriented? (y/n): ") != 'y'){
  for (k in 1:length(wormPostures))  
    wormPostures[[k]] = reversePosture(wormPostures[[k]])
  intensityMatrix = t(sapply(wormPostures[i.4], function(x) if(is.null(x$intensity)) rep(NA,ncol(intBasis)) else x$intensity))
  image(x=i.4/frameRate, z=intensityMatrix, xlab='time (s)', ylab='body index (nose to tail)', col=gray(0:100/100))
}



# remove extraneous data about pharynx from tail end (we only saved it because we werent sure which end was the head)
for (j in 1:length(wormPostures)){
  wormPostures[[j]]$pharIntensity = wormPostures[[j]]$pharIntensity[1:(length(pharynxPointsAlongBackbone)/2)]
}


save.image(file=jpegDirectory %+% "../analyzedFrames.RData")


























#####################################################################################
########################### Calculate locomotion metrics ############################

# calculate worm's x and y position (in pixels). Here we'll take the image center to be (0,0)
# so that it matches up with the lawnboundary data (where we try to keep the lawn boundary in the center).
stageData$x = stageData$xstage*umPerStageUnit/umPerPx + stageData$cx-fullImageWidth/2   # pixels. 
stageData$y = stageData$ystage*umPerStageUnit/umPerPx + stageData$cy-fullImageHeight/2  # pixels. 

# translate so that origin is the center of the lawn
if (exists('lawnCenter')){
  stageData$x = stageData$x - lawnCenter['x']
  stageData$y = stageData$y - lawnCenter['y']
}
  

# calculate displacements on various timescales
dx1 = paddedDiff(stageData$x, frameRate)
dy1 = paddedDiff(stageData$y, frameRate)
dx5 = paddedDiff(stageData$x, 5*frameRate)
dy5 = paddedDiff(stageData$y, 5*frameRate)
dx10 = paddedDiff(stageData$x, 10*frameRate)
dy10 = paddedDiff(stageData$y, 10*frameRate)

# calculate speed based on absolute displacement of center of mass
comSpeed1 =  sqrt( dx1^2  + dy1^2  )
comSpeed10 = sqrt( dx10^2 + dy10^2 ) / 10

# calculate speed projected onto body axis. This could probably be refined. 
meanAngle = sapply(wormPostures, function(x) x$meanAngle) + pi # add pi because angle is originally head->tail and we want tail->head
bodyAxisSpeed1 =  (dx1 * cos(meanAngle)  + dy1 *  sin(meanAngle))
bodyAxisSpeed10 = (dx10 * cos(meanAngle) + dy10 * sin(meanAngle)) / 10


# calculate movement direction and angular velocity on a couple different timescales
movementDirection1 =  atan2(dy1,  dx1)
movementDirection10 = atan2(dy10, dx10)

angularVelocity1 = paddedDiff(movementDirection1, 1*frameRate)
angularVelocity1 = angularVelocity1 %% (2*pi)
angularVelocity1[which(angularVelocity1 > pi)] = angularVelocity1[which(angularVelocity1 > pi)] - 2*pi 

angularVelocity10 = paddedDiff(movementDirection10, 10*frameRate) 
angularVelocity10 = angularVelocity10 %% (2*pi)
angularVelocity10[which(angularVelocity10 > pi)] = angularVelocity10[which(angularVelocity10 > pi)] - 2*pi 
angularVelocity10 = angularVelocity10/10


# figure out the worm's head position
xHead = sapply(1:length(wormPostures), function(i) wormPostures[[i]]$endPoints[3] +  # px relative to jpeg boundary
                 stageData$bbc1[i] - # px of jpeg boundary relative to full image
                 fullImageWidth/2 +  # px of full image boundary relative to stage position
                 stageData$xstage[i]*umPerStageUnit/umPerPx  # stage position (px)
)                
yHead = sapply(1:length(wormPostures), function(i) wormPostures[[i]]$endPoints[1] +  # px relative to jpeg boundary
                 (fullImageHeight-stageData$bbr2[i]) - # px of jpeg boundary relative to full image
                 fullImageHeight/2 +  # px of full image boundary relative to stage position
                 stageData$ystage[i]*umPerStageUnit/umPerPx  # stage position (px)
)

# figure out when the worm was on food and how far from the lawn boundary it was
onFoodHead = rep(NA,length(xHead))
onFoodHead[!is.na(xHead)] = in.out(as.matrix(lawnBoundaryData), cbind(xHead,yHead)[!is.na(xHead),])
distToFoodHead = sapply(1:nrow(stageData), function(i) minDist(pt=c(xHead[i], yHead[i]), ptSet=lawnBoundaryData)) #10 min for a 6h-video
distToFoodHead = distToFoodHead*sign(0.5-onFoodHead)

onFoodBody = in.out(as.matrix(lawnBoundaryData), as.matrix(stageData[,c('x','y')]))
distToFoodBody = sapply(1:nrow(stageData), function(i) minDist(pt=stageData[i,c('x','y')], ptSet=lawnBoundaryData)) #10 min for a 6h-video
distToFoodBody = distToFoodBody*sign(0.5-onFoodBody)


# make plots of speed,angular speed, and position relative to the lawn
i = i.4 # downsample data for plotting purposes
layout(matrix(c(1,2,3,4,4,5),ncol=2), widths=c(2,1))
par(mar=c(4,2,1,1))
# plot absolute speed (displacement of center of mass)
plot(i/frameRate/60, comSpeed1[i]*umPerPx/1000, type='l', col=adjustcolor('black',0.5), xlab='time (minutes)')
lines(i/frameRate/60, comSpeed10[i]*umPerPx/1000,col='red')
points(i[!onFoodBody[i]]/frameRate/60, rep(max(comSpeed1[i]*umPerPx/1000, na.rm=TRUE), sum(!onFoodBody[i], na.rm=TRUE)), col='darkred', pch="|")
abline(h=0)
# plot speed along body axis
plot(i/frameRate/60, bodyAxisSpeed1[i]*umPerPx/1000, type='l', xaxs='i',col=adjustcolor('black',0.5), xlab='time (minutes)')
lines(i/frameRate/60, bodyAxisSpeed10[i]*umPerPx/1000,col='red')
abline(h=0)
# plot angular velocity
plot(i/frameRate/60, angularVelocity1[i], type='l', xaxs='i',col=adjustcolor('black',0.5), xlab='time (minutes)')
lines(i/frameRate/60, angularVelocity10[i],col='red')
abline(h=0)
# plot x-y coordinates of worm over time, color-coded by bodyAxisSpeed
cols = adjustcolor(colorRampPalette(c('cyan', 'blue','black','red','orange'))(100), 0.5)
speedIndices = 1+99*(bodyAxisSpeed1-min(bodyAxisSpeed1,na.rm=TRUE))/diff(range(bodyAxisSpeed1,na.rm=TRUE))
plot(stageData$x[i]*umPerPx/1000, stageData$y[i]*umPerPx/1000, type='o', pch=16, cex=.3, xlab='distance (mm)', ylab='distance (mm)')#,  col=cols[speedIndices[i]])

# overlay the lawn boundary onto the trajectory and show when the worm was in or outside the food
lines(lawnBoundaryData*umPerPx/1000, col='blue')
points(stageData$x[i][!onFoodHead[i]]*umPerPx/1000, stageData$y[i][!onFoodHead[i]]*umPerPx/1000, col='darkred', pch=16,cex=0.2)

# plot distance to food
plot(i/frameRate/60, distToFoodHead[i], type='l', ylab='dist to food', xlab='time (minutes)')
abline(h=0)


#################################################################################
#################### Detect defecation motor program ############################

L = sapply(wormPostures, function(x) x$dL) * nBbAngles * umPerPx # worm length, in microns
L = zoo::na.approx(L)                                            # interpolate to avoid missing data
L.lp = runmed(L, k=1+50*frameRate)                                 # calculate lowpass with moving median filter
L.hp = (L - L.lp)/L.lp                                                 # hp is the high-pass filtered estimate of fractional change in body length

# normalize DMP-finding kernel such that input with peak amplitude 1 yields output 1.
kernel = -dnorm(seq(0, 40, by=1/frameRate), 20, 0.9) - dnorm(seq(0, 40, by=1/frameRate), 23.5, 0.9)
kernel = kernel/diff(range(kernel))
kernel = kernel/sum(kernel^2)

bp = stats::filter(L.hp, kernel)                                           # filter with DMP-finding kernel
aboveThreshold = which(bp > 0.03)                                        # find where body is contacted by roughly greater than 3% of the worm's total length
DMPs = aboveThreshold[is.localmax(aboveThreshold, bp, k=1+10*frameRate)] # find the local maxima of the contractions
DMPs = na.omit(DMPs)

i = i.4c # for downsampling 
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,5,6),nrow=4, byrow=TRUE))
par(mar=c(4,4,1,1))
plot(i/frameRate/60, L[i], type='l',ylab='body length (um)', bty='n', yaxt='n', xlab='time (min)')
axis(2,las=2)
abline(v=DMPs/frameRate/60, col=adjustcolor('red',0.2), lwd=2)

plot(bp[i], type='l')
abline(h=0.025)

periods = diff(DMPs)/frameRate # in seconds
plot(DMPs[-1]/frameRate/60, periods, xlab='minutes', ylab='defecation cycle period (seconds)', ylim=c(0,100),pch=16,cex=0.7)
lines(lowess(periods ~ I(DMPs[-1]/frameRate/60), f=.3), col='red',lwd=2)

# histogram of time between defecation events
hist(periods,breaks=500, xlim=c(0,100), main='Defecation cycle periods', xlab='time (sec)', ylab='count')

# let's look at the DMP average signal
contractionMatrix = sapply(DMPs, function(a) {L.hp[a + -(length(kernel)/2):(length(kernel)/2) ]})
plot(1:nrow(contractionMatrix)/frameRate, rowMeans(contractionMatrix), xlab='time (s)', ylab='change in length (%)', type='l')

# make a heatmap of the length during each DMP event
image(contractionMatrix)


###############################################################################
##################### identify pumping ########################################

# interpolate to single pixel sampling
n = length(pharynxPointsAlongBackbone)/2
intensity = t(sapply(wormPostures, function(x) {
  x_interp = seq(1,160,by=1);
  if(all(is.na(x$pharIntensity))) rep(NA,length(x_interp)) 
  else approx(x=1:n*x$dL, x$pharIntensity, xout=x_interp)$y
}))

# Remove potentially erroneous frames based on worm length (especially head-to-egg contact).
# Estimate envelope as 90th percentile over 10 minute windows.
envelope = apply(matrix(L, nrow=12000), 2, quantile, probs=c(0.9), na.rm=TRUE)
envelope = approx(x=(1:length(envelope)-1)*12000+6000, y=envelope, xout=1:length(L), rule=2)$y # Upsample
exceedsEnvelope = which(L - envelope > 30) # identify frames that exceed the length envelope by 30um. 
intensity[exceedsEnvelope,] = NA           # remove the frames (by setting their intensities to NA)

# plot the pharynx intensity
par(mfrow=c(3,1),mar=c(4,4,2,2))
i = seq(1,10000, length.out=2000)
image(x=i/frameRate, y=1:ncol(intensity), z=intensity[i,], main='(pharynx) intensity along centerline', xlab='time (seconds)', col=gray(0:100/100))

intensity.hp = t(t(intensity) - runmean2(t(intensity), k=15))
image(x=i/frameRate, y=1:ncol(intensity.hp), z=intensity.hp[i,], main='di', col=gray(0:100/100))

getPharynxPoints = function(x, k, n=4, us=1, distRange=c(8,17)){
  if (all(is.na(x))) return(c(NA,NA))
  
  if (us>1){
    x = spline(x,n=us*length(x))$y
    k=k*us
  }
  minima = which(is.localmax(i = k:(length(x)-k), context=x, k=k, FUN=min)) + (k-1)
  o = order(x[minima])
  m = minima[o][1:min(n, length(minima))]
  i=1
  while(!any(abs(m[i]-m) %within% (distRange*us)) & i < length(m))
    i=i+1
  if (i==length(m))
    return(c(NA,NA))
  
  return(sort(m[c(i, min(which(abs(m[i]-m) %within% (distRange*us))))])/us)
}


pharynxPoints = t(apply(intensity.hp, 1, function(x,k) getPharynxPoints(x, k=5, us=5) )) # slowwww. sometimes >10-20 minutes for a 6-hour dataset
pharynxDist = as.numeric(diff(t(pharynxPoints)))
points(i/frameRate, pharynxPoints[i,1], col=adjustcolor('red',1),pch='-',cex=1)
points(i/frameRate, pharynxPoints[i,2], col=adjustcolor('blue',1),pch='-',cex=1)

# for every contiguous segment greater than 1s
#   filter it, detect pumps
getContiguousSegments = function(x){ # start/end are inclusive
  bad = is.na(x)
  starts = which(c(TRUE, bad) & c(!bad,FALSE))
  ends = which(c(FALSE, !bad) & c(bad,TRUE))-1
  cbind(starts,ends)
}

segs = getContiguousSegments(pharynxDist)
butterFilt = butter(n=1, type='pass',W=c(0.2,0.8))        # for 20hz acquisition

pumpingEvents = NULL
pumpingValidFrame = rep(FALSE, length(pharynxDist))

par(mar=c(1,1,1,1),mfrow=c(4,4))
for (r in 1:nrow(segs)){
  i=segs[r,1]:segs[r,2]

  if (length(i) < 20) next
  filteredPharynx = as.numeric(signal::filter(butterFilt, pharynxDist[i]-median(pharynxDist[i])))
  belowThresholdIndices = which(filteredPharynx < -0.6)
  if (length(belowThresholdIndices)>0){
    belowThresholdLocalMinima = is.localmax(belowThresholdIndices, filteredPharynx, k=2, FUN=min,na.rm=TRUE)
    pe = belowThresholdIndices[belowThresholdLocalMinima] # for 20 Hz acquisition
    pumpingEvents=c(pumpingEvents, na.omit(pe)+i[1]-1)
  }
  pumpingValidFrame[i]=TRUE
  
  #plot(1:length(filteredPharynx)/frameRate,filteredPharynx*umPerPx,type='o',
  #     pch=16,cex=0.6,xlab='time (s)',ylab='Grinder displacement (um)'); #title(main=r,line=0)
  #abline(h=0)
  #if (length(belowThresholdIndices)>0){
  #  points(pe/frameRate, filteredPharynx[pe]*umPerPx, col='blue',pch=16)
    #abline(v=pe,col='red')
  #}
}

pumpOccurred = rep(0,length(pharynxDist))
pumpOccurred[pumpingEvents] = 1
rawPumpRate    = as.numeric(signal::filter(filt = c(rep(1,frameRate),rep(0,frameRate)), a=1, x=pumpOccurred))
validFrameRate = as.numeric(signal::filter(filt = c(rep(1,frameRate),rep(0,frameRate)), a=1, x=pumpingValidFrame))
pumpingRate = rawPumpRate/validFrameRate*frameRate
pumpingRate[validFrameRate < frameRate] = NA

plot(pumpingRate, type='l', xlab='time (min)', ylab='pumping rate (Hz)')


############################################################################
############################### Detect nose lifts ##########################
#L = sapply(wormPostures, function(x) x$dL) * nBbAngles * umPerPx # worm length, in microns
#L.lp = runmed(L, k=1+5*frameRate)                               # calculate lowpass with moving median filter
#L.hp = (L - L.lp)/L.lp                                           # hp is the high-pass filtered estimate of fractional change in body length
#L.bp = runmean(L.hp, k=frameRate/2)

#noseTipIntensity = sapply(wormPostures, function(x) mean(x$pharIntensity[10:40]) )
#noseTipIntensity.lp = runmed(noseTipIntensity,k=1+10*frameRate)
#noseTipIntensity.hp = (noseTipIntensity-noseTipIntensity.lp)/noseTipIntensity.lp
#noseTipIntensity.bp = runmean(noseTipIntensity.hp, k=frameRate/2)

#noseWidths = sapply(wormPostures, function(x) mean(x$widthSplineWeights[c(1:3,1:3+nrow(widthFit))]) )
#noseWidths.hp = noseWidths - runmed(noseWidths, k=1+5*frameRate)

#plot(L.bp, noseTipIntensity.bp, pch='.'); abline(v=0); abline(h=0);
#abline(v=-0.02, col='red')
#abline(h=-0.35, col='red')
#noseLiftFrames = which(noseTipIntensity.bp < -0.35 & L.bp < -0.02 & stageData$autofocusing==0 & errors==0)

#z = rep(NA, length(wormPostures)); z[noseLiftFrames]=1
#noseLiftPeriods = getContiguousSegments(z)



####################################################################################
########################### Calculate curvature ####################################
diffBbBasis = diff(bbBasis)
curvature = t(sapply(wormPostures, function(x) {
  a = diffBbBasis %*% x$bbSplineWeights
  curvature = tan(a/2)/(x$dL*umPerPx/1000)   # in units of 1/mm
  if (is.null(a))
    return(rep(NA, nBbAngles-1))
  return(curvature[1+bbKnots*(nBbAngles-2)])
}))

i=i.4c
par(mfrow=c(3,1))
cols = colorRampPalette(c('blue','blue','white','black','black'))(100)
# first, show speed
plot(i/frameRate/60, comSpeed1[i+framesToProcess[1]-1]*umPerPx/1000,type='l',xaxs='i',col=adjustcolor('black',0.5), xlab='time (minutes)',ylab='comSpeed (mm/s)')
lines(i/frameRate/60, comSpeed10[i+framesToProcess[1]-1]*umPerPx/1000,col='red')
# now add kymograph
image(x=(i-i[1])/frameRate, z=curvature[i,], col=cols, zlim=c(-10,10), yaxt='n', xlab='time (s)')

par(mar=c(4,2,1,1))
plot(i/frameRate/60, bodyAxisSpeed1[i]*umPerPx/1000,type='l',xaxs='i',col=adjustcolor('black',0.5), xlab='time (minutes)')
lines(i/frameRate/60, bodyAxisSpeed10[i]*umPerPx/1000,col='red')



#####################################################################################
###################  Look at where omega turns happened #############################
errors = sapply(wormPostures, function(x) x$error!=0)
errors = mmand::erode(mmand::dilate(errors, rep(1,9)), rep(1,9))
#omegaStateStarts = c(which(errors & c(TRUE, !errors))-1, nFramesToProcess)
#omegaStateEnds = c(1, which(c(FALSE,errors)  & !errors))



####################################################################
################### Save progress to here ##########################
save(list=setdiff(ls(),c("intensity","intensity.hp")), file=jpegDirectory %+% "../analyzedFramesComplete.RData")



############################################################################
#################### Detect egg laying (semi-automated) ####################

vulvalPositions = c(15:20, 32 + 15:20)
i = seq(1,length(wormPostures))
W = t(sapply(wormPostures[i], function(x) {
  if (is.null(x$widthSplineWeights)) as.integer(rep(NA,2*ncol(widthBasis)))
  else x$widthSplineWeights
  } ))
dW = diff(W, lag=frameRate/10)
area = sapply(wormPostures[i], function(x) x$area)
dArea = diff(area, lag=frameRate/10)
len = sapply(wormPostures[i], function(x) x$dL)*nBbAngles*umPerPx
dLen = diff(len, lag=frameRate/10)

prevSectionWidth = rowSums(W[,vulvalPositions-6])
vulvalSectionWidth = rowSums(W[,vulvalPositions])
postSectionWidth = rowSums(W[,vulvalPositions+6])
eggLayingness = c(0,rowSums(diff(runmean(cbind(0.5*prevSectionWidth, vulvalSectionWidth, 0.5*postSectionWidth), k=5),lag=5)))

## first, identify possible egg-laying events based on a crude heuristic
possibleEggEvents = which(dArea > -100 &
                          eggLayingness > 10 &
                          apply(dW[,vulvalPositions] %within% c(6,40), 1, any) &
                          apply(dW[,vulvalPositions], 1, function(x) !any(is.na(x))) &
                          (dLen/len) %within% c(-0.05,0.05)) + floor(frameRate/10/2)
positionEggDetected = apply(dW[possibleEggEvents-floor(frameRate/10/2), vulvalPositions], 1, which.max)

plot(possibleEggEvents/10/60, positionEggDetected+rnorm(length(positionEggDetected),0,0.15), xlab='time (min)', ylab='egg position')
abline(h=6.5)
p = sum(positionEggDetected < 6.5)/length(positionEggDetected)
title(main="side 1 has proportion " %+% round(p,2))

#Steph determines V/D side
par(mfrow=c(1,2), mar=c(4,4,2,2))
i=1
for (k in possibleEggEvents){
  processFrame(jpegFiles[framesToProcess[k]],TRUE)
  title(line=0, main = framesToProcess[k])
  processFrame(jpegFiles[framesToProcess[k]+1],TRUE)
  title(line=0, main = framesToProcess[k+1])
  title(main=i,line=-2)
  
  result=readline()
  
  if(result=='u'){
    i=i+1
  }
  if(result=='v'){
    if(positionEggDetected[i]<6.5){
      possibleEggEvents = possibleEggEvents[positionEggDetected < 6.5]
    }
    else{
      possibleEggEvents = possibleEggEvents[positionEggDetected > 6.5]
    }
    break
  }
  if(result=='q'){
    break
  }
  
}

#if (p > 0.5) possibleEggEvents = possibleEggEvents[positionEggDetected < 6.5]
#if (p < 0.5) possibleEggEvents = possibleEggEvents[positionEggDetected > 6.5]

eggEventAmplitude = apply(dW[possibleEggEvents,vulvalPositions],1,max)
positionEggDetected = apply(dW[possibleEggEvents,vulvalPositions],1,which.max)


## now begin manual scoring of possible egg events
numEggs = vector('character', length(possibleEggEvents))
par(mfrow=c(1,2), mar=c(4,4,2,2))
i=1
for (j in possibleEggEvents){
  processFrame(jpegFiles[framesToProcess[j]],TRUE)
  title(line=0, main = framesToProcess[j])
  processFrame(jpegFiles[framesToProcess[j]+1],TRUE)
  title(line=0, main = framesToProcess[j+1])
  title(main=i,line=-2)

  numEggs[i] = readline()
  if (numEggs[i]=='q') break

  i=i+1
}

numEggs=as.numeric(numEggs)
eggEvents = possibleEggEvents[which(numEggs >= 1)]
positionEggDetected = apply(dW[eggEvents-floor(frameRate/10/2), vulvalPositions], 1, which.max)
numEggs = numEggs[which(numEggs >= 1)]


# plot density estimates
d = t(sapply(1:length(eggEvents), function(i) matrix(c(eggEvents[i]+c(-1,0,1), c(0,numEggs[i],0)), ncol=2)))
x = as.numeric(t(d[,1:3]))
y = as.numeric(t(d[,4:6]))
plot(x, y, type='l')





########################################################################################################
########## save the analyzed behavior in a spreadsheet (also save the global environment) ##############

allBehaviorMatrix = cbind(
  area = round(sapply(wormPostures, function(x) x$area)*umPerPx^2),           # area in um^2
  length = round(sapply(wormPostures, function(x) x$dL)*nBbAngles*umPerPx),   # length in um
  bodyAngles = round(t(sapply(wormPostures, function(x) x$bbSplineWeights+x$meanAngle)),3),
  curvature = round(curvature,2),                                          # curvature (1/mm) at 13 points (each of the bbKnots), defined as 1/r where r is the radius of curvature
  x=round(stageData$x*umPerPx),                                            # worm center-of-mass location in um
  y=round(stageData$y*umPerPx),                                            # worm center-of-mass location in um
  xHead = round(xHead*umPerPx),
  yHead = round(yHead*umPerPx),
  distToFoodHead = round(distToFoodHead*umPerPx),                          # distance to food boundary in um (from worm head)
  distToFoodBody = round(distToFoodBody*umPerPx),                          # distance to food boundary in um (from worm body)
  DMPevent = eventsToTimeseries(DMPs, length(wormPostures)),               # binary, whether this frame is a DMP event
  eggEvent = eventsToTimeseries(eggEvents, length(wormPostures), numEggs), # how many eggs were laid in this frame
  eggLoc = eventsToTimeseries(positionEggDetected, length(wormPostures), numEggs),   # which of 12 (6 on each side) body segments was the egg detected at 
  noseLift = 1:length(wormPostures),                        # hacking to skip this
  pumpingRate = pumpingRate,                                               # pumping rate, in Hz
  comSpeed1 = round(comSpeed1*umPerPx),                                    # absolute speed in um/s, using displacement over 1 second
  comSpeed10 = round(comSpeed10*umPerPx),                                  # absolute speed in um/s, using displacement over 10 seconds
  bodyAxisSpeed1 = round(bodyAxisSpeed1*umPerPx),                          # speed projected onto body axis in um/s, using displacement over 1 second
  bodyAxisSpeed10 = round(bodyAxisSpeed10*umPerPx),                        # speed projected onto body axis in um/s, using displacement over 10 seconds
  omegaTurn = sapply(wormPostures, function(x) x$error),                   # binary, was the worm in a self-intersecting posture in this frame?
  autofocusing = stageData$autofocusing,                                   # binary, was an autofocus occurring in this frame?
  laserOn = stageData$laserOn                                              # binary, was the laser on in this frame?
)

experimentName = strsplit(jpegDirectory,"/")[[1]]
experimentName = experimentName[length(experimentName)-1]

write.table(allBehaviorMatrix, file=jpegDirectory %+% "../"%+%experimentName%+%"_fullDataTable.csv", sep='\t', quote=FALSE)

save(list=setdiff(ls(),c("intensity","intensity.hp","dW")), file=jpegDirectory %+% "../analyzedFramesCompleteWithEggs.RData")


