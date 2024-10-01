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
  xFilt = caTools::runmean(x, k=k, ...)
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

