########################################################################################################################
############################Basic geometric functions for Fastpop R script##############################################
########################################################################################################################
extractPCAMean <- function(knownPCA, numVector, rowNames){
  centroidCoord <- data.frame("ceu"=rep(NA, 3), "chb"=rep(NA, 3), 
                              "yri"=rep(NA, 3), "na"=rep(NA, 3), row.names = rowNames)
  for (j in 1:length(centroidCoord)){ #iterate over all populations for x, y, z coordinates
    for (i in 1:length(rowNames)){
      centroidCoord[i, j] <- mean(knownPCA[(numVector[j]):((numVector[j+1])-1), i])
    }
  }
  return(centroidCoord)
}
crossProduct <- function(a, b){
  return(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]))
}

dotProduct <- function(sampleVectorA, sampleVectorB){ #returns dot product for 2 vectors
  return((sampleVectorA[1] * sampleVectorB[1]) + (sampleVectorA[2] * sampleVectorB[2]) + 
           (sampleVectorA[3] * sampleVectorB[3]))
}

checkDistBetweenPoints<-function(point1, point2){ #finds distance between two 3D points
  return(sqrt((point2[1] - point1[1])^2 + (point2[2] - point1[2])^2 + 
                (point2[3] - point1[3])^2 )) 
}
findMagnitude <- function(sampleVector){
  return(sqrt(sampleVector[1]^2 + sampleVector[2]^2 + sampleVector[3]^2))
}

#find the distance between orig point and projected point on line
#(see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html)
distanceFromPointToLine <-function(linePt1, linePt2, samplePt){
 d <- findMagnitude(crossProduct(samplePt-linePt1,samplePt-linePt2))/
      findMagnitude(linePt2-linePt1)
}

definePlane<-function(vector1, vector2, vector3){
  #find the coordinates of the normal vector (a,b,c) to the side 
  #(m,p,q) = (x2-x1,y2-y1,z2-z1)
  mpq <- c("m" = vector2[1] - vector1[1], "p" = vector2[2] - vector1[2], "q" = vector2[3] - vector1[3])
  #(alpha,beta,gamma) = (x3-x1,y3-y1,z3-z1)
  abg <- c("alpha"=vector3[1] - vector1[1],"beta"=vector3[2] - vector1[2], "gamma"=vector3[3] - vector1[3])
  #use <m,p,q>/<alpha,beta,gamma> cross product and scalar plane equation, a(x-x0) + b(y-yo) + c(z-z0) - d = 0
  ijk <- crossProduct(mpq, abg)
  d0 <- -((ijk[1] * vector1[1]) + (ijk[2]*vector1[2]) + (ijk[3]*vector1[3]))
  #divide by c (ijk[3]);  c/c = 1
  a <- ijk[1]/ijk[3]; b <- ijk[2]/ijk[3]; c <- 1; d <- d0/ijk[3]
  return(c(a,b,c,d))
}

projectPointToPlane<-function(planeName, origPoint){
  #find the projection of the point on the closest side by solving for t 
  #t = parameter for the line equation, its direction vector being the normal vector from the side
  f.t <- function(t) { 
    t*(plane.vectors.line[[planeName]][1]^2 + 
         plane.vectors.line[[planeName]][2]^2 + 
         plane.vectors.line[[planeName]][3]^2) + 
      plane.vectors.line[[planeName]][1]*(origPoint["x"]) + 
      plane.vectors.line[[planeName]][2]*(origPoint["y"]) + 
      plane.vectors.line[[planeName]][3]*(origPoint["z"]) + 
      plane.vectors.line[[planeName]][4] 
  }
  t <- uniroot(f.t, c(-1000,1000))$root
  tempFellow <- c(
    origPoint["x"] + plane.vectors.line[[planeName]][1]*t, 
    origPoint["y"] + plane.vectors.line[[planeName]][2]*t, 
    origPoint["z"] + plane.vectors.line[[planeName]][3]*t)
  names(tempFellow)<-c("x","y","z")
  return(tempFellow)
}

createProjPointToLine<-function(directionVectors, vectorName, tempCentroidVector, samplePoint){ 
  #creates a point on a line where samplePoint is perpindicular (minimum distance)
  #use parametric form to find line equation
  a <- directionVectors[[vectorName]][1]
  b <- directionVectors[[vectorName]][2] 
  c <- directionVectors[[vectorName]][3] 
  #find unknown in equation, t
  f.t <- function(t) { 
    t*(a^2 + b^2 + c^2) + 
      a*(tempCentroidVector[1] - samplePoint[1]) + 
      b*(tempCentroidVector[2] - samplePoint[2]) + 
      c*(tempCentroidVector[3] - samplePoint[3]) 
  }
  t <- uniroot(f.t, c(-10000,10000))$root
  
  #use (x0,y0,z0), (a,b,c), and t to calculate new point
  #x = x0 + ta, y = y0 + tb, z = z0 + tc
  tempFellow <- c(tempCentroidVector[1] + a*t, tempCentroidVector[2] + b*t, tempCentroidVector[3] + c*t) 
  names(tempFellow)<-xyzVector
  return(tempFellow) #returns (x,y,z)
}

calcDirectionVectors <- function(lineList, directionVectorsList){
  for (n in 1:length(lineList)){
    for (m in 1:length(row.names(directionVectorsList))){
      tempPointsInLine <- unlist(strsplit(lineList[n], split = "\\."))
      directionVectorsList[m,n] <- centroidCoord[m, tempPointsInLine[2]] -centroidCoord[m, tempPointsInLine[1]] 
    }
  }
  return(directionVectorsList)
}

checkDistanceFromPlane<-function(allPlaneEquations, planeLine, samplePoint){  
  #distance = |dotproduct of v and n| = |a(x-x0)+b(y-y0)+c(z-z0)+d|/sqrt(a^2 + b^2 + c^2)
  #v = vector between samplePoint and point on plane, n =  unit normal vector
  #minimum distance from plane is length of line perpindicular to plane that passes through samplePoint
  a<-allPlaneEquations[[planeLine]][1]
  b<-allPlaneEquations[[planeLine]][2]
  c<-allPlaneEquations[[planeLine]][3]
  d<-allPlaneEquations[[planeLine]][4]
  return (abs( (a*samplePoint["x"]) + (b*samplePoint["y"]) + (c*samplePoint["z"]) + d) / (sqrt(a^2 + b^2 + c^2)))
}

#use determinant sign comparison to test if point is inside/outside tetrahedron 
#(http://steve.hollasch.net/cgindex/geometry/ptintet.html)
checkForInsideOrOutsideTetrahedron<-function(samplePoint, sampleDetMatrix){
  detMatrix4 <- detMatrix3 <- detMatrix2 <- detMatrix1 <- sampleDetMatrix
  detMatrix1[1, 1:3] <- samplePoint[1:3] #fills first 3 values of first row with (x,y,z) of sample
  detMatrix2[2, 1:3] <- samplePoint[1:3] #fills second row of second matrix
  detMatrix3[3, 1:3] <- samplePoint[1:3]
  detMatrix4[4, 1:3] <- samplePoint[1:3]
  
  detValues <- c(det(sampleDetMatrix), det(detMatrix1), det(detMatrix2),det(detMatrix3),det(detMatrix4))
  signOfdetValues <- sign(detValues) #reduces each determinant value to either -1 for neg, 0 for zero, or 1 for positive
  
  if (signOfdetValues[1] == signOfdetValues[2] && signOfdetValues[2] == signOfdetValues[3] &&
      signOfdetValues[3] == signOfdetValues[4] && signOfdetValues[4] == signOfdetValues[5]) {
    return('inside') #only if signs are either all positive or all negative
  }
  else{return('outside')}
}

#use barycentric coordinates to test if projected point is inside/outside of triangular plane
#(http://blackpawn.com/texts/pointinpoly/)
checkForInsideOrOutsideTriangle <- function(samplePoint, samplePlane){
  trianglePoints <- list("a" = c(0,0,0), "b" = c(0,0,0), "c" = c(0,0,0))
  pointList <- unlist(strsplit(x = samplePlane, split = "\\.")) #splits samplePlane("pt1.pt2.pt3") into separate strings
  for (i in (1:length(trianglePoints))){
    trianglePoints[[i]] <- centroidCoord[[pointList[i]]] #places x,y,z variables of specified point 
  }
  v0 <- trianglePoints[["c"]] - trianglePoints[["a"]] #create 3 vectors from point "a"
  v1 <- trianglePoints[["b"]] - trianglePoints[["a"]]
  v2 <- samplePoint - trianglePoints[["a"]]
  #find 5 possible dot products
  dot00 <- dotProduct(v0, v0); dot01 <- dotProduct(v0, v1); 
  dot02 <- dotProduct(v0, v2); dot11 <- dotProduct(v1, v1); dot12 <- dotProduct(v1, v2)
  #compute barycentric coordinates and determine if inside/outside
  invDenom <- 1 / (dot00 * dot11 - dot01 * dot01)
  u <- (dot11 * dot02 - dot01 * dot12) * invDenom
  v <- (dot00 * dot12 - dot01 * dot02) * invDenom
  
  if (u > 0 && v > 0 && (u + v <1)){
    return('inside')}
  else{return('outside')}  #will return outside if the point is on vertex or on edge
}

checkForInsideorOutsideSegment <- function(linePt1Name, linePt2Name, samplePt){
  linePt1 <- unlist(centroidCoord[linePt1Name])
  linePt2 <- unlist(centroidCoord[linePt2Name])
  lineDistance <- checkDistBetweenPoints(linePt1, linePt2)
  sampleToPt1 <- checkDistBetweenPoints(linePt1, samplePt)
  sampleToPt2 <- checkDistBetweenPoints(linePt2, samplePt)

  if (lineDistance < sampleToPt1){ #projPt is outside segment, closer to pt2
    return(linePt2Name)}
  if (lineDistance < sampleToPt2){ #projPt is outside segment, closer to pt1
    return(linePt1Name)}
  else{ # projPt is inside segment
    return("inside")}}

getProjPointsOnLines <- function(planeOfInterest, projPointOnPlane){
  projPointsOnLines <- data.frame(rep(NA,3), rep(NA,3), rep(NA,3), row.names = xyzVector)
  lineNames <- linesInPlanes[[planeOfInterest]]
  names(projPointsOnLines) <- lineNames
  for(i in 1:ncol(projPointsOnLines)){
    refCentroid <- unlist(strsplit(lineNames[[i]], split = "\\."))[1] #arbitrarily takes the first point from 2-point line
    projPointsOnLines[[i]] <- createProjPointToLine(directionVectors, lineNames[[i]], centroidCoord[[refCentroid]], projPointOnPlane)
  }
  return(projPointsOnLines)
}