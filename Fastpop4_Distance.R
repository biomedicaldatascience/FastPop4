#########################################################################################################
#   Probability Inference on Distance-based method using Principal Components     - updated  4/24/2019  #
#   created by Jinyoung Byun, Kevin Xu, and Chris Amos                                                  #
#########################################################################################################

options(digits = 15) #increase the default number of decimal digits
source("Fastpop4_functions.R")

ttt <- read.table(file = "fastpop4_pca.discovery.txt", header=FALSE, sep = " ")

pcaBoundaries <- c(1, 166, 303, 506, 549)#denotes lower boundaries of ceu, chb, yri, na, and unknown samples
numKnown <- tail(pcaBoundaries,1) - 1
knownPca <- ttt[1:numKnown,1:3]
scores.svd <- read.table(file="new.pca.txt", header=FALSE, sep=" ")
k<- length(scores.svd[,1]);k
rownames(scores.svd)=scores.svd[,2]


################################# final result variable #####################################
results <- data.frame("ID" = (1:k),"ceu"=NA,"chb"=NA,"yri"=NA,"na"=NA) 

####################################### coordinates extraction ##############################
#############################################################################################
xyzVector <- c('x', 'y', 'z')

#call function from "Fastpop4_functions.R"
centroidCoord <- extractPCAMean(knownPca, pcaBoundaries, xyzVector) #returns dataframe with columns = populations, rows = xyz
centroidNames <- names(centroidCoord) #CEU = Europe, CHB = Chinese, YRI = African, NA = Native American

#######create dataframe to hold all possible points (alphabetized) in all possible planes########
pointsInPlanes <- data.frame("ceu.chb.na"=rep(NA,3),"ceu.na.yri"=rep(NA,3),"ceu.chb.yri"=rep(NA,3),
                             "chb.na.yri"=rep(NA,3),row.names = c("pt1","pt2","pt3"))
pointsInPlanes[1:3,] <- c("ceu.chb.na" = unlist(strsplit("ceu.chb.na", split = "\\.")),
                          "ceu.na.yri"= unlist(strsplit("ceu.na.yri", split = "\\.")),
                          "ceu.chb.yri" = unlist(strsplit("ceu.chb.yri", split = "\\.")),
                          "chb.na.yri" = unlist(strsplit("chb.na.yri", split = "\\."))) 

#############create dataframe to hold all possible lines in all possible planes##################
linesInPlanes <- data.frame("ceu.chb.na"=rep(NA,3),"ceu.na.yri"=rep(NA,3),"ceu.chb.yri"=rep(NA,3),
                            "chb.na.yri"=rep(NA,3),row.names = c("line1","line2","line3"))
for (i in (1:ncol(pointsInPlanes))){
  linesInPlanes[1:3,i] <- c(paste(pointsInPlanes[1,i],pointsInPlanes[2,i],sep="."), #centroids 1/2
                            paste(pointsInPlanes[1,i],pointsInPlanes[3,i],sep="."), #centroids 1/3
                            paste(pointsInPlanes[2,i],pointsInPlanes[3,i],sep="."))} #centroids 2/3

############create dataframe to hold the 2 points in all possible lines##################
lineNames <- unlist(linesInPlanes)
pointsInLines <- as.data.frame(matrix(NA, ncol = 12, nrow = 2)) #12 lines, 2 points in each line
colnames(pointsInLines) <- lineNames
for (i in 1:ncol(pointsInLines)){
  pointsInLines[[i]] <- unlist(strsplit(lineNames[i], split = "\\."))
}
######################list to hold distance between the centroids ######################
distance.centroids <- list (
  "ceu.chb"=checkDistBetweenPoints(centroidCoord[["ceu"]], centroidCoord[["chb"]]),
  "ceu.yri"=checkDistBetweenPoints(centroidCoord[["ceu"]], centroidCoord[["yri"]]),
  "ceu.na"=checkDistBetweenPoints(centroidCoord[["ceu"]], centroidCoord[["na"]]),
  "chb.na"=checkDistBetweenPoints(centroidCoord[["chb"]], centroidCoord[["na"]]),
  "na.yri"=checkDistBetweenPoints(centroidCoord[["na"]], centroidCoord[["yri"]]),
  "chb.yri"=checkDistBetweenPoints(centroidCoord[["chb"]], centroidCoord[["yri"]]))

###################### direction vectors of those line ######################
directionVectors <- data.frame("ceu.chb"=rep(NA,3), "ceu.yri"=rep(NA,3), "ceu.na"=rep(NA,3),
                               "chb.na"=rep(NA,3), "na.yri"=rep(NA,3), "chb.yri"=rep(NA,3), row.names=c("a","b","c"))
#return filled dataframe using function from "Fastpop_functions.R"
directionVectors <- calcDirectionVectors(names(directionVectors), directionVectors)

################################# Define the planes ######################################
#find the one-line coordinates of each plane and put them into a dataframe
plane.vectors.line <- data.frame("ceu.chb.na"=rep(NA,4),"ceu.na.yri"=rep(NA,4),"ceu.chb.yri"=rep(NA,4),"chb.na.yri"=rep(NA,4),row.names = c("a","b","c","d"))
planeNames <- names(plane.vectors.line) #stores names of planes

#planes: ceu/chb/na, ceu/na/yri, ceu/chb/yri, chb/na/yri
plane.vectors.line["ceu.chb.na"] <- definePlane(centroidCoord[["ceu"]], centroidCoord[["chb"]], centroidCoord[["na"]])
plane.vectors.line["ceu.na.yri"] <- definePlane(centroidCoord[["ceu"]], centroidCoord[["na"]], centroidCoord[["yri"]])
plane.vectors.line["ceu.chb.yri"] <- definePlane(centroidCoord[["ceu"]], centroidCoord[["chb"]], centroidCoord[["yri"]])
plane.vectors.line["chb.na.yri"] <- definePlane(centroidCoord[["chb"]], centroidCoord[["na"]], centroidCoord[["yri"]])

##################create template matrix for use in checking inside/outside tetrahedron#######
detMatrix0 <- matrix(c(centroidCoord[["ceu"]], 1, centroidCoord[["chb"]], 1, 
                       centroidCoord[["na"]], 1, centroidCoord[["yri"]], 1), 
                       nrow = 4, ncol = 4, byrow = TRUE)

################################# Set up some helper funcs ######################################

getPointsInPlaneCoordinates <- function(sampleCentroidsInPlane, samplePlaneName){ #returns 3 sets of x,y,z coordinates for POI (plane of interest)
  tempcentroidsInPOICoord <- data.frame(rep(NA,3), rep(NA,3), rep(NA,3))
  colnames(tempcentroidsInPOICoord)<-pointsInPlanes[[samplePlaneName]][1:3] #get point names from dataframe that holds 3 corresponding points for each of the 4 planes
  rownames(tempcentroidsInPOICoord)<-xyzVector
  for (i in (1:ncol(tempcentroidsInPOICoord))){
    tempcentroidsInPOICoord[,i] <-centroidCoord[[toString(sampleCentroidsInPlane[i])]] #places only the 3 pertinent points' (x,y,z) coordinates
  }
  return(tempcentroidsInPOICoord)
}

checkClosestPlane<-function(samplePoint){
  distToPlanes <- data.frame("ceu.chb.na"=-1,"ceu.na.yri"=-1,"ceu.chb.yri"=-1,"chb.na.yri"=-1,row.names = "distance")
  #find distances from sample to centroids and sort
  centroidsDist <- c("ceu"=NA, "chb" = NA, "na" = NA, "yri" = NA)
  
  #for loop that indexes distance between sample and each centroid
  for (i in 1:length(centroidsDist)){
    centroidsDist[i] <- checkDistBetweenPoints(samplePoint, centroidCoord[[names(centroidsDist[i])]])
  }

  sortedCentroidsDist <- sort(centroidsDist)
  
  #closest plane will have closest two centroids, and one other centroid
  namesCloseCentroids<-names(head(sortedCentroidsDist, 2))
  namesFarCentroids<-names(tail(sortedCentroidsDist, 2))
  namesPlanes <- data.frame("plane1" = sort(append(namesCloseCentroids, namesFarCentroids[1])), 
                           "plane2" = sort(append(namesCloseCentroids, namesFarCentroids[2])))
  #check distance to plane for both possible planes, return closest plane
  namePlane1 <- paste(namesPlanes[1,1], namesPlanes[2,1], namesPlanes[3,1], sep = ".") # plane named as "pt1.pt2.p3"
  namePlane2 <- paste(namesPlanes[1,2], namesPlanes[2,2], namesPlanes[3,2], sep = ".")
  distToPlanes[namePlane1] <- checkDistanceFromPlane(plane.vectors.line,namePlane1, samplePoint) #will return positive distance
  distToPlanes[namePlane2] <- checkDistanceFromPlane(plane.vectors.line,namePlane2, samplePoint) 
  return(colnames(sort(distToPlanes)[3])) #since it's sorted, first two will be -1, 3rd will be closest plane
}

#if inside tetrahedron and inside closest plane, compute probabilities of 3 points on closest plane
checkInsideProbability<-function(planeName, pointOnPlane){
  pDivisor<-0
  #find the 3 individual points and 3 lines associated with closest plane
  centroidsInPOI <- c(pointsInPlanes[[planeName]][1:3]) #POI = Plane of Interest
  linesInPOI <- c(linesInPlanes[[planeName]][1:3]) #3 possible lines in POI (points 1.2, 1.3, 2.3)
  names(centroidsInPOI) <- names(linesInPOI) <- planeName
  #dataframe to hold coordinates of the plane's 3 centroids
  centroidsInPOICoord<- getPointsInPlaneCoordinates(centroidsInPOI, planeName)

  prob4Pop <- c("ceu" = 0, "chb" = 0, "yri" = 0, "na" = 0) #3 will be filled with values, 1 will be left 0
  
  #hold (x,y,z) coordinates for the 3 points projected from pointOnPlane to each edge of triangle
  #projPointsOnLines <- data.frame("samplePoint12"= rep(NA,3), "samplePoint13"=rep(NA,3), "samplePoint23"=rep(NA,3), row.names = xyzVector)
  projPointsOnLines <- data.frame(rep(NA,3), rep(NA,3), rep(NA,3), row.names = xyzVector)
  names(projPointsOnLines) <- linesInPOI
  
  #hold distance between proj. points on plane and point on triangle edges
  distBetweenProjPtsOnTri <- c("distToProjPt1" = NA, "distToProjPt2" = NA, "distToProjPt3" = NA)
  
  #get lengths of the plane's 3 lines
  distanceBetween12 <- distance.centroids[[toString((linesInPOI[1]))]] 
  distanceBetween13 <- distance.centroids[[toString((linesInPOI[2]))]]
  distanceBetween23 <- distance.centroids[[toString((linesInPOI[3]))]] 

  #loop to fill in projPointsOnLines and distBetweenProjPtsOnTri
  for (i in (1:length(projPointsOnLines))){ 
    pt <- unlist(strsplit(linesInPlanes[[planeName]][i], split = "\\."))[1] #arbitrarily takes the first point from 2-point line
    projPointsOnLines[[i]] <- createProjPointToLine(directionVectors, linesInPOI[i], centroidsInPOICoord[[pt]], pointOnPlane)
    distBetweenProjPtsOnTri[i] <- checkDistBetweenPoints(projPointsOnLines[[i]], pointOnPlane)
    
    pDivisor <- pDivisor + (1/distBetweenProjPtsOnTri[i]) #pDivisor + 1/distance
  }
  # it is *still* possible that projecting a point from inside the triangle to a vector results in an
  # intersection that is outside the triangle. In that case, l1, l3, and l5 will be larger than the
  # distance between the centroid pairs. To find l2, l4, l6 when that is true, we use the difference 
  # between the external point and the length of the vector is question.
  
  #l1 = distance between 1st projected samplePoint (on triangle edge) and centroid 1
  l1 <- checkDistBetweenPoints(projPointsOnLines[[1]], centroidsInPOICoord[,1]) 
  #l2 = line "pt1.pt2" length - l1
  l2 <- max(c(l1, distanceBetween12)) - min(c(l1, distanceBetween12)) #max-min
  
  l3 <- checkDistBetweenPoints(projPointsOnLines[[2]], centroidsInPOICoord[,1])
  l4 <- max(c(l3, distanceBetween13)) - min(c(l3, distanceBetween13))
  
  l5 <- checkDistBetweenPoints(projPointsOnLines[[3]], centroidsInPOICoord[,2])
  l6 <- max(c(l5, distanceBetween23)) - min(c(l5, distanceBetween23))

  if (l1 > distanceBetween12 || l3 > distanceBetween13 ||l5 > distanceBetween23) {
    return(probsForTwoClosestLinesOnly(planeName, pointOnPlane, projPointsOnLines))  
  }
  distTo_12 <- distBetweenProjPtsOnTri["distToProjPt1"]
  distTo_13 <- distBetweenProjPtsOnTri["distToProjPt2"]
  distTo_23 <- distBetweenProjPtsOnTri["distToProjPt3"]
  
  #assign probabilities to appropriate 3 centroids; probabilities are weighted by distance of samplePoint to each line
  prob4Pop[toString(centroidsInPOI[1])]<-(1/distTo_12*((1/l1)/(1/l1+1/l2)) + 1/distTo_13*((1/l3)/(1/l3+1/l4))) / pDivisor
  prob4Pop[toString(centroidsInPOI[2])]<-(1/distTo_12*((1/l2)/(1/l1+1/l2)) + 1/distTo_23*((1/l5)/(1/l5+1/l6))) / pDivisor
  prob4Pop[toString(centroidsInPOI[3])]<-(1/distTo_13*((1/l4)/(1/l3+1/l4)) + 1/distTo_23*((1/l6)/(1/l5+1/l6))) / pDivisor
  return(prob4Pop)
}
#if the point is inside the tetrahedron but projects outside of a triangle face, 
#find probability for only the two closest lines
probsForTwoClosestLinesOnly<-function(planeName, pointOnPlane, projPointsOn3Lines){
  probabilityFrame <- data.frame("ceu" = rep(0,3), "chb" = rep(0,3), "na" = rep(0,3), 
                                 "yri" = rep(0,3),row.names= c("line1P", "line2P", "finalP"))
  centroidsInPOI <- c(pointsInPlanes[[planeName]][1:3]) #POI = Plane of Interest
  centroidsInPOICoord<- getPointsInPlaneCoordinates(centroidsInPOI, planeName)
  
  position2Points<- c("NA", "NA")#stores if proj pt is within segment or not
  distanceTo3Lines <- c(rep(NA,3))
  names(distanceTo3Lines) <- names(projPointsOn3Lines)
  
  for (i in 1:length(distanceTo3Lines)){
    distanceTo3Lines[i] <- checkDistBetweenPoints(projPointsOn3Lines[[i]], pointOnPlane)
  }
  distanceTo3Lines <- sort(distanceTo3Lines) #keeps distance values only of two closest lines
  closestTwoLines <- names(head(distanceTo3Lines,2)) #names of two closest lines
  
  projPoints2Lines <- data.frame(rep(NA,5), rep(NA,5), row.names = c("x", "y", "z","l1","l2"))
  names(projPoints2Lines) <- closestTwoLines
  line1Coord <- c(NA,NA)
  line2Coord <- list(c("x"=NA, "y"=NA, "z"=NA), c("x"=NA, "y"=NA, "z"=NA))
  
  TwoLinesPtsNames<- data.frame(rep(NA,2),rep(NA,2)) #dataframe to hold names of points in line 1/2
  names(TwoLinesPtsNames) <- closestTwoLines
  for (i in 1:length(TwoLinesPtsNames)){
    TwoLinesPtsNames[[i]]<-pointsInLines[[closestTwoLines[i]]]

    Centroid1Coord <- centroidCoord[[TwoLinesPtsNames[1,i]]]; 
    Centroid2Coord <- centroidCoord[[TwoLinesPtsNames[2,i]]];

    #arbitrary selection of first point for coord
    projPoints2Lines[1:3,i] <-createProjPointToLine(directionVectors, closestTwoLines[i], 
                                                            Centroid1Coord, pointOnPlane)
    #check if point is inside or outside segment; return "inside" or centroid name it's closer to
    position2Points[i] <- checkForInsideorOutsideSegment(TwoLinesPtsNames[1,i], TwoLinesPtsNames[2,i],pointOnPlane)
    #if inside, calculate l1 and l2
    #l1 = distance between projected point on 2 closest lines and 1 of line's centroids (arbitrary)
    #l2 = difference between l1 and entire segment distance
    if(position2Points[i] =="inside"){
      refPt1 <- TwoLinesPtsNames[[i]][1] #use an arbitrary point in line (first centroid coord)
      otherPt2 <- TwoLinesPtsNames[[i]][2]
      projPoints2Lines['l1',i]<- checkDistBetweenPoints(projPoints2Lines[1:3,i], Centroid1Coord)
      minMax <- sort(c(projPoints2Lines['l1',i], distance.centroids[[closestTwoLines[1]]]))
      projPoints2Lines['l2',i]<- minMax[2]- minMax[1]
      probabilityFrame[i,refPt1] <- (1/projPoints2Lines['l1',i]) / (1/projPoints2Lines['l1',i] + 1/projPoints2Lines['l2',i])
      probabilityFrame[i,otherPt2] <- (1- probabilityFrame[i, refPt1])
    }
    else {
      probabilityFrame[i, position2Points[i]] <- 1 #set probability of closest centroid to 100%
    }
  }
  #weight = (1/l1)(1/l1 + 1/l2)
  line1DWeight <- (1/(distanceTo3Lines[1])/(1/(distanceTo3Lines[1])+1/(distanceTo3Lines[2])))
  line2DWeight <- (1/(distanceTo3Lines[2])/(1/(distanceTo3Lines[1])+1/(distanceTo3Lines[2])))
  for(i in 1:length(probabilityFrame)){
    probabilityFrame["finalP",i] <- (probabilityFrame[1,i]*line1DWeight) + (probabilityFrame[2,i]*line2DWeight)
  }
  return(unlist(probabilityFrame["finalP",]))
}

probsForTwoClosestCentroidsOnly<-function(planeName, pointOnPlane){
  probabilityVector <- c("ceu" = 0, "chb" = 0, "yri" = 0, "na" = 0)
  inOrOutSegment <- NA
  #find two closest centroids 
  distanceCentroidVector <- c("chb" = -1, "ceu" = -1, "yri" = -1, "na" = -1)
  
  #centroidsInPOI <- c(pointsInPlanes[[planeName]][1:3])
  #centroidsInPOICoord<- getPointsInPlaneCoordinates(centroidsInPOI, planeName)
  
  for (i in (1:ncol(centroidsInPOICoord))){ #for each point in POI, store distances from sample point to centroids
    tempCoord <- centroidsInPOICoord[,i]; names(tempCoord) <-xyzVector 
    distanceCentroidVector[names(centroidsInPOICoord[i])] <- checkDistBetweenPoints(pointOnPlane, tempCoord)
  }
  distanceCentroidVector<- sort(distanceCentroidVector) #1st value will be pt4 = -1, 4th value will be max 
  closestTwo <- sort(names(distanceCentroidVector[2:3])) #alphabetizes 2nd and 3rd values
  closestLine <- paste(closestTwo[1],".",closestTwo[2], sep = "") #creates "pt2.pt3"
  
  #take name of closest line, 1 of the line's 2 centroids' coordinates, and samplePoint;
  #find distance between the centroid and a projected point on the line where samplePoint is perpindicular
  perpPointOnLine <- createProjPointToLine(directionVectors, closestLine, centroidsInPOICoord[[toString(closestTwo[1])]], pointOnPlane)
  
  #check if point is inside or outside segment; return "inside" or centroid name it's closer to
  inOrOutSegment <- checkForInsideorOutsideSegment(closestTwo[1], closestTwo[2], perpPointOnLine)
  
  if (inOrOutSegment=="inside"){
    #l1 = distance between projected samplePoint on closest line and centroid 1 (arbitrary)
    l1 <- checkDistBetweenPoints(perpPointOnLine, centroidsInPOICoord[[toString(closestTwo[1])]])
    l2 <- max(c(l1, distance.centroids[[toString(closestLine)]])) - min(c(l1, distance.centroids[[toString(closestLine)]]))
    
    probabilityVector[toString(closestTwo[1])] <- (1/l1) / (1/l1 + 1/l2)
    probabilityVector[toString(closestTwo[2])] <- 1 - probabilityVector[toString(closestTwo[1])]
  }
  else{ #inOrOutSegment == name of closest centroid
    probabilityVector[inOrOutSegment] <- 1
  }
  return(probabilityVector)
}

for(j in 1:k){
  #load the sample coordinates (the three first scores from the PCA)
  sample <- c(scores.svd[j,3],scores.svd[j,4],scores.svd[j,5]) 
  names(sample) <- xyzVector
  
  #a few variables neet to be set by default
  inorout <- "NA" ; inorout2 <- "NA"
  P <- c("ceu"=NA,"chb"=NA,"yri"=NA,"na"=NA)
  sample.on.closest.plane <- NA
  
  #check inside/outside tetrahedron
  inorout<-checkForInsideOrOutsideTetrahedron(sample, detMatrix0)
  
  ############## if outside, find the closest plane and apply FastPop3d to the projected point ############
  if (inorout == "outside") {
    closest.plane<-checkClosestPlane(sample)
    centroidsInPOI <- c(pointsInPlanes[[closest.plane]][1:3]) #gets centroid names in plane of interest (POI)
    centroidsInPOICoord<- getPointsInPlaneCoordinates(centroidsInPOI, closest.plane) #gets coordinates
    
    #find the projection of the point on the closest side 
    sample.on.closest.plane<-projectPointToPlane(closest.plane, sample)
    #then find the position of the projected point to the plane to determine if it's inside the triangle or not.
    #projPointsOnLines <- data.frame()
    inorout2<-checkForInsideOrOutsideTriangle(sample.on.closest.plane, closest.plane)
    #if outside both tetrahedron and closest plane, find probabilities of only two closest centroids
    if(inorout2 == 'outside'){
      P <-probsForTwoClosestCentroidsOnly(closest.plane, sample.on.closest.plane)
    }
    #else, find probabilities of 3 centroids that comprise plane
    else{
      P <- checkInsideProbability(closest.plane, sample.on.closest.plane)
    }
    for (i in (1:length(P))){
      results[[centroidNames[i]]][j] <- P[centroidNames[i]]
    }
  }
 ####################if inside : apply FastPop3d with the projected of the point on each side of the tetrahedron #####
  else{ 
    #create dataframe to hold probability output (4 values) for each of the 4 planes
    p4planes <- data.frame("P1" = rep(NA,4), "P2" = rep(NA,4), "P3" = rep(NA,4), "P4"= rep(NA,4), row.names = centroidNames)
    #create variable to hold distance of point to each of 4 planes
    d4planes <- c("D1" = NA, "D2" = NA, "D3" = NA, "D4"= NA)
    for (i in (1:length(p4planes))){ 
      projPointOnPlane <- projectPointToPlane(planeNames[i], sample)
      test2 <- checkInsideProbability(planeNames[i], projPointOnPlane)[1:4]
      p4planes[[i]] <- checkInsideProbability(planeNames[i], projPointOnPlane)
      d4planes[i] <- checkDistBetweenPoints(sample, projPointOnPlane)
    }
    P <- c("ceu"=0,"chb"=0, "yri"=0, "na"=0)
    divisor4Pop <- (1/d4planes[1] + 1/d4planes[2] + 1/d4planes[3] + 1/d4planes[4])  
    for(i in (1:length(P))){
      centroid <-centroidNames[i]
      #add up probability for each centroid, weighted by distance of point to respective plane
      P[[centroid]] <- ((p4planes[centroid,1])*(1/d4planes[1]) + (p4planes[centroid,2])*(1/d4planes[2])
                        + (p4planes[centroid,3])*(1/d4planes[3]) + (p4planes[centroid,4])*(1/d4planes[4]))/divisor4Pop
      results[[centroid]][j] <- P[[centroid]]
    }
  }

}
write.table(results, file="fastpop4_ancestry.ot.csv",row.names=T,sep=',')