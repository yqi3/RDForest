#######################################
#
# Functions to calculate difference distance metrics
# Based on Banerjee, 2005, Biometrika, "Geodetic Computations in Spatial Modeling"
#
######################################

# Calculates Euclidean distance
eudis <- function (x, y) {
  return(sqrt(sum((x - y)^2)))
}

# Converts degrees to radians
rad <- function(x) {
  return(x * pi/180)
}

# Remember:
# degree = rad * 180/pi
# rad    = degree * pi/180

# Calculates Euclidean, chord and geodetic distance between two points in longitud-latituted scale
dis <- function(lat1, lon1, lat2, lon2, R=6371) {
  
  # create three-dimensional vectors from the latitud-longitud information
  u1 <- c(R * cos(rad(lat1)) * cos(rad(lon1)), R * cos(rad(lat1)) * sin(rad(lon1)), R * sin(rad(lat1)))
  u2 <- c(R * cos(rad(lat2)) * cos(rad(lon2)), R * cos(rad(lat2)) * sin(rad(lon2)), R * sin(rad(lat2)))

  # Naive Euclidean distance between both points
  naive = eudis(c(lon1,lat1), c(lon2,lat2)) * pi * R/180

  # Chord distance
  chord = eudis(u1,u2)

  geodetic = R * acos(sin(rad(lat1))*sin(rad(lat2))+cos(rad(lat1))*cos(rad(lat2))*cos(rad(lon2-lon1)))

  return(list(naive = naive, chord = chord, geodetic=geodetic))
}

disvec <- function(lat1, lon1, lat2, lon2, R=6371) {
  # lat1: scalar
  # lon1: scalar
  # lat2: 1xn vector
  # lon2: 1xn vector
  library(fields)
  
  # create three-dimensional vectors from the latitud-longitud information
  u1 <- matrix(c(R * cos(rad(lat1)) * cos(rad(lon1)), R * cos(rad(lat1)) * sin(rad(lon1)), R * sin(rad(lat1))),
               nrow = 1, ncol = 3)
  u2 <- matrix(c(R * cos(rad(lat2)) * cos(rad(lon2)), R * cos(rad(lat2)) * sin(rad(lon2)), R * sin(rad(lat2))),
                  nrow = length(lat2), ncol = 3, byrow = FALSE)

  # Naive Euclidean distance between a point (lon1,lat1) and n different points, stored in matrix
  # where each row is the location of each point
  A = matrix(c(lon1, lat1), nrow=1,ncol=2, byrow=TRUE)
  B = matrix(c(lon2, lat2), nrow=length(lon2),ncol=2, byrow=FALSE)
  naive = rdist(A, B) * pi * R/180 # 1xn matrix

  # Chord distance
  chord = rdist(u1,u2)  #1xn matrix

  geodetic = R * acos(sin(rad(lat1))*sin(rad(lat2))+cos(rad(lat1))*cos(rad(lat2))*cos(rad(lon2-lon1))) 

  return(list(naive = naive, chord = chord, geodetic=geodetic))
}

latlon.to.u1u2 <- function(lat1, lon1, lat2, lon2, R=6371) {

  # lat1, lon1: scalars
  # lat2, lon2: vectors
  
  # create three-dimensional vectors from the latitud-longitud information
  u1 <- matrix(c(R * cos(rad(lat1)) * cos(rad(lon1)), R * cos(rad(lat1)) * sin(rad(lon1)), R * sin(rad(lat1))),
               nrow = 1, ncol = 3)
  u2 <- matrix(c(R * cos(rad(lat2)) * cos(rad(lon2)), R * cos(rad(lat2)) * sin(rad(lon2)), R * sin(rad(lat2))),
                  nrow = length(lat2), ncol = 3, byrow = FALSE)
  return(list(u1=u1, u2=u2))
}


# Triangular or Edge kernel function
Kt = function(u) {
  (1-abs(u)) * (abs(u)<= 1)
}
