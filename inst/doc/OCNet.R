## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----lib, echo = FALSE--------------------------------------------------------
suppressPackageStartupMessages(library(OCNet))
suppressPackageStartupMessages(library(spam))
suppressPackageStartupMessages(library(igraph))


## ---- fig.width=7, fig.height=2-----------------------------------------------
set.seed(1)
OCN <- create_OCN(30, 20, outletPos = 1)
OCN <- aggregate_OCN(landscape_OCN(OCN), thrA = 3)
par(mfrow = c(1, 3), mai = c(0, 0, 0.2, 0.2))
draw_simple_OCN(OCN, thrADraw = 3)
title("Optimal Channel Network")
draw_elev3D_OCN(OCN, drawRiver = FALSE, addColorbar = FALSE, expand = 0.2, theta = -30)
title("Elevation")
draw_thematic_OCN(OCN$AG$streamOrder, OCN, discreteLevels = TRUE, colPalette = rainbow(4))
title("Strahler stream order")


## ----ex-net, echo=FALSE, fig.cap="Representation of the different aggregation levels at which the network is defined (excluding the null levels N4 and N8). The example is obtained from a single-outlet 8x8 lattice. Letter 'O' identifies the outlet pixel. Arrows on the other pixels identify flow directions. Numbers represent the cumulative drainage area (in number of pixels). At the FD level, all 64 pixels belong to the network^[Note that this pattern of flow directions was not derived by an OCN search algortihm, but rather drawn manually for illustration purposes.]. To obtain the RN level, a threshold value of 5 on drainage area is applied to distinguish pixels belonging to the river network. The network at the AG level consists of 9 nodes. The SC level is obtained by splitting the lattice into portions whose pixels drain into the AG nodes and edges. In this example, because there is only one outlet, all pixels belong to a single node at the CM level.", out.width = '70%'----
knitr::include_graphics("example_networks.png")

## ----ex3-net, echo=FALSE, fig.cap="Aggregation of the previous network with `A_thr` equal to 5 pixels and `maxReachLength` equal to 3 pixel sides. Note that the length of the diagonal segments of the edges is equal to `sqrt(2)`.", out.width = '70%'----
knitr::include_graphics("example_networks3.png")

## ----ex-ind, echo=FALSE, fig.cap="Top-left panel: arrows display flow directions, numbers identify contributing areas. Bottom-left panel: the network is aggregated by imposing a threshold equal to 2 pixels. The other panels display indices of nodes at the different aggregation levels. ", out.width = '90%'----
knitr::include_graphics("indices_example.png")

## -----------------------------------------------------------------------------
ex <- aggregate_OCN(landscape_OCN(OCN_4), thrA = 2)

ex$FD$toRN
ex$FD$toSC

ex$RN$toFD
ex$RN$toAG
ex$RN$toAGReach

ex$AG$toFD
ex$AG$ReachToFD
ex$AG$toRN
ex$AG$ReachToRN

ex$SC$toFD


## ---- fig.width=4, fig.height=4-----------------------------------------------
set.seed(1)
OCNwe <- create_OCN(20, 20, outletPos = 3, cellsize = 500)
par(mai=c(0,0,0,0))
draw_simple_OCN(OCNwe)


## ---- fig.width=7, fig.height=2.8---------------------------------------------
OCNwe <- landscape_OCN(OCNwe, slope0 = 0.01)
par(mai=c(0,0,0,0.5))
draw_elev3D_OCN(OCNwe, drawRiver = FALSE)


## -----------------------------------------------------------------------------
thr <- find_area_threshold_OCN(OCNwe)
# find index corresponding to thr$Nnodes ~= 20
indThr <- which(abs(thr$nNodesAG - 20) == min(abs(thr$nNodesAG - 20)))
indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
thrA20 <- thr$thrValues[indThr] # corresponding threshold area

## ---- fig.width=7, fig.height=2.8---------------------------------------------
OCNwe <- aggregate_OCN(OCNwe, thrA = thrA20)
par(mai=c(0.1,0,0.1,0))
draw_subcatchments_OCN(OCNwe)
points(OCNwe$AG$X,OCNwe$AG$Y, pch = 21, col = "blue", bg = "blue")


## ---- fig.width=7, fig.height=2.8---------------------------------------------
OCNwe <- paths_OCN(OCNwe, includePaths = TRUE)
par(mai=c(0.1,0,0.1,0))
draw_thematic_OCN(OCNwe$RN$downstreamPathLength[ , OCNwe$RN$outlet], OCNwe, 
                  backgroundColor = "#606060")


## -----------------------------------------------------------------------------
## Input data
OCNwe <- rivergeometry_OCN(OCNwe, widthMax = 5)   # evaluate river width 
K <- 10*OCNwe$RN$width                             # calculate carrying capacity 
pop0 <- 2*mean(K)*runif(OCNwe$RN$nNodes)           # initial random population vector
nTimestep <- 100                                   # number of timesteps
r <- 1.05                                          # proliferation rate
pd <- 0.5                                          # probability to move downstream
pu <- 1 - pd                                       # probability to move upstream
Go <- 5                                            # parameter controlling mobility 
# (no. individuals exiting from outlet node at carrying capacity is pu*Go) 

## -----------------------------------------------------------------------------
## Weights for upstream movement
Y <- rep(1,OCNwe$RN$nNodes)                    
for (i in 1:OCNwe$RN$nNodes){
  if (i != OCNwe$RN$outlet){
    Y[i] <- OCNwe$RN$A[i]/(OCNwe$RN$W[ , OCNwe$RN$downNode[i]] %*% OCNwe$RN$A)
  }
}

## -----------------------------------------------------------------------------
## Evaluate expected number of individuals moving at carrying capacity
GKK <- rep(0, OCNwe$RN$nNodes)
for (i in (1:OCNwe$RN$nNodes)){
  path <- OCNwe$RN$downstreamPath[[i]][[OCNwe$RN$outlet]] # select path
  prod <- Go                                                # initialize product of Y 
  for (j in path){
    prod <- prod*Y[j]
  }
  GKK[i] <- (pu/pd)^(length(path))*prod  
}

## -----------------------------------------------------------------------------
## Run metapopulation model
pop <- matrix(data=0,ncol=nTimestep,nrow=OCNwe$RN$nNodes)  # metapopulation matrix
pop[,1] <- pop0                                              # initialization
for (t in 2:nTimestep){
  for (i in 1:OCNwe$RN$nNodes){
    pop[i, t] <- 
      # Beverton-Holt growth model
      r*pop[i, t-1]/(1 + pop[i, t-1]*(r-1)/K[i]) +
      # individuals exiting from node i
                - (pu*(sum(OCNwe$RN$W[ , i])>0) + pd*(sum(OCNwe$RN$W[i, ])>0)) * 
      GKK[i] * (pop[i,t-1]/K[i]) +
      # individuals entering in i from the upstream nodes
                + pd * OCNwe$RN$W[ , i] %*% (GKK*pop[ , t-1]/K) +
      # individuals entering in i from the downstream node
                + pu * Y[i] * OCNwe$RN$W[i, ] %*% (GKK*pop[ , t-1]/K) 
    
  }
}

## ---- fig.width=7, fig.height=2.8---------------------------------------------
par(mfrow = c(1, 2))
plot(pop[OCNwe$RN$outlet, ], type = "l", ylim = c(0, 1.05*K[OCNwe$RN$outlet]), col = "red", 
     xlab = "Time", ylab = "Population", lwd = 2)
title("Evolution of local pop. size")
lines(c(1, nTimestep),c(K[OCNwe$RN$outlet], K[OCNwe$RN$outlet]), col = "red", lty = 2)
farthestNode <- which(OCNwe$RN$downstreamPathLength[ , OCNwe$RN$outlet]
                      == max(OCNwe$RN$downstreamPathLength[ , OCNwe$RN$outlet]))[1]
lines(pop[farthestNode, ], type="l", col="blue",lwd=2)
lines(c(1, nTimestep), c(K[farthestNode], K[farthestNode]), col = "blue", lty = 2)

plot(colSums(pop), type = "l", xlab = "Time", ylab = "Population", lwd = 2, ylim = c(0, 1.05*sum(K)))
lines(c(1, nTimestep), c(sum(K),sum(K)), lty = 2)
title("Evolution of metapop. size")

## ---- fig.width=7, fig.height=5-----------------------------------------------
par(mfrow = c(2, 2), mai = c(0.1, 0, 0.2, 0))
draw_thematic_OCN(pop[,1], OCNwe, colLevels = c(0, max(K), 1000),
                  drawNodes = TRUE)
title("Time = 1")
draw_thematic_OCN(pop[,5], OCNwe, colLevels = c(0, max(K), 1000),
                  drawNodes = TRUE)
title("Time = 5")
draw_thematic_OCN(pop[,20], OCNwe, colLevels = c(0, max(K), 1000),
                  drawNodes = TRUE)
title("Time = 20")
draw_thematic_OCN(pop[,100], OCNwe, colLevels = c(0, max(K), 1000),
                  drawNodes = TRUE)
title("Time = 100")


## ---- fig.width=7, fig.height=5-----------------------------------------------
par(mfrow = c(2, 3), mai = c(0, 0, 0.2, 0))
peano0 <- create_peano(0)
draw_simple_OCN(peano0)
title("Iteration: 0 - Lattice size: 2x2")

peano1 <- create_peano(1)
draw_simple_OCN(peano1)
title("Iteration: 1 - Lattice size: 4x4")

peano2 <- create_peano(2)
draw_simple_OCN(peano2)
title("Iteration: 2 - Lattice size: 8x8")

peano3 <- create_peano(3)
draw_simple_OCN(peano3)
title("Iteration: 3 - Lattice size: 16x16")

peano4 <- create_peano(4)
draw_simple_OCN(peano4)
title("Iteration: 4 - Lattice size: 32x32")

peano5 <- create_peano(5)
draw_simple_OCN(peano5)
title("Iteration: 5 - Lattice size: 64x64")


## ---- fig.width=5, fig.height=3.5---------------------------------------------
par(mai = c(0, 0, 0, 0))
peano5 <- landscape_OCN(peano5)
draw_elev2D_OCN(peano5)


## ---- fig.width=5, fig.height=4-----------------------------------------------
par(mai=c(0.1,0.1,0.1,0.1))
g <- OCN_to_igraph(OCNwe, level = "AG")
plot.igraph(g, vertex.color = rainbow(OCNwe$AG$nNodes), 
     layout = matrix(c(OCNwe$AG$X,OCNwe$AG$Y),ncol = 2, nrow = OCNwe$AG$nNodes))

## ---- fig.width=5, fig.height=3-----------------------------------------------
par(mai=c(0,0,0,0))
draw_thematic_OCN(c(1:OCNwe$AG$nNodes), OCNwe, discreteLevels = TRUE, drawNodes = TRUE,
                  colPalette = rainbow,  cex = 3, riverColor = "#999999",
                  backgroundColor = "#00000000", addLegend = FALSE)
text(OCNwe$AG$X, OCNwe$AG$Y)


