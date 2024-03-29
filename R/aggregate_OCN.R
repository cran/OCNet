
aggregate_OCN <- function(OCN,
                          thrA=0.002*OCN$FD$nNodes*OCN$cellsize^2,
                          streamOrderType="Strahler",
                          maxReachLength=Inf,
                          equalizeLengths=FALSE,
                          breakpoints=NULL,
                          displayUpdates=FALSE){
  
  if (!("slope" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to aggregate_OCN.')
  }
  
  if (maxReachLength < OCN$cellsize*sqrt(2)){
    stop("maxReachLength cannot be smaller than OCN$cellsize*sqrt(2).")
  }
  #t1 <- Sys.time()
  
  if (thrA==0) maxReachLength <- OCN$cellsize*sqrt(2)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # BUILD NETWORK AT RN LEVEL ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if(displayUpdates){message("Calculating network at RN level...      \r", appendLF = FALSE)}
  
  #print('Crop data at FD level to RN level...',quote=FALSE); 
  RN_mask <- as.vector(OCN$FD$A >= thrA)# RN_mask allows to sample RN-level values from matrices/vectors at FD level   
  RN_to_FD <- which(RN_mask) # RN_to_FD[i] is the pixel ID at the FD level of the pixel whose ID at the RN level is i
  FD_to_RN <- RN_mask*cumsum(as.numeric(RN_mask)) # FD_to_RN[i] is the pixel ID at the RN level of the pixel whose ID at the FD level is i
  # if pixel i at FD level doesn't belong to RN, then FD_to_RN[i]=0
  
  Nnodes_RN <- length(RN_to_FD)
  
  W_RN <- OCN$FD$W[RN_mask,,drop=FALSE]
  W_RN <- W_RN[,RN_mask,drop=FALSE]
  
  Outlet_RN <- FD_to_RN[OCN$FD$outlet]
  Outlet_RN <- Outlet_RN[Outlet_RN!=0] # remove outlets if the corresponding catchment size is lower than threshold
  DownNode_RN <- numeric(Nnodes_RN)
  # for (i in 1:Nnodes_RN){
  #   if (!(i %in% Outlet_RN)){
  #     DownNode_RN[i] <- which(W_RN[i,]==1)
  #   }}
  tmp <- W_RN@rowpointers
  NotOutlet <- which((tmp[-1] - tmp[-length(tmp)])==1)
  DownNode_RN[NotOutlet] <- W_RN@colindices
  
  # reverse downNode_RN
  DownNode_RN_rev <- vector("list",Nnodes_RN)
  for (i in 1:Nnodes_RN){
    d <- DownNode_RN[i]
    if (d!=0){DownNode_RN_rev[[d]] <- c(DownNode_RN_rev[[d]],i)  }}
  
  A_RN <- OCN$FD$A[RN_mask]
  X_RN <- OCN$FD$X[RN_mask]
  Y_RN <- OCN$FD$Y[RN_mask]
  Z_RN <- OCN$FD$Z[RN_mask]
  Length_RN <- OCN$FD$leng[RN_mask]
  
  # Drainage density
  DrainageDensity_RN <- sum(Length_RN)/(OCN$dimX*OCN$dimY*OCN$cellsize^2)
  
  # Connectivity indices at pixel level
  DegreeIn <- colSums(W_RN)
  DegreeOut <- rowSums(W_RN)
  Confluence <- DegreeIn>1
  Source <- DegreeIn==0
  SourceOrConfluence <- Source|Confluence
  ConfluenceNotOutlet <- Confluence&(DownNode_RN!=0)
  ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet
  
  OutletNotChannelHead <- (DownNode_RN==0)&(!ChannelHeads)
  IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
  IsNodeAG[breakpoints] <- TRUE
  whichNodeAG <- which(IsNodeAG)
  
  # Calculate slope for each pixel of the river network 
  Slope_RN <- OCN$FD$slope[RN_mask]
  
  # sort nodes in downstream direction
  ind_sort <- sort(A_RN, index.return=TRUE)
  ind_sort <- ind_sort$ix
  Upstream_RN <- vector("list",Nnodes_RN)
  Nupstream_RN <- numeric(Nnodes_RN)
  for (i in 1:Nnodes_RN){
    ups <- as.numeric(DownNode_RN_rev[[ind_sort[i]]])
    nodes <- numeric(0)
    for (u in ups){ nodes <- c(nodes, Upstream_RN[[u]])}
    Upstream_RN[[ind_sort[i]]] <- c(nodes, ind_sort[i])
    Nupstream_RN[ind_sort[i]] <- length(Upstream_RN[[ind_sort[i]]])
  }
  
  # RN_to_CM[i] indicates outlet to which reach i drains
  RN_to_CM <- numeric(Nnodes_RN)
  for (i in 1:OCN$nOutlet){
    RN_to_CM[Upstream_RN[[Outlet_RN[i]]]] <- i
  }
  if (displayUpdates){message("Calculating network at RN level... 100.0%\n", appendLF = FALSE)}
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # BUILD NETWORK AT AG LEVEL ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if(displayUpdates){message("Calculating network at AG level...      \r", appendLF = FALSE)}
  
  # Vector that attributes reach ID to all river network pixels
  #print('Define nodes of aggregated network...',quote=FALSE); 
  Nnodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(Nnodes_AG)
  RN_to_AG <- numeric(Nnodes_RN)
  AG_to_RN <- vector("list", Nnodes_AG)
  reachID <- 1
  X_AG <- NaN*numeric(Nnodes_AG)
  Y_AG <- NaN*numeric(Nnodes_AG)
  Z_AG <- NaN*numeric(Nnodes_AG)
  A_AG <- NaN*numeric(Nnodes_AG)
  
  # new version
  new_maxLength <- maxReachLength
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID
    AG_to_RN[[reachID]] <- i
    X_AG[reachID] <- X_RN[i]
    Y_AG[reachID] <- Y_RN[i]
    Z_AG[reachID] <- Z_RN[i]
    A_AG[reachID] <- A_RN[i]
    
    if (equalizeLengths){
      if (IsNodeAG[i]){ # if the node is pre-established AG node, re-estimate new_maxReachLength; otherwise, use previously calculated value
        j <- DownNode_RN[i]
        tmp_length <-  Length_RN[i]
        while (!IsNodeAG[j] && j!=0) {
          tmp_length <-  tmp_length + Length_RN[j]
          j_old <- j
          j <- DownNode_RN[j]}
        
        n_splits <- max(1, ceiling((tmp_length+1.5*OCN$cellsize)/maxReachLength)) # it works even if tmp_length < maxReachLength. min is 1, to prevent issues when i=outlet
                                                                                  # +1.5 cellsize to avoid unsplittable reaches
        new_maxLength <- min(tmp_length/n_splits + 1.5*OCN$cellsize, maxReachLength) # sum +1.5 cellsize to avoid creating additional reaches due to rounding
        new_maxLength <- max(new_maxLength, 1.5*OCN$cellsize)
      }
    } else {new_maxLength <- maxReachLength}
    
    lengReach <- Length_RN[i]
    j <- DownNode_RN[i] 
    j0 <- j
    tmp <- NULL
    
    while (!IsNodeAG[j] && j!=0 && lengReach <= new_maxLength) {
        tmp <- c(tmp, j)
        lengReach <-  lengReach + Length_RN[j]
        j_old <- j
        j <- DownNode_RN[j]} 
      if (lengReach > new_maxLength){
        j <- j_old
        whichNodeAG <- c(whichNodeAG[1], j, whichNodeAG[-1]) # new channel head is added at the front, so that same new_maxLength is applied
                                                             # second position (first is to be deleted)
        ChannelHeads[j] <- 1
        lengReach <- lengReach - Length_RN[j]
        tmp <- tmp[-length(tmp)]
      }
      Length_AG[reachID] <- lengReach
      RN_to_AG[tmp] <- reachID
      AG_to_RN[[reachID]] <- c(AG_to_RN[[reachID]], tmp)
      reachID <- reachID + 1
      whichNodeAG <- whichNodeAG[-1]
  }
  
  # # old version
  # if (equalizeLengths){
  #   while (length(whichNodeAG) != 0){ # explore all AG Nodes
  #     i <- whichNodeAG[1] # select the first
  #     RN_to_AG[i] <- reachID
  #     AG_to_RN[[reachID]] <- i
  #     j <- DownNode_RN[i]
  #     X_AG[reachID] <- X_RN[i]
  #     Y_AG[reachID] <- Y_RN[i]
  #     Z_AG[reachID] <- Z_RN[i]
  #     A_AG[reachID] <- A_RN[i]
  #     Length_AG[reachID] <- Length_RN[i]
  #     tmp_length <- Length_RN[i]
  #     tmp <- NULL
  #     j0 <- j
  #     while (!IsNodeAG[j] && j!=0) {
  #       tmp_length <-  tmp_length + Length_RN[j]
  #       j_old <- j
  #       j <- DownNode_RN[j]}
  # 
  #     if (tmp_length > maxReachLength){
  #       n_splits <- ceiling(tmp_length/maxReachLength)
  #       new_maxLength <- tmp_length/n_splits
  #       new_maxLength <- max(new_maxLength, 1.5*OCN$cellsize)
  #       j <- j0
  #       while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= new_maxLength) {
  #         tmp <- c(tmp, j)
  #         RN_to_AG[j] <- reachID
  #         AG_to_RN[[reachID]] <- c(AG_to_RN[[reachID]], tmp)
  #         Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
  #         j_old <- j
  #         j <- DownNode_RN[j]}
  #       if (Length_AG[reachID] > new_maxLength){
  #         j <- j_old
  #         Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
  #         ChannelHeads[j] <- 1
  #         whichNodeAG <- c(whichNodeAG,j)}
  # 
  #     } else {
  #       RN_to_AG[tmp] <- reachID
  #       AG_to_RN[[reachID]] <- c(AG_to_RN[[reachID]], tmp)
  #       Length_AG[reachID] <- tmp_length
  #     }
  # 
  #     reachID <- reachID + 1
  #     whichNodeAG <- whichNodeAG[-1]
  #   }
  # } else {
  #   while (length(whichNodeAG) != 0){ # explore all AG Nodes
  #     i <- whichNodeAG[1] # select the first
  #     RN_to_AG[i] <- reachID
  #     AG_to_RN[[reachID]] <- i
  #     j <- DownNode_RN[i] 
  #     X_AG[reachID] <- X_RN[i]
  #     Y_AG[reachID] <- Y_RN[i]
  #     Z_AG[reachID] <- Z_RN[i]
  #     A_AG[reachID] <- A_RN[i]
  #     #Length_AG[reachID] <- Length_RN[i]
  #     tmp_length <- Length_RN[i]
  #     tmp <- NULL
  #     j0 <- j # useless?
  #     while (!IsNodeAG[j] && j!=0 && tmp_length <= maxReachLength) {
  #       tmp <- c(tmp, j)
  #       tmp_length <-  tmp_length + Length_RN[j]
  #       j_old <- j
  #       j <- DownNode_RN[j]} 
  #     if (tmp_length > maxReachLength){
  #       j <- j_old
  #       whichNodeAG <- c(whichNodeAG, j)
  #       ChannelHeads[j] <- 1
  #       tmp_length <- tmp_length - Length_RN[j]
  #       tmp <- tmp[-length(tmp)]
  #     }
  #     Length_AG[reachID] <- tmp_length
  #     RN_to_AG[tmp] <- reachID
  #     AG_to_RN[[reachID]] <- c(AG_to_RN[[reachID]], tmp)
  #     reachID <- reachID + 1
  #     whichNodeAG <- whichNodeAG[-1]
  #   }
  # }

  Nnodes_AG <- length(X_AG)
  
  # FD_to_SC: vector of length OCN$FD$nNodes containing subcatchmentID for every pixel of the catchment
  # AG_to_FD: list containing FD indices of pixels belonging to a given reach 
  # SC_to_FD: list containing FD indices of pixels belonging to a given subcatchment 
  FD_to_SC <- numeric(OCN$FD$nNodes)
  
  # initialize FD_to_SC by attributing SC values to pixels belonging to AG level
  FD_to_SC[RN_mask] <- RN_to_AG 
  
  # attribute new SC values to pixels corresponding to outlets of catchments without reaches (because the drained area of the catchment is < thrA)
  Nnodes_SC <- Nnodes_AG + sum(OCN$FD$A[OCN$FD$outlet]<thrA)
  FD_to_SC[OCN$FD$outlet[OCN$FD$A[OCN$FD$outlet] < thrA]] <- (Nnodes_AG+1):Nnodes_SC 
  IndexHeadpixel <- which(OCN$FD$A==OCN$cellsize^2) # find FD pixels corresponding to headwaters
  
  
  AG_to_FD <- vector("list", Nnodes_AG)
  for(i in 1:Nnodes_AG) { # attribute river network pixels to fields of the AG_to_FD list 
    AG_to_FD[[i]] <- RN_to_FD[AG_to_RN[[i]]]
  }
  SC_to_FD <- AG_to_FD[1:Nnodes_AG] # initialize SC_to_FD by attributing the pixels that belong to reaches
  
  # add pixels corresponding to outlets of catchments without reaches
  if (Nnodes_SC > Nnodes_AG){
    for (i in (Nnodes_AG+1):Nnodes_SC){
      SC_to_FD[[i]] <- OCN$FD$outlet[OCN$FD$A[OCN$FD$outlet]<thrA][i-Nnodes_AG]
    }}
  
  # for (i in 1:length(IndexHeadpixel)){ # i: index that spans all headwater pixels
  #   p <- IndexHeadpixel[i] # p: ID of headwater pixel
  #   pNew <- p; # pNew: pixel downstream of p 
  #   k <- 0; # k: SC value of pixel pNew
  #   sub_p <- integer(0) # sub_p is the subset of pixels downstream of pixel p
  #   while (k==0){ # continue downstream movement until a pixel to which the SC has already been attributed is found
  #     k <- FD_to_SC[pNew]
  #     if (k==0){
  #       sub_p <- c(sub_p,pNew)
  #       pNew <- OCN$FD$downNode[pNew]
  #     }}
  #   FD_to_SC[sub_p] <- k
  #   SC_to_FD[[k]] <- c(SC_to_FD[[k]],sub_p)
  # }
ll <- continue_FD_SC(IndexHeadpixel, FD_to_SC, SC_to_FD, OCN$FD$downNode)
FD_to_SC <- ll$FD_to_SC
SC_to_FD <- ll$SC_to_FD
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # CALCULATE PROPERTIES AT AG LEVEL ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  #print('W matrix at AG level...',quote=FALSE); 
  # Adjacency matrix at reach level
  DownNode_AG <- numeric(Nnodes_AG)
  W_AG <- spam(0,Nnodes_AG,Nnodes_AG)
  ind <- matrix(0,Nnodes_AG,2)
  reachID <- sum(ChannelHeads) + 1
  for (i in 1:Nnodes_RN){ 
    if (DownNode_RN[i] != 0 && RN_to_AG[DownNode_RN[i]] != RN_to_AG[i]) {
      DownNode_AG[RN_to_AG[i]] <- RN_to_AG[DownNode_RN[i]]
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],DownNode_AG[RN_to_AG[i]])
    }
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  Outlet_AG <- RN_to_AG[Outlet_RN]
  
  # reverse downNode_AG
  DownNode_AG_rev <- vector("list",Nnodes_AG)
  for (i in 1:Nnodes_AG){
    d <- DownNode_AG[i]
    if (d!=0){DownNode_AG_rev[[d]] <- c(DownNode_AG_rev[[d]],i)  }}
  
  # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
  # sort nodes in downstream direction
  ind_sort <- sort(A_AG, index.return=TRUE)
  ind_sort <- ind_sort$ix
  Upstream_AG <- vector("list",Nnodes_AG)
  Nupstream_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    ups <- as.numeric(DownNode_AG_rev[[ind_sort[i]]])
    nodes <- numeric(0)
    for (u in ups){ nodes <- c(nodes, Upstream_AG[[u]])}
    Upstream_AG[[ind_sort[i]]] <- c(nodes, ind_sort[i])
    Nupstream_AG[ind_sort[i]] <- length(Upstream_AG[[ind_sort[i]]])
  }
 
  # AG_to_CM[i] indicates outlet to which reach i drains
  AG_to_CM <- numeric(Nnodes_AG)
  for (i in 1:OCN$nOutlet){
    AG_to_CM[Upstream_AG[[Outlet_AG[i]]]] <- i
  }
  
  ind_sort <- sort(A_AG, index.return=T)
  ind_sort <- ind_sort$ix
  if (streamOrderType=="Strahler"){
    # calculate Strahler stream order
    StreamOrder_AG <- numeric(Nnodes_AG)
    for (j in ind_sort){
      tmp <- DownNode_AG_rev[[j]] # set of reaches draining into j
      if (length(tmp)>0){
        IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
        if (IncreaseOrder > 1) {
          StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
        } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
      } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
    }
  } else if (streamOrderType=="Shreve"){
    # calculate Shreve stream order
    StreamOrder_AG <- numeric(Nnodes_AG)
    for (j in ind_sort){
      tmp <- DownNode_AG_rev[[j]]  # set of reaches draining into j
      if (length(tmp)>0){
        StreamOrder_AG[j] <- sum(StreamOrder_AG[tmp])
      } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
    } 
  }
  
  # Calculate slopes of reaches
  Slope_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    if (!(i %in% Outlet_AG))
      Slope_AG[i] <- (Z_AG[i] - Z_AG[DownNode_AG[i]])/Length_AG[i]
  }

  if(displayUpdates){message("Calculating network at AG level... 100.0%\n", appendLF = FALSE)}
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # CALCULATE PROPERTIES AT SC LEVEL ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if(displayUpdates){message("Calculating network at SC level...      \r", appendLF = FALSE)}
  
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE) 
  #t1 <- Sys.time()
  
  #print('Subcatchment properties...',quote=FALSE) 
  # calculate subcatchment properties: Local Elevation, Local Drained Area, Upstream Area
  Z_SC <- numeric(Nnodes_SC)
  Alocal_SC <- numeric(Nnodes_SC)
  for (i in 1:Nnodes_SC) {
    Z_SC[i] <- mean(OCN$FD$Z[SC_to_FD[[i]]])
    Alocal_SC[i] <- length(SC_to_FD[[i]])*OCN$cellsize^2
  }
  # drained area at AG level: note that the first Nnodes_AG elements of Alocal_SC correspond to subcatchments with reaches
  
  #  Areach_AG: includes the areas drained by the reaches
  Areach_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG) {
    Areach_AG[i] <- sum(Alocal_SC[Upstream_AG[[i]]])  
  }

  # coordinates of AG nodes considered at the downstream end of the respective edge
  XReach <- numeric(Nnodes_AG)
  YReach <- numeric(Nnodes_AG)
  ZReach <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    tmp <- AG_to_RN[[i]]
    ind <- which(A_RN[tmp]==max(A_RN[tmp]))
    node <- tmp[ind]
    XReach[i] <- X_RN[node]
    YReach[i] <- Y_RN[node]
    ZReach[i] <- Z_RN[node]
  }
  XReach[Outlet_AG] <- NaN
  YReach[Outlet_AG] <- NaN
  ZReach[Outlet_AG] <- NaN
  
  # build neighbouring nodes at FD level
  # find list of possible neighbouring pixels
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  if (length(OCN$typeInitialState)!=0 | OCN$FD$nNodes==OCN$dimX*OCN$dimY){ # all OCNs
    NeighbouringNodes <- NN_OCN(OCN$dimX, OCN$dimY, OCN$periodicBoundaries, movement)
    # NeighbouringNodes <- vector("list", OCN$dimX*OCN$dimY)
    # cont_node <- 0
    # for (cc in 1:OCN$dimX) {
    #   for (rr in 1:OCN$dimY) {
    #     cont_node <- cont_node + 1
    #     neigh_r <- rep(rr,8)+movement[1,]
    #     neigh_c <- rep(cc,8)+movement[2,]
    #     if (OCN$periodicBoundaries == TRUE){
    #       neigh_r[neigh_r==0] <- OCN$dimY
    #       neigh_c[neigh_c==0] <- OCN$dimX
    #       neigh_r[neigh_r>OCN$dimY] <- 1
    #       neigh_c[neigh_c>OCN$dimX] <- 1
    #     }
    #     NotAboundary <- neigh_r>0 & neigh_r<=OCN$dimY & neigh_c>0 & neigh_c<=OCN$dimX # only effective when periodicBoundaries=FALSE
    #     NeighbouringNodes[[cont_node]] <- neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*OCN$dimY
    #   }}
  } else {
    NeighbouringNodes <- NN_river(OCN$dimX, OCN$dimY, OCN$periodicBoundaries, movement, OCN$FD$toDEM, OCN$FD$nNodes)
    # NeighbouringNodes <- vector("list", OCN$dimX*OCN$dimY)
    # for (i in 1:OCN$FD$nNodes){
    #   nodeDEM <- OCN$FD$toDEM[i]
    #   cc <- (nodeDEM %% OCN$dimX); if (cc==0) cc <- OCN$dimX
    #   rr <- (nodeDEM - cc)/OCN$dimX + 1
    #   neigh_r <- rep(rr,8)+movement[1,]
    #   neigh_c <- rep(cc,8)+movement[2,]
    #   if (OCN$periodicBoundaries == TRUE){
    #     neigh_r[neigh_r==0] <- OCN$dimY
    #     neigh_c[neigh_c==0] <- OCN$dimX
    #     neigh_r[neigh_r>OCN$dimY] <- 1
    #     neigh_c[neigh_c>OCN$dimX] <- 1
    #   }
    #   NotAboundary <- neigh_r>0 & neigh_r<=OCN$dimY & neigh_c>0 & neigh_c<=OCN$dimX # only effective when periodicBoundaries=FALSE
    #   NeighbouringNodes[[nodeDEM]] <- (neigh_r[NotAboundary]-1)*OCN$dimX + neigh_c[NotAboundary]
    # }
  }

  if (OCN$FD$nNodes < OCN$dimX*OCN$dimY){ # general contour OCNs and real rivers
    NeighbouringNodes <- NN_FD(OCN$FD$nNodes, OCN$dimX, OCN$dimY, NeighbouringNodes, OCN$FD$toDEM)
    # NeighbouringNodes_FD <- vector("list", OCN$FD$nNodes)
    # DEM_to_FD <- numeric(OCN$dimX*OCN$dimY)
    # DEM_to_FD[OCN$FD$toDEM] <- 1:OCN$FD$nNodes
    # for (i in 1:OCN$FD$nNodes){
    #   indDEM <- OCN$FD$toDEM[i]
    #   tmp <- DEM_to_FD[NeighbouringNodes[[indDEM]]]
    #   NeighbouringNodes_FD[[i]] <- tmp[tmp != 0]
    # }
    # NeighbouringNodes <- NeighbouringNodes_FD
  }

  
  # Subcatchment adjacency matrix: find which subcatchments have borders in common
  # W_SC <- spam(0,Nnodes_SC,Nnodes_SC)
  # indices <- matrix(0,Nnodes_SC*20,2)
  # k <- 1
  # for (i in 1:Nnodes_SC){
  #   set <- SC_to_FD[[i]]
  #   nodes <- numeric(0)
  #   for (s in set){ nodes <- union(nodes, FD_to_SC[NeighbouringNodes[[s]]])}
  #   NeighSubcatch <- setdiff(nodes, i)
  #   indices[k:(k+length(NeighSubcatch)-1),1] <- i
  #   indices[k:(k+length(NeighSubcatch)-1),2] <- NeighSubcatch
  #   k <- k + length(NeighSubcatch)
  #   if (displayUpdates){
  #     if ((i %% max(1,round(Nnodes_SC*0.01)))==0){
  #       message(sprintf("Calculating network at SC level... %.1f%%\r",i/Nnodes_AG*100), appendLF = FALSE)}}
  # }
  # indices <- indices[1:(k-1),]
  # W_SC[indices] <- 1  
  
  ll <- WSC(Nnodes_SC,SC_to_FD,FD_to_SC,NeighbouringNodes)
  W_SC <- spam(0,Nnodes_SC,Nnodes_SC) 
  W_SC[cbind(ll[[1]],ll[[2]])] <- 1
  
  # X,Y of subcatchment centroids
  X_SC <- numeric(Nnodes_SC)
  Y_SC <- numeric(Nnodes_SC)
  for (i in 1:Nnodes_SC){
    X_SC[i] <- mean(OCN$FD$X[SC_to_FD[[i]]])
    Y_SC[i] <- mean(OCN$FD$Y[SC_to_FD[[i]]])
  }
  
  if(displayUpdates){message("Calculating network at SC level... 100.0%\n", appendLF = FALSE)}
  
  #%%%%%%%%%%%%%%%%%%%%%#
  # EXPORT VARIABLES ####
  #%%%%%%%%%%%%%%%%%%%%%#
  
  #FD level
  OCN$FD[["toRN"]] <- FD_to_RN
  OCN$FD[["toSC"]] <- FD_to_SC
  
  # RN level
  OCN$RN[["A"]] <- A_RN
  OCN$RN[["W"]] <- W_RN
  OCN$RN[["downNode"]] <- DownNode_RN
  OCN$RN[["drainageDensity"]] <- DrainageDensity_RN
  OCN$RN[["leng"]] <- Length_RN
  OCN$RN[["nNodes"]] <- Nnodes_RN
  OCN$RN[["nUpstream"]] <- Nupstream_RN
  OCN$RN[["outlet"]] <- Outlet_RN
  OCN$RN[["slope"]] <- Slope_RN
  OCN$RN[["toFD"]] <- RN_to_FD
  OCN$RN[["toAGReach"]] <- RN_to_AG
  OCN$RN[["toCM"]] <- RN_to_CM
  OCN$RN[["upstream"]] <- Upstream_RN
  OCN$RN[["X"]] <- X_RN
  OCN$RN[["Y"]] <- Y_RN
  OCN$RN[["Z"]] <- Z_RN
  
  # AG level
  OCN$AG[["A"]] <- A_AG
  OCN$AG[["AReach"]] <- Areach_AG
  OCN$AG[["W"]] <- W_AG
  OCN$AG[["downNode"]] <- DownNode_AG
  OCN$AG[["leng"]] <- Length_AG
  OCN$AG[["nNodes"]] <- Nnodes_AG
  OCN$AG[["nUpstream"]] <- Nupstream_AG
  OCN$AG[["outlet"]] <- Outlet_AG
  OCN$AG[["slope"]] <- Slope_AG
  OCN$AG[["streamOrder"]] <- StreamOrder_AG
  OCN$AG[["ReachToFD"]] <- AG_to_FD
  OCN$AG[["ReachToRN"]] <- AG_to_RN
  OCN$AG[["toCM"]] <- AG_to_CM
  OCN$AG[["upstream"]] <- Upstream_AG
  OCN$AG[["X"]] <- X_AG
  OCN$AG[["XReach"]] <- XReach
  OCN$AG[["Y"]] <- Y_AG
  OCN$AG[["YReach"]] <- YReach
  OCN$AG[["Z"]] <- Z_AG
  OCN$AG[["ZReach"]] <- ZReach
  
  # SC level
  OCN$SC[["ALocal"]] <- Alocal_SC
  OCN$SC[["W"]] <- W_SC
  OCN$SC[["nNodes"]] <- Nnodes_SC
  OCN$SC[["toFD"]] <- SC_to_FD
  OCN$SC[["X"]] <- X_SC
  OCN$SC[["Y"]] <- Y_SC  
  OCN$SC[["Z"]] <- Z_SC
  
  # other
  OCN$thrA <- thrA
  OCN$streamOrderType <- streamOrderType
  OCN$maxReachLength <- maxReachLength
  
  # re-define AG_to_RN, AG_to_FD, RN_to_AG considering AG nodes as pixels and not reaches
  AG_to_FDnode <- numeric(Nnodes_AG)
  AG_to_RNnode <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    tmpFD <- AG_to_FD[[i]]
    AG_to_FDnode[i] <- tmpFD[OCN$FD$A[tmpFD]==min(OCN$FD$A[tmpFD])]
    tmpRN <- AG_to_RN[[i]]
    AG_to_RNnode[i] <- tmpRN[OCN$RN$A[tmpRN]==min(OCN$RN$A[tmpRN])]
  }
  RN_to_AGnode <- numeric(Nnodes_RN)
  for (i in 1:Nnodes_AG){
    RN_to_AGnode[AG_to_RNnode[i]] <- i
  }
 
  OCN$RN[["toAG"]] <- RN_to_AGnode
  OCN$AG[["toFD"]] <- AG_to_FDnode
  OCN$AG[["toRN"]] <- AG_to_RNnode
  
  invisible(OCN)
}