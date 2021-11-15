# Topological data analysis of a cloud of points in N dimensions. Takes the matrix of
# the point coordinates with the coordinates in the columns and the points in the rows and
# an optional upper limit of the number of cliques in the graph constructed by connecting
# points at a distance bellow a threshold. The threshold is determined automatically.


stda = function(m, maxCliqN = 1000, dv = 2) {
  library(igraph)
  nervD = dist(m)
  dr = range(nervD)
  nervD = as.matrix(nervD)
  thr = dr[1] + diff(dr) / dv
  repeat {
    G = dTG(nervD, thr)
    Cli = cliques(G)
    l = length(Cli)
    if (l>10000) thr=thr/2 else {
    if (abs(l - maxCliqN) > (0.1 * maxCliqN)) {
      thr = thr - thr * 0.01 * (l - maxCliqN) / l
    }
    else
      break
    print(paste("threshold", thr, "No of cliques", l, collapse = " "))
  }}
  print(paste(c("Threshold used -", thr, "from range", dr), collapse = " "))
  D = ChC(m, Cli)
  Cqe = D[[1]]
  Cqi = D[[2]]
  D = D[[3]]
  print("Boundary matrix completed.")
  Drs = EH(D, Cqi)
  Dr = Drs[[1]]
  Betti = Drs[[3]]
  res = img(Dr, Cqe, Cqi, t = paste("Cut off", thr, collapse = " "))
  j = res[[1]]
  lres = res[[2]]
  cycles = res[[3]]
  peres = cycles[cycles[, 2] == 0, ]
  colnames(peres) = c("Birth Radius", "End Radius")
  # return names of persistent cycles
  return(list(Persistent_Cycles = peres, Betti = Betti))
}

dTG = function(dm, coff, pl = T) {
  require(igraph)
  require(reshape2)

  # distance matrix to graph
  # for each point take only the nearest neighbors with distances to that point less than coff
  AL = sapply(seq_along(dm[, 1]), function(i) {
    x = dm[i, -i]
    x[x < coff]
  })
  names(AL) = rownames(dm)
  L0 = lapply(AL, as.list)

  # remove singlets
  L = L0[lengths(L0) > 0]
  if (length(L) < 4) {
    print("Less than 4 verteces")
    return(NULL)
  }

  # transform to edge list - only couples of points which are connected
  EL = melt(L)
  EL = EL[, c(3, 2, 1)]

  rownames(EL) = seq_along(EL[, 1])
  colnames(EL) = c("From", "To", "Weight")
  EL = as.data.frame(EL, stringsAsFactors = F)
  EL[, 3] = as.double(EL[, 3])
  EL = EL[EL[, 3] != 0, ]

  # covert distance to weight
  EL[, 3] = 1 / (EL[, 3]) ^ 2

  G = graph_from_data_frame(EL, directed = F)
  G=simplify(G)
  if (pl)
    plot(G,
         vertex.size = 1,
         main = paste("Graph at Cut Off", coff, collapse = " "))
  return(G)
}

ChC = function(m, cc) {
  require(parallel)
  require(Rfast)

  # Builds a Cech complex from a n-dim points matrix and a clique list of a graph from
  # the distance matrix of those points (igraph) at a given cut off.
  # The names of the vertices are numbers. The clique list has to be calculated at high epsilon
  # of the coverage so that all cliques of interest are included and the algorithm can work \
  # backwards reconstituting the Cech complexes at all lower epsilon values.
  # Outputs D[i,j] general matrix as described
  # in Edelsbrunner and Harer, "Persistance Homology - a survey"

  mxm = max(lengths(cc))
  L = length(cc)
  Cqi = lapply(cc,function(cl) as.numeric(names(cl)))
  k = ncol(m) + 2

  # A separate list with a structure parallel to that of Cqi which holds the epsilons
  # (the diameter of the spheres) at which each clique turns into a simplex

  cl <- makeCluster(4)
  ex <-
    Filter(function(x)
      is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  clusterExport(cl, ex)
  clusterEvalQ(cl, library(pdist))
  clusterExport(cl, list("Cqi", "m"), envir = environment())

  Cqe = parSapply(cl, Cqi, function(l) {
    n = length(l)
    if (n == 1)
      res = 0
    if (n == 2)
      res = dist(m[l, ])
    if (n > 2 & n < k) {
      mx = m[l, ]
      res = 2 * CHM(mx)
    }
    if (n >= k)
      res = -Inf
    return(res)
  })
  stopCluster(cl)
  clovr = which(is.infinite(Cqe))
  clovR = sapply(clovr, function(i) {
    faces = which(sapply(Cqi, function(c) {
      all(c %in% Cqi[[i]]) & length(c) == (length(Cqi[[i]]) - 1)
    }))
    #print(Cqe[faces])
    return(max(Cqe[faces]))
  })
  if (length(clovr) > 0)
    Cqe[clovr] = clovR
  print(sum(is.infinite(Cqe)))
  print("Epsilons calculated...")
  Cqi = Cqi[order(Cqe)]
  Cqe = Cqe[order(Cqe)]
  #Cqi=Cqi[order(lengths(Cqi),Cqe)]

  x = Cqi[lengths(Cqi) == 2]
  y = rank(Cqe[lengths(Cqi) == 2])
  x = rbind(
    data.frame(
      Name = sapply(x, function(c) {
        c[1]
      }),
      Rank = y,
      stringsAsFactors = F
    ),
    data.frame(
      Name = sapply(x, function(c) {
        c[2]
      }),
      Rank = y,
      stringsAsFactors = F
    )
  )
  x = aggregate(Rank ~ Name, data = x, FUN = min)
  x = x[order(x[, "Rank"]), "Name"]
  y = sapply(Cqi[lengths(Cqi) == 1], function(l) {
    l %in% x
  })
  Cqi[lengths(Cqi) == 1][y] = lapply(1:sum(lengths(Cqi) == 1), function(i) {
    Cqi[[i]] = x[i]
  })

  # Form the combined boundary matrix
  print("Calculating boundary matrix...")
  D = sapply(Cqi, function(sj) {
    sapply(Cqi, function(si) {
      if (sum(!(sj %in% si)) == 1)
        res = all(si %in% sj)
      else
        res = 0
      return(res)
    })
  })
  s = length(Cqi[lengths(Cqi) == 1])
  D = rbind(rep(0, ncol(D)), D)
  D = cbind(rep(0, nrow(D)), D)
  Cqe = c(0, Cqe)            # adding a dummy first cycle of NOTHING that dies with the first vertex
  D[1, 2:(s + 1)] = 1
  D = D * 1
  colnames(D) = c(0, sapply(Cqi, paste, collapse = "|"))
  rownames(D) = c(0, sapply(Cqi, paste, collapse = "|"))

  return(list(Cqe, Cqi, D))
}

EH = function(D, Cqi) {
  # Based on Edelbrunner&Harer Computational Topology, takes joint boudary matrix
  # and the list of cliques and returns the reduced boundary matrix, the number of cycles
  # the Betti numbers

  Dr = D
  print(c("Boundary matrix dimensions ", dim(Dr)))
  cndr = colnames(Dr)[1]
  lows = apply(Dr, 2, function(c) {
    if (max(c) == 1)
      max(which(c == 1))
    else
      0
  })
  for (j in 2:ncol(Dr)) {
    cndr[j] = colnames(Dr)[j]
    loL = which(lows[1:(j - 1)] == lows[j] & lows[1:(j - 1)] != 0)
    while (length(loL) > 0) {
      j1 = min(loL)
      Dr[, j] = Dr[, j] + Dr[, j1]
      Dr[, j] = Dr[, j] %% 2
      cndr[j] = list(unlist(cndr[c(j1, j)]))
      lows = apply(Dr, 2, function(c) {
        if (max(c) == 1)
          max(which(c == 1))
        else
          0
      })
      loL = which(lows[1:(j - 1)] == lows[j] & lows[1:(j - 1)] != 0)
    }
  }
  print("RBM ready.")
  cndr = sapply(cndr, function(n) {
    x = table(n)
    x = names(x)[x %% 2 == 1]
    x = paste(x, collapse = "+")
    return(x)
  })
  colnames(Dr) = cndr
  psimp = lengths(Cqi) - 1
  Zc = colSums(Dr) == 0
  Zc = Zc[-1]

  ZLp = t(sapply(unique(psimp), function(p) {
    z = sum(Zc[psimp == p])
    l = sum(psimp == p) - z
    return(c(z, l))
  }))
  print(ZLp)
  betti = c(ZLp[1:(nrow(ZLp) - 1), 1] - ZLp[2:nrow(ZLp), 2], ZLp[nrow(ZLp), 2])

  return(list(Dr, ZLp, betti))
}

# Imaging of reduced boundary matrix and the persistent homology plot
# takes the reduced boundary matrix, list of epsilon (radius) corresponding
# to the dissapearance of cycles and emergence of simplices
# as well as the list of cliques

img = function(m, eps, clq, t = NULL) {
  lows = apply(m, 2, function(c) {
    if (max(c) == 1)
      return(max(which(c == 1)))
    else return(0)
  })

  d = c(0, lengths(clq))
  lowz = which(lows == 0)
  lows = cbind(lows, 1:ncol(m))   #generate the reduced boundary matrix enhanced image
  lowss = lows[lows[, 1] > 0, ]
  m[lowss] = m[lowss] + 1

  for (i in lowz) {
    if (i %in% lowss[, 1]) {
      m[i, i:(which(lows[, 1] == i) - 1)] = 3
    }
    else
      m[i, i:ncol(m)] = 3
  }
  cols = max(m) + 1
  mt = t(m[nrow(m):1, ])
  bw = colorRampPalette(c("white", "grey", "black", "red"))
  if (ncol(m) > 30) {
    rownames(mt) = NULL
    colnames(mt) = NULL
  }
  print(c("max eps", max(eps)))
  image(mt,axes = FALSE, xlab = "", ylab = "", col = bw(cols), main = t)
  axis(2, at = ((ncol(mt):1) - 1) / (ncol(mt) - 1), labels = F)
  text(par("usr")[1] - 0.02,  ((ncol(mt):1) - 1) / (ncol(mt) - 1),  cex = 0.5,  srt = 45, labels = colnames(mt)[ncol(mt):1], pos = 2, xpd = NA)
  axis(3, at = ((nrow(mt):1) - 1) / (nrow(mt) - 1), labels = F)
  text(((nrow(mt):1) - 1) / (nrow(mt) - 1), par("usr")[4] + 0.03, cex = 0.5, srt = 45, labels = rownames(mt)[nrow(mt):1], pos = 3, xpd = NA )
  box(which = "plot")
  if (ncol(m) < 30) grid(nrow(mt), ncol(mt))

  # generate the persistent homology plot
  eps[eps == 0] = eps[eps == 0] + runif(length(eps[eps == 0]), 0, min(eps[eps > 0] / 5)) # stagger the zerro cycle lines to make them all visible
  x = eps[unlist(rep(lowz[-1], each = 2))]
  ye = sapply(lowz[-1], function(i) {
                if (2 %in% m[i, ])
                  return(which(m[i, ] == 2))
                else  return(0)
      })
  ye[ye == 0] = ncol(m)
  y = eps[unlist(c(rbind(lowz[-1], ye)))]
  col = d[lowz[-1]]
  z=1:(length(lowz) - 1)
  z=z[order(col)]
  for (i in z) {
    j = i * 2 - 1
    plot(x[j:(j + 1)], y[j:(j + 1)], xlim = c(0, max(x)), ylim = c(0, max(y)), cex = 0, xlab = "", ylab = "")
    lines(x[j:(j + 1)], y[j:(j + 1)], xlim = c(0, max(x)), ylim = c(0, max(y)), col = col[i])
    par(new = T)
  }

  xa = x[(1:(length(x) / 2)) * 2]
  ya = eps[ye]
  cx = 1 * (ye != ncol(m))
  plot(xa[z], ya[z], xlim = c(0, max(x)), ylim = c(0, max(y)), cex = cx[z], col = col[z], xlab = "Epsilon at Emergence of Cycle", ylab = "Epsilon at Disappearance of Cycle", main = t, pch=16)

  for (i in z) {
    j = i * 2 - 1
    plot(y[j:(j + 1)], c(i,i), xlim = c(0, max(y)), ylim = c(1,max(z)), cex = 0, xlab = "", ylab = "")
    lines(y[j:(j + 1)], c(i,i), xlim = c(0, max(y)), ylim = c(1,max(z)), col = col[i])
    par(new = T)
  }
  plot(ya[z], z, xlim = c(0, max(y)), ylim = c(1,max(z)), cex = cx[z]*0.7, col = col[z], xlab = "Epsilon", ylab="Cycles", main = t, pch=16)

  return(list(lowz[-1], data.frame(Simplices = rownames(m), Cycles = colnames(m), stringsAsFactors = F), data.frame(Beginning = xa, End = ya * cx, row.names = colnames(m)[lowz][-1], stringsAsFactors = F)))
}

CHM = function(m) {
  # Algorithm for finding the circumscribed ball or smallest enclosing m-ball
  # for m-simplices following Cheng/Hu/Martin's approach by projecting the n-dim points from n to
  # the appropriate m dimensions and the radius back to n by prcomp, or taking the radius of the last
  # n-dimensional simplex which is a subset of the m-dimensional simplex viewed when m>n
  # m is the kxn matrix of points.

  require(pdist)
  k = nrow(m)
  n = ncol(m)
  pc = prcomp(m)
  mk = pc$x[, 1:(k - 1)]
  A = t(sapply(2:k, function(i) {
    mk[i, ] - mk[i - 1, ]
  }))
  B = rowSums(t(sapply(2:k, function(i) {
    mk[i, ] ^ 2 - mk[i - 1, ] ^ 2
  }))) / 2

  C = solve(A) %*% B
  c = t(C) %*% t(pc$rotation[, 1:(k - 1)]) + colMeans(m)
  d = as.numeric(pdist(c, m)@dist)
  dm = mean(d)
  return(dm)
}
