newgac <- function(amat, x, y, z, type = "pag") 
{
  if (type %in% c("dag", "cpdag", "pdag")) {
    if (!isValidGraph(amat = amat, type = type)) {
      message("The input graph is not a valid ", type, 
              ". See function isValidGraph() for details.\\n")
    }
  }
  res <- rep(NA, 3)
  f <- NULL
  if (type %in% c("dag", "pdag", "cpdag")) {
    res[1] <- pcalg:::isAmenable(m = amat, x = x, y = y, type = type)
    f <- pcalg:::bforbiddenNodes(m = amat, x = x, y = y)
    res[2] <- (length(intersect(f, z)) == 0)
    res[3] <- newcond3fast(x = x, y = y, z = z, m = amat)
  }
  else {
    res[1] <- pcalg:::isAmenable(amat, x = x, y = y, type = type)
    f <- pcalg:::forbiddenNodes(amat, x = x, y = y)
    res[2] <- (length(intersect(f, z)) == 0)
    res[3] <- pcalg:::cond2(x = x, y = y, z = z, m = amat, type = type)
  }
  list(gac = all(res), res = res, f = f)
}

newcond3fast <- function(x = x, y = y, z = z, m = m)
{
  oneDag <- pdag2dag(as(t(m), "graphNEL"))
  dagAmat <- t(as(oneDag$graph, "matrix"))
  gb <- newgbg(dagAmat, x, y)
  msep(t(gb), x, y, z)
}


##needed to delete i <- 1 from gbg
newgbg <- function(m, x, y) 
{
  tmp <- m
  for (i in 1:length(x)) {
    Desc <- pcalg:::bPossibleDeProper(m, x[i], x[-i])
    if (length(intersect(y, Desc)) != 0) {
      ch <- as.vector(which(m[x[i], ] == 0 & m[, x[i]] == 
                              1))
      cand <- intersect(ch, Desc)
      j <- 0
      while (j < length(cand)) {
        j <- j + 1
        pathOK <- ((length(intersect(y, pcalg:::bPossibleDeProper(m, 
                                                          cand[j], x))) != 0) | (cand[j] %in% y))
        if (pathOK) {
          tmp[cand[j], x[i]] <- 0
        }
      }
    }
  }
  tmp
}

