# This code implements vintage sparse PCA (vsp)
# allows for symmetric or directed or bipartite (i.e. rectangular A).
# you may input a matrix, or a Matrix, or an igraph.


# TODO:
# - Currently has adjMat as input.  Should could allow igraph object and tidy text.
# - automatic diagnostics for k selection
# - bind u and v for first rotation, then rotate u and v.  then do u and v individually... hopefully, this makes B more diagonal.
# - allow for random restarts of varimax.

# flow:
#  1) compute vsp with k as large as you think it could be.
#  2) inspect eigengap / other diagnostics.
#  3) after selecting desired k, input the output of step 1 back into vsp + selected k.


# Version 1.0:  January 26, 2020  karlrohe@stat.wisc.edu

# library(irlba)
library(rARPACK)
library(Matrix)
library(tidyverse)


vsp = function(data, k = 5, tau = -1, normalization = T, centering = T, recenter = T, rescale = T, sym = NULL, fast = TRUE){
  # data is a matrix or a Matrix or an igraph
  # k is the number of clusters.
  # tau is regularization parameter.  it can be a 2 vectors (one for row, another for column).  The default is pretty good
  # normalization = T divides element A_{i,j} by sqrt (row sum * col sum); tau regularization is used in that step.
  # centering =T does row&column centering.  it is performed implicitly in the SVD.
  # For the SVD, vsp uses rARPACK with a custom vector-matrix multiply function (fcent defined below)
  # if A is square & rowSums(A) = colSums(A) & A is not symmetric, set sym = FALSE (to save time)
  # fast ==T only uses largest 10% of leverage scores in computing varimax rotation.


  ######
  # First, compute the Principal Components:
  ######

  # if data is an igraph, do PCA on adjacency matrix
  if(is.igraph(data)) pcs= bigPCA(A = as_adjacency_matrix(data), k, tau, normalization, centering,sym)


  # if data is some sort of matrix or Matrix, then do PCA.
  if(!(data %>% dim %>% is.null() )) pcs= bigPCA(A = data, k, tau, normalization, centering,sym)
  # data= bigPCA(data,k = k, tau = tau, normalization = normalization, centering = centering, sym = sym, tw = tw)


  ######
  # Second, compute the varimax rotation. the function also recenters and rescales (rescaling not yet implemented)
  ######

  vsp_factors = vsp_rotate(pcs, k, fast, recenter, rescale)

  return(vsp_factors)
  # print("input 'data' should be either a matrix or the output of a previous run of vsp. input not read as one of these forms. no output. ")
}
bigPCA = function(A, k = 5, tau = -1, normalization = T, centering = T, sym = NULL){
  # some defaults if needed for testing: k = 5, tau = -1, normalization = T, centering = T, recenter = T, rescale = T, sym = NULL

  rs = rowSums(A)
  cs = colSums(A)

  n = nrow(A)
  d = ncol(A)

  if(n!=d) sym = FALSE
  if(is.null(sym)) if(sum(!(rs!=cs))>1) sym = FALSE
  if(is.null(sym)) sym = isSymmetric(A)



  # the next lines define matrix L.
  # in memory constrained settings, this could be put into
  #   matrix-vector multiplication function fcent.

  # if A contains ANY negative elements, turn off normaliation and print that warning.
  if(normalization){
    if(is(A, 'sparseMatrix')) if(min(A@x) < 0){
      print("Normalization turned off because data matrix has negative elements. Normalization should not be used in such cases.")
      normalization=F
    }
    if(is.matrix(A)) if(min(A) < 0){
      print("Normalization turned off because data matrix has negative elements. Normalization should not be used in such cases.")
      normalization=F
    }
  }

  if(normalization){

    if(tau <0)  tau = mean(rs)
    D = Diagonal(n = n, x = 1/sqrt(rs + tau))
    tmp = D%*%A
    D = Diagonal(n = d, x = 1/sqrt(cs + mean(cs)))
    L = tmp %*% D

  }

  if(!normalization) L = A



  # this is a fancy way of doing the eigendecomposition for the centered
  #  version of L... that is, the matrix that is like L - 1 d^T - d 1^T + 11^T
  #  the reason it is fancy is that multiplying by that matrix can be fast, but can't
  #  define the full nxn (dense) matrix...
  if(!centering){
    rsL = rep(0,nrow(L))
    csL = rep(0, ncol(L))
    grandMean = 0
  }
  if(centering){
    rsL = rowSums(L)
    csL = colSums(L)
    grandMean = mean(L)
  }
  args = list(A = L,At=t(L), rs = rsL, cs = csL, n = nrow(L), d = ncol(L),
              meanA = grandMean, onesR = rep(1, nrow(L)), onesC = rep(1, ncol(L)))
  ei = eigs_sym(fcent, n=nrow(L), k=k, args=args, which = "LA")
  U = ei$vectors
  # compute localization statistic
  loc_U_stat = colSums(U^4)
  loc_V_stat = NULL

  V = NULL
  if(!sym){
    ss = t(U)%*%L %>% t %>% as.matrix %>% svd
    V = ss$u
    loc_V_stat = colSums(V^4)
  }

  scree =sqrt(ei$values)



  out = list(U = U, V= V, scree = scree, sym = sym, centered = centering, normalized = normalization, loc_U_stat=loc_U_stat, loc_V_stat=loc_V_stat, rsA = rs, csA = cs, rsL = rsL, csL = csL, bigPCAoutput = T)


  class(out) = "bigpc"
  return(out)
}

fcent = function(x, args){
  # multiplies centered(args$A)%*%x, quickly for sparse A.

  # args$ contains: A, rs, n, meanA, onesR, At
  # if problem is symmetric is_null(At) is true
  # if problem is not symmetric, then must also contain cs, d, onesC
  #  A is dgcMatrix
  #  rs is rowSums(A)
  #  n is nrow(A)
  #  meanA = mean(A)
  #  onesR = rep(1, n) # R stands for row (length n)
  #  At = t(A) (or At = NULL)
  #  cs = colSums(A)
  #  d = ncol(A)
  #  onesC = rep(1,d)  # C stands for column (length d)

  # this code stores both A and At because we presume that At%*%x is faster than t(A)%*%x.
  #   this is the case when A is of the type dgCMatrix and At is also dgCMatrix.

  # here is test code for Ax and Atx which do left and right multiplication.
  # A = matrix(rexp(1000*100, 10), nrow = 1000)
  # n = nrow(A)
  # d = ncol(A)
  # rs = rowSums(A)
  # cs = colSums(A)
  # meanA = mean(A)
  #
  # Acent = A - outer(rep(1/n,n), cs) - outer(rs, rep(1/d,d)) + outer(rep(1,n), rep(1,d))*meanA
  #
  #
  # x = rnorm(d)
  # args = list(A = A, rs = rowSums(A), cs = colSums(A), n = nrow(A), d = ncol(A), meanA = mean(A), onesR = rep(1, nrow(A)), onesC = rep(1, ncol(A)))
  # sd(Ax(x,args) - (Acent%*%x))
  #
  # x = rnorm(1000)
  # args = list(A = A,At=t(A), rs = rowSums(A), cs = colSums(A), n = nrow(A), d = ncol(A), meanA = mean(A), onesR = rep(1, nrow(A)), onesC = rep(1, ncol(A)))
  # sd(Atx(x,args) - (t(Acent)%*%x))

  if(is_null(args$At)){ # if tA is not stored, then the problem is symmetric.
    # mx = mean(x)
    return(
      Ax(x, args)
      # as.vector(args$A %*%x - args$onesR*(as.numeric((t(args$rs)%*%x)/args$n) - mx*args$meanA)  - args$rs*mx )
    )
  }
  # otherwise, the problem is asymmetric.
  return(Ax(Atx(x, args), args))
}

Ax = function(x, args){
  # mx = mean(x)
  return(
    as.vector(args$A %*%x - args$onesR*(as.numeric((t(args$cs)%*%x)/args$n) - sum(x)*args$meanA)  - args$rs*mean(x) )
  )
}

Atx = function(x, args){
  # mx = mean(x)
  return(
    as.vector(args$At %*%x - args$onesC*(as.numeric((t(args$rs)%*%x)/args$d) - sum(x)*args$meanA)  - args$cs*mean(x) )
  )
}

vsp_keep = function(vsp_factors, keepThese, fast = T, recenter = NULL, rescale = NULL){
  # #  this function is written such that you can do the following:

  # # Are any of the singular vectors localized on a few documents?
  # #   Lets find them. remove them.  find an eigengap.  then rotate.
  # #   vsp_keep makes sure that we don't have to do another svd.

  # U = v_cent$pcs$U
  # l4 = colSums(U^4)^(1/4)  # this is a measure of localization.
  # localized = l4>.15  # .15 was a manually chosen cutoff.
  # hist(U[,localized], breaks = 1000)  # notice how the localized eigenvectors are highly highly skewed.
  # hist(U[,!localized], breaks = 1000  # while the others are only moderately skewed
  #
  # # remove the singular values that correspond to localized
  # #   singular vectors.  then, make a screeplot

  # v_cent$pcs$scree[which(!localized)] %>% plot
  # # there is an eigengap at 8
  # # so, these are the singular vectors we want to keep and rotate:
  # good = which(!localized)[1:8]

  # #  you can pass this to vsp_keep, with the old vsp object,
  # system.time({ vsp_factors = vsp_keep(v_cent, keepThese = good, fast = F, recenter=F)})


  if(length(keepThese) == 1) keepThese = 1:keepThese
  vsp_factors$pcs$U = vsp_factors$pcs$U[,keepThese]
  vsp_factors$pcs$V = vsp_factors$pcs$V[,keepThese]
  vsp_factors$pcs$scree =  vsp_factors$pcs$scree[keepThese]
  vsp_factors$pcs$loc_U_stat = vsp_factors$pcs$loc_U_stat[keepThese]
  vsp_factors$pcs$loc_V_stat = vsp_factors$pcs$loc_V_stat[keepThese]

  if(is.null(recenter)) recenter = vsp_factors$recentered
  if(is.null(rescale)) rescale = vsp_factors$rescaled

  return(vsp_rotate(vsp_factors$pcs, k  = length(keepThese), fast, recenter = recenter, rescale = rescale))
}

vsp_rotate = function(rotate_input, k, fast = T, recenter = NULL, rescale = NULL){

  if(class(rotate_input) == "vsp"){
    if(is.null(recenter)) recenter = rotate_input$recentered
    if(is.null(rescale)) rescale = rotate_input$rescaled
    pcs = rotate_input$pcs
  }
  if(class(rotate_input) == "bigpc") pcs = rotate_input

  #
  # k can be numeric (top k) or a vector of logicals or a vector of indices

  # some defaults if needed for testing:
  # k = 5, fast = T, recenter = NULL, rescale = NULL


  if(length(k) == 1){
    U = pcs$U[,1:k]
    V = pcs$V[,1:k]
    scree =  pcs$scree[1:k]
  }
  if(length(k) > 1){
    U = pcs$U[,k]
    V = pcs$V[,k]
    scree = pcs$scree[k]
  }

  if(fast & nrow(U) < 1000){
    print("fast rotation has been turned off because standard varimax will be fast already")
    fast = F
  }

  Rz = eiv(U, fast = fast)  # compute varimax rotation and ensure each factor will skew + after rotation.
  Z = U%*%Rz*sqrt(nrow(U))



  # if symmetric analysis, then compute B with Rz's and return
  muZ = NULL
  muY= NULL

  if(pcs$sym | is.null(V)){
    pcs$sym = T
    Y = NULL
    Ry = NULL

    B = t(Rz)%*%diag(scree[1:k])%*%Rz
    rho = mean(B)
    B = B/rho
    if(recenter){
      # recenter Z
      # note that this uses pcs$U because this is a symmetric decomposition:
      muZ = vsp_recenter(means = pcs$csA/nrow(Z), U, scree, Rz)/sqrt(nrow(Z))
      Z = Z + matrix(muZ, nrow = nrow(Z), ncol = ncol(Z), byrow=T)
    }
    if(rescale){}  #TODO. also change output to rescaled = rescale
  }

  # if asymmetric analysis, then compute Ry, Y, and then B.
  if(!pcs$sym){
    Ry = eiv(V, fast = fast)
    B = t(Rz)%*%diag(scree)%*%Ry
    rho = mean(B)
    B = B/rho
    Y = V%*%Ry*sqrt(nrow(V))
    if(recenter){
      muZ = vsp_recenter(means = pcs$csA/nrow(Z), V, scree, Rz)/sqrt(nrow(Z))
      muY = vsp_recenter(means = pcs$rsA/nrow(Y), U, scree, Ry)/sqrt(nrow(Y))
      Z = Z + matrix(muZ, nrow = nrow(Z), ncol = ncol(Z), byrow=T)
      Y = Y + matrix(muY, nrow = nrow(Y), ncol = ncol(Y), byrow=T)
    }
    if(rescale){}  #TODO. also change output to rescaled = rescale
  }
  # clean_factors(B)  TODO.
  vsp_factors = list(Z = Z, Y = Y, B=B, recentered = recenter, rescaled = F, Rz = Rz, Ry = Ry, muZ = muZ, muY= muY, pcs = pcs)
  class(vsp_factors) = "vsp"
  return(vsp_factors)
}

vsp_recenter = function(means, V, scree, Rz) return(means %*% V %*% diag(1/scree) %*% Rz)

# clean_factors = function(B){  # TODO
#   #permute rows to make B diagonal-ish
#   apply(B, 1, function(x) return(sign(mean((x-mean(x))^3))))
#   permC = apply(abs(B),1,which.max)
#   permR = apply(abs(B),2,which.max)
#   # if each row has a unique column,
#   if(length(unique(permC)) == length(permC)){
#     B = B[, permC]
#   }
#   B =
#
#
# }

eiv = function(U, fast=TRUE){
  # compute varimax rotation of PCA output U.

  if(fast) {
    lev = sqrt(rowSums(U^2))
    rotHatU = varimax(U[lev>quantile(lev,.9),], normalize=F)$rot
  }
  if(!fast) rotHatU = varimax(U, normalize=F)$rot
  # vv = U %>% apply(1, function(x) return(x*sqrt(sqrt(sum(x^2))))) %>% varimax(normalize=F)
  # rotHatU = vv$rot

  # switch signs to ensure each column of Zhat has positive third moment.
  signss = sign(colSums((U%*%rotHatU)^3))
  rotHatU%*%diag(signss)
  # list(Z = U%*%rotHatU, R = rotHatU)

}







my_line <- function(x,y,...){
  points(x,y,pch = ".", ...)
  lines(c(0,1000), c(0,0), col = "red", lwd = 3, ...)
  lines(c(0,0),c(0,1000),  col = "red", lwd = 3, ...)
  points(x,y,pch = ".",  ...)
  points(0,0, col = "grey", cex = 2, lwd = 2)

}


paneltxt <- function(x, y, labels, cex, font, ...)
{
  lab = labels
  text(.5, .5, lab,  cex=cex, font=font)
}



plot.vsp = function(v, k=NULL, nsamp = 1000, levSamp = TRUE, plotScree = FALSE, whichPlot = "z", transformSqrt= F){
  # kmax = length(v$scree)
  # if(is.null(k))k = kmax
  # if(k > kmax) k = kmax
  # if(plotScree){
  #   plot(2:kmax, v$scree[2:kmax])
  #   abline(9999, k)
  # }
  #
  if(tolower(whichPlot) == "z") X = v$Z
  if(tolower(whichPlot) == "y") X = v$Y
  if(tolower(whichPlot) == "u") X = v$pcs$U
  if(tolower(whichPlot) == "v") X = v$pcs$V

  if(is.null(k)) k = ncol(X)
  n = nrow(X)
  if(nsamp > n){
    samp = 1:n
  }else{
    if(levSamp){
      lev = X[,1:k]^2 %>% rowSums %>% sqrt
      samp = sample(n, nsamp, F, lev)
    }else{
      samp = sample(n, nsamp, F)
    }
  }

  if(plotScree) tmp= scan()

  # colnames(X)[1:k]  = paste0("Lepto:" , round(colMeans(X[,1:k]^4)/colMeans(X[,1:k]^2)^2))
  colnames(X)[1:k]  = round(colMeans((X[,1:k] - colMeans(X[,1:k]))^4)/colMeans((X[,1:k] - colMeans(X[,1:k]))^2)^2)
  # X = as.data.frame(X)
  if(transformSqrt) (X[samp,1:k]) %>% apply(2, function(x) return(sqrt(abs(x))*sign(x))) %>% pairs(lower.panel = my_line, upper.panel = my_line)
  if(!transformSqrt) pairs(X[samp,1:k], labels = colnames(X)[1:k], lower.panel = my_line, upper.panel = my_line, text.panel = paneltxt)


}