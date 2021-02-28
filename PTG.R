## ---------------------------------------------------------------------------------------------------------------------------------
## Author: Minoru Kusaba
## Affiliation: School of Multidisciplinary Sciences Department of Statistical Science SOKENDAI
## Contact: kusaba@ism.ac.jp
## File name: PTG.R
## Task: define functions for PTG algorithm
## Last update: 2019/12/16
## For more information on this algorithm, please see our paper : https://www.nature.com/articles/s41598-021-81850-z (Published: 26 February 2021)
## ---------------------------------------------------------------------------------------------------------------------------------

## require PTG_length_parameter.cpp file and dir_PTG (directory containing "PTG.R" and "PTG_length_parameter.cpp" files)

## add library
library("mvtnorm") #multinormal distribution
library("lpSolve") #to solve bipartite graph matching
library("proxy") #to solve bipartite graph matching
library("tmvtnorm") #truncated multinormal distribution
library("MASS") #to use ginv
library("pheatmap") #for heatmap drawing 
library("scatterplot3d") #to plot 3d points
library("rgl") #interactive 3d plot
library("plot3D") #to plot 3d points

## load Rcpp file (Set directory containing "PTG.R" and "PTG_length_parameter.cpp" files on dir_PTG)
Rcpp::sourceCpp(file = paste(dir_PTG,"PTG_length_parameter.cpp",sep = ""), showOutput = TRUE, rebuild = FALSE) 

## functions for creating grid-points on 1~3D space

## line grid
gtm_pts = function (M) {
  N = M - 1
  matrix(((-N/2):(N/2))/(N/2), ncol = 1)
}

## rectangle grid
gtm_rctg = function (xDim, yDim) {
  if ((xDim < 2) || (yDim < 2) || (yDim != floor(yDim)) || 
      (xDim != floor(xDim))) {
    stop("Invalid grid dimensions")
  }
  r1 = 0:(xDim - 1)
  r2 = (yDim - 1):0
  fx = function(x, y) return(x)
  fy = function(x, y) return(y)
  X = outer(r1, r2, fx)
  Y = outer(r1, r2, fy)
  grid = cbind(matrix(X, ncol = 1), matrix(Y, ncol = 1))
  maxVal = max(grid)
  grid = grid * (2/maxVal)
  maxXY = apply(grid, 2, max)
  grid[, 1] = grid[, 1] - maxXY[1]/2
  grid[, 2] = grid[, 2] - maxXY[2]/2
  grid
}

## cube grid (when scale_cor = T,distances of each grid-points are scaled to be equal)
gtm_cube = function (xDim, yDim, zDim, scale_cor = F) {
  if ((xDim < 2) || (yDim < 2) || (zDim < 2) || (zDim != floor(zDim)) || (yDim != floor(yDim)) || 
      (xDim != floor(xDim))) {
    stop("Invalid grid dimensions")
  }
  Z = gtm_pts(zDim)
  XY = gtm_rctg(xDim, yDim)
  XYZ = matrix(0, ncol = 3, nrow = xDim*yDim*zDim)
  for(i in 1:zDim){
    layer_i = cbind(XY, Z[i])
    XYZ[(1+xDim*yDim*(i-1)):(xDim*yDim*i),] = layer_i
  }
  if(scale_cor){
    dis_lyr = 2/(max(xDim,yDim)-1)
    dis_z = 2/(zDim-1)
    XYZ[,3] = XYZ[,3]*(dis_lyr/dis_z)
  }
  XYZ
}

## cylinder grid (cDim define number of grid on one layer:zDim define number of the layers)
gtm_cld = function (cDim, zDim, scale_cor = F) {
  if ((cDim < 2) || (zDim < 2) || (zDim != floor(zDim)) || (cDim != floor(cDim))) {
    stop("Invalid grid dimensions")
  }
  Z = gtm_pts(zDim)
  angles = seq(from = 0, to = 2*pi-(2*pi/cDim), by = 2*pi/cDim)
  XY = cbind(cos(angles), sin(angles))
  XYZ = matrix(0, ncol = 3, nrow = cDim*zDim)
  for(i in 1:zDim){
    layer_i = cbind(XY, Z[i])
    XYZ[(1+cDim*(i-1)):(cDim*i),] = layer_i
  }
  if(scale_cor){
    dis_cir = 2*pi/cDim
    dis_z = 2/(zDim-1)
    XYZ[,3] = XYZ[,3]*(dis_cir/dis_z)
  }
  XYZ
}

## cone grid (when scale_cor = T,distances of each grid-points are scaled to be equal:zDim define number of the layers)
gtm_cone = function(zDim,scale_cor = T){
  lens = seq(from = 1/zDim,to = 1,by = 1/zDim)
  cnum = 4*(seq(1:zDim))
  z = seq(from = 1-2/zDim,to = -1,by = -2/zDim)
  grid = matrix(c(0,0,1),ncol = 3,nrow = 1)
  for(i in 1:length(cnum)){
    num = cnum[i]
    XYZ = cbind(lens[i]*cos((seq(from = 0,to = 2*pi - 2*pi/num,by = 2*pi/num))),lens[i]*sin(seq(from = 0,to = 2*pi - 2*pi/num,by = 2*pi/num)),rep(z[i],num))
    grid = rbind(grid,XYZ)
  }
  if(scale_cor){
    scale = sqrt((pi/2)^2-1)/2
    grid[,3] = grid[,3]*scale
  }
  grid
}

## general function for creating various shape of grid points
crt_grid = function(grid_str, xDim = NA, yDim = NA, zDim = NA, scale_cor = F, cDim = NA){
  if((grid_str != "line") && (grid_str != "rctg") && (grid_str != "cube") && (grid_str != "cld") && (grid_str != "cone")){
    stop("Structure name is Invalid : valid names are line,rctg,cube,cld,and,cone")
  }
  if(grid_str == "line"){
    grid = gtm_pts(xDim)
  }
  if(grid_str == "rctg"){
    grid = gtm_rctg(xDim,yDim)
  }
  if(grid_str == "cube"){
    grid = gtm_cube(xDim,yDim,zDim,scale_cor)
  }
  if(grid_str == "cld"){
    grid = gtm_cld(cDim,zDim,scale_cor)
  }
  if(grid_str == "cone"){
    grid = gtm_cone(zDim = zDim,scale_cor = scale_cor)
  }
  return(list(grid = grid, grid_str = grid_str, xDim = xDim, yDim = yDim, zDim = zDim, scale_cor = scale_cor, cDim = cDim))
}

## general function for interpolating various shape of grid points (n_sep is number of separation between grid-points:round_num is used to avoid numerical error)
itp_grid = function(grid_sum, n_sep = 1,round_num = 12){
  # obtain values
  xDim = grid_sum$xDim
  yDim = grid_sum$yDim
  zDim = grid_sum$zDim
  cDim = grid_sum$cDim
  scale_cor = grid_sum$scale_cor
  grid_str = grid_sum$grid_str
  # create igrid
  if(grid_str == "line"){
    xDim = (n_sep + 1)*xDim - n_sep
    igrid = gtm_pts(xDim)
  }
  if(grid_str == "rctg"){
    xDim = (n_sep + 1)*xDim - n_sep
    yDim = (n_sep + 1)*yDim - n_sep
    igrid = gtm_rctg(xDim = xDim,yDim = yDim)
  }
  if(grid_str == "cube"){
    xDim = (n_sep + 1)*xDim - n_sep
    yDim = (n_sep + 1)*yDim - n_sep
    zDim = (n_sep + 1)*zDim - n_sep
    igrid = gtm_cube(xDim = xDim,yDim = yDim,zDim = zDim,scale_cor = scale_cor)
  }
  if(grid_str == "cld"){
    cDim = (n_sep + 1)*cDim
    zDim = (n_sep + 1)*zDim - n_sep
    igrid = gtm_cld(cDim = cDim,zDim = zDim,scale_cor = scale_cor)
  } 
  if(grid_str == "cone"){
    zDim = (n_sep+1)*zDim
    igrid = gtm_cone(zDim = zDim,scale_cor = scale_cor)
  }
  # get original positions of the grid in igrid
  grid = grid_sum$grid
  L = ncol(grid)
  K = nrow(grid)
  obs = numeric(K)
  if(L == 1){
    for(i in 1:K){
      obs[i] = which((round(abs(igrid[,1] - grid[i,1]),round_num) == 0) == 1)
    }
  }
  if(L == 2){
    for(i in 1:K){
      obs[i] = which(((round(abs(igrid[,1] - grid[i,1]),round_num) == 0) * (round(abs(igrid[,2] - grid[i,2]),round_num) == 0)) == 1)
    }
  }
  if(L == 3){
    for(i in 1:K){
      obs[i] = which(((round(abs(igrid[,1] - grid[i,1]),round_num) == 0) * (round(abs(igrid[,2] - grid[i,2]),round_num) == 0) * (round(abs(igrid[,3] - grid[i,3]),round_num) == 0)) == 1)
    }
  }
  return(list(grid = grid,igrid = igrid,obs = obs, grid_str = grid_str, xDim = xDim, yDim = yDim, zDim = zDim, scale_cor = scale_cor, cDim = cDim))
}

## define kernel function 
gpcov = function(d,v,l){
  v*exp(-d/(2*l))
}

## initialize parameters from prior distribution (epsilon is used to make covariance matrix Cst positive-semidefinite:Ch_epsilon is Ch version)
initialize_parameters_random = function(grid_sum,data,shape,rate,epsilon,v_g,l_g,v_r,l_r,Ch_epsilon){
  grid = grid_sum$grid
  # set hyper parameters
  D = ncol(data) #dimension of data
  N =  nrow(data) #number of data
  L = ncol(grid) #dimension of nodes
  K = nrow(grid) #number of nodes
  
  d0 = shape #shape parameter of gamma distribution
  s0 = rate #rate parameter of gamma distribution
  d = d0 + N*D/2
  
  # derive Cst covariance matrix
  sqdist = dist(x = grid,y = grid)^2 #squared distance matrix for nodes
  class(sqdist) = "matrix"
  pCst_g = gpcov(sqdist,v_g,l_g) #fixed
  Cst_g = (pCst_g + v_g * epsilon * diag(K)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_g = solve(Cst_g)
  inCst_g = (inCst_g + t(inCst_g))/2
  
  pCst_r = gpcov(sqdist,v_r,l_r) #fixed
  Cst_r = (pCst_r + v_r * epsilon * diag(K)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_r = solve(Cst_r)
  inCst_r = (inCst_r + t(inCst_r))/2
  
  # set initial parameters
  beta0 = rgamma(n = 1,shape = d0,rate = s0) #initial variance
  g0 = rtmvnorm(n = 1, mean = rep(0,K),sigma = Cst_g,lower=rep(0, length = K),upper=rep( Inf, length = K)) #initial variance parameters
  r0 = rmvnorm(n = 1,mean = rep(0,K),sigma = Cst_r) #initial length scale parameters
  # generate H0
  Ch = get_Ch(r0,sqdist,epsilon = Ch_epsilon,L = L)
  H0 = rmvnorm(D,mean = rep(0,K),sigma = Ch) #initial H
  return (list(beta0 = beta0,g0 = g0,H0 = H0,r0 = r0,Cst_r = Cst_r,inCst_r = inCst_r,sqdist = sqdist,d = d,Cst_g = Cst_g,inCst_g = inCst_g))
}

## main body of algorithm
## function for GTM_LDLV
GTM_LDLV = function(grid_sum,data,shape,rate,epsilon,Iter,burn.in.samples,maxit_BFGS,maxit_SANN,par_scale,v_r,l_r,v_g,l_g,Ch_epsilon){
  # get initial parameters  
  init = initialize_parameters_random(grid_sum = grid_sum,data = data,shape = shape,rate = rate,epsilon = epsilon,v_r = v_r,l_r = l_r,v_g = v_g,l_g = l_g,Ch_epsilon=Ch_epsilon)
  # get hyper parameters
  grid = grid_sum$grid
  D = ncol(data) #dimension of data
  N =  nrow(data) #number of data
  L = ncol(grid) #dimension of nodes
  K = nrow(grid) #number of nodes
  d = init$d #shape parameter of gamma distribution
  s0 = rate #rate parameter of gamma distribution
  # set initial parameters
  beta = init$beta0
  g = init$g0
  H = init$H0
  r = init$r0
  Cst_r = init$Cst_r
  inCst_r = init$inCst_r
  Cst_g = init$Cst_g
  inCst_g = init$inCst_g
  sqdist = init$sqdist
  # to memorize parameters thought algorithm
  C_memory = array(0,dim = c(K,N,Iter)) # Connection matrix between data-points and projected-grid
  g_memory = matrix(0,ncol = K,nrow = Iter)
  r_memory = matrix(0,ncol = K,nrow = Iter)
  Y_memory = array(0,dim = c(D,K,Iter))
  H_memory = array(0,dim = c(D,K,Iter))
  beta_memory = numeric(Iter)
  Hesse_memory = numeric(Iter)
  # set frequentry used values
  rep_D = rep(0,D)
  rep_K = rep(0,K)
  diag_D = diag(D)
  lower=rep(0, length = K)
  upper=rep( Inf, length = K)
  parscale = rep(par_scale,K)
  
  for(iter in 1:Iter){
    # sampling of C
    g = as.vector(g)
    Y = H %*% diag(g)
    Y_memory[,,iter] = Y
    tY = t(Y)
    
    # sampling of beta
    dis = dist(x = tY, y = data)^2 #squared distance matrix between data and prototypes
    edis = exp(-(beta/2)*dis)
    sums = matrix(colSums(edis),ncol = N,nrow = K,byrow = T)
    probs = edis/sums
    rows_K = 1:K
    vec_C = apply(X = probs,MARGIN = 2,FUN = function(x){sample(rows_K,size = 1,prob = x)})
    C = matrix(0,ncol = N,nrow = K)
    C[cbind(vec_C,1:N)] = 1
    # sum of square error
    sse = sum(C * dis)
    C_memory[,,iter] = C
    s = s0 + sse/2
    beta = rgamma(n = 1,shape = d,rate = s) #sampling of beta
    beta_memory[iter] = beta
    
    # sampling of g
    vecG = rowSums(C)
    G = diag(vecG)
    diagh = diag(colSums(H * H))
    sigmag = ginv((beta * G %*% diagh) + inCst_g)
    sigmag = (sigmag + t(sigmag))/2
    vecCtXH = diag(C %*% data %*% H)
    myug = as.vector(beta * sigmag %*% vecCtXH)
    g = rtmvnorm(n = 1, mean = myug,sigma = sigmag,lower=lower,upper=upper,algorithm = "gibbs",burn.in.samples = burn.in.samples,start.value = g)
    g = as.vector(g)
    if(sum(is.na(g)) == 0){
      g = g
    }else{
      g = abs(myug)
    }
    g_memory[iter,] = g
    diagg = diag(g)
    
    # sampling of H
    Ch = get_Ch(r,sqdist,epsilon = Ch_epsilon,L = L)
    inCh = solve(Ch)
    sigmah = ginv((beta * G %*% (diagg^2))+inCh)
    sigmah = (sigmah + t(sigmah))/2
    myuhcol = beta * sigmah %*% diagg %*% C %*% data
    myuh = t(myuhcol)
    Hpre = rmvnorm(D,mean = rep_K,sigma = sigmah)
    H = Hpre + myuh
    H_memory[,,iter] = H
    
    # sampling of r
    r = as.vector(r)
    
    loglike_r = function(r){
      give_logDens(r = r,H = H,inCst = inCst_r,d = sqdist,epsilon = Ch_epsilon,L = L)
    }
    
    gr_loglike = function(r){
      give_dr(r = r,d = sqdist,H = H,inCst = inCst_r,epsilon = Ch_epsilon,L = L)
    }
    
    if(sum(is.na(gr_loglike(r))) != 0){
      r = rmvnorm(1,rep_K,Cst_r)
    }
    
    value = try(optim(par = r,fn = loglike_r,gr = gr_loglike,method = "BFGS",control = list(fnscale = -1,maxit = maxit_BFGS,parscale = parscale))$par)
    if(class(value) == "try-error"){
      r = optim(par = r,fn = loglike_r,method = "SANN",control = list(fnscale = -1,maxit = maxit_SANN))$par
    }else{
      r = value
    }
    # r is optimized
    # hessian part
    hesse = give_ddr(r = r,d = sqdist,H = H,inCst = inCst_r,epsilon = Ch_epsilon,L = L)
    s_r = give_dr(r = r,d = sqdist,H = H,inCst = inCst_r,epsilon = Ch_epsilon,L = L)
    V_r = solve(-hesse)
    m_r = r + V_r %*% s_r
    # sampling new r
    r_new = rmvnorm(1,m_r,V_r)
    unif = runif(1,min = 0,max= 1)
    Bool = (exp(loglike_r(r_new) - loglike_r(r) + dmvnorm(r,mean = m_r,sigma = V_r,log = T) - dmvnorm(r_new,mean = m_r,sigma = V_r,log = T))>unif)
    Hesse_memory[iter] = Bool
    r_n = Bool * (r_new - r) + r
    if(sum(is.na(r_n)) != 0){
      r = r
    }else{
      r = r_n
    }
    # get randomized r
    r = as.vector(r)
    r_memory[iter,] = r
  }
  return(list(accepted_r = Hesse_memory,beta = beta_memory,g = g_memory,H = H_memory,r = r_memory,Y = Y_memory,C = C_memory,Cst_r = Cst_r,Cst_g = Cst_g))
}

## get ensemble of GTM_LDLV
get_sample_means = function(result,Cut_point,data){
  Iter = nrow(result$g)
  K = ncol(result$r)
  D = nrow(result$Y)
  
  Y_sum = matrix(0,D,K)
  r_sum = numeric(K)
  g_sum = numeric(K)
  
  for(i in Cut_point:Iter){
    Y_sum = result$Y[,,i] + Y_sum
    r_sum = result$r[i,] + r_sum
    g_sum = result$g[i,] + g_sum
  }
  
  Y_sam = Y_sum/(Iter-Cut_point+1)
  r_sam = r_sum/(Iter-Cut_point+1)
  g_sam = g_sum/(Iter-Cut_point+1)
  beta_sam = mean(result$beta[Cut_point:Iter])
  
  return(list(r = r_sam,g = g_sam,Y = Y_sam,beta = beta_sam))
}

## GTM_LDLV_ens (combine GTM_LDLV and ensamble function:This is the first step of PTG:Constrain the estimation by lambda)
GTM_LDLV_ens = function(grid_sum,data,lambda = 1,shape = 1,rate = 1,epsilon = 0.01,Iter = 10000,burn.in.samples = 10,maxit_BFGS = 10,maxit_SANN = 100,par_scale = 0.001,Cut_point = 5000,Ch_epsilon = 0.0001){
  result = GTM_LDLV(grid_sum = grid_sum,data = data,shape = shape,rate = rate,epsilon = epsilon,Iter = Iter,burn.in.samples = burn.in.samples,maxit_BFGS = maxit_BFGS,maxit_SANN = maxit_SANN,par_scale = par_scale,v_r = 1/lambda,l_r = lambda,v_g = 1/lambda,l_g = lambda,Ch_epsilon = Ch_epsilon)
  result_ens = get_sample_means(result = result,Cut_point = Cut_point,data=data)
  
  tY = t(result_ens$Y)
  dis = dist(x = tY, y = data)^2 #squared distance matrix between data and prototypes
  grid = grid_sum$grid
  beta = result_ens$beta 
  edis = exp(-(beta/2)*dis)
  sums = matrix(colSums(edis),ncol = nrow(data),nrow = nrow(grid),byrow = T)
  probs = edis/sums
  x = t(probs) %*% grid
  
  return(list(result = result,result_ens = result_ens,map = x,grid_sum = grid_sum))
}

## interpolate the parameters learned in GTM-LDLV by gaussian process regression (second step of PTG) 
interpolation = function(result_ens,grid_sum,v_g,l_g,v_r,l_r,epsilon,n_sep,epsilon_Cf = 0.0001,round_num = 12){
  # get grid and igrid
  igrid_sum = itp_grid(grid_sum = grid_sum,n_sep = n_sep,round_num = round_num)
  grid = igrid_sum$grid
  igrid = igrid_sum$igrid
  obs = igrid_sum$obs
  # get values
  FULL = nrow(igrid) #full number
  N = nrow(grid) #number of grid-points(obs)
  K = FULL - N #number of positions to be interpolated
  L = ncol(grid) #dimension of grid
  # get matricies(r and g)
  sqdist = dist(x = igrid,y = igrid)^2 #squared distance matrix for nodes
  class(sqdist) = "matrix"
  pCst_g = gpcov(sqdist,v_g,l_g) #fixed
  Cst_g = (pCst_g + v_g * epsilon * diag(FULL)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_g = solve(Cst_g[obs,obs])
  inCst_g = (inCst_g + t(inCst_g))/2
  
  pCst_r = gpcov(sqdist,v_r,l_r) #fixed
  Cst_r = (pCst_r + v_r * epsilon * diag(FULL)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_r = solve(Cst_r[obs,obs])
  inCst_r = (inCst_r + t(inCst_r))/2
  
  # get interpolated r and g
  new_r = numeric(FULL)
  new_r[obs] = result_ens$r
  new_r[-obs] = Cst_r[-obs,obs] %*% inCst_r %*% result_ens$r
  
  new_g = numeric(FULL)
  new_g[obs] = result_ens$g
  new_g[-obs] = Cst_g[-obs,obs] %*% inCst_g %*% result_ens$g
  
  # get interpolated Cf
  Ch = get_Ch(new_r,sqdist,epsilon = 0,L = L)
  G = new_g %*% t(new_g)
  Cf = G * Ch
  
  #get interpolated Y
  Y_obs = t(result_ens$Y)
  inv_iCf = solve(Cf[obs,obs] + diag(length(obs)) * epsilon_Cf)
  inv_iCf = (inv_iCf + t(inv_iCf))/2
  Y_itp = t(Cf[-obs,obs] %*% inv_iCf %*% Y_obs)
  iY = matrix(0,ncol = ncol(Cf),nrow = ncol(Y_obs))
  iY[,obs] = t(Y_obs)
  iY[,-obs] = Y_itp
  
  return(list(igrid_sum = igrid_sum,ir = new_r,ig = new_g,iCf = Cf,iY = iY,Cst_g = Cst_g,Cst_r = Cst_r,beta = result_ens$beta))
}

## tuning the interpolated parameters with one-to-one constraints on latent variables (third step of PTG)
itp_tuning = function(itp,data,shape = 1,rate = 1,Iter = 30,maxit_BFGS = 10,maxit_SANN = 100,par_scale = 0.001,Ch_epsilon = 0.0001){
  # get initial parameters  
  Y = itp$iY
  # get hyper parameters
  grid = itp$igrid_sum$igrid
  D = ncol(data) #dimension of data
  N =  nrow(data) #number of data
  L = ncol(grid) #dimension of grid
  K = nrow(grid) #number of grid
  d0 = shape #shape parameter of gamma distribution
  d = d0 + N*D/2
  s0 = rate #rate parameter of gamma distribution
  beta = itp$beta
  g = itp$ig
  r = itp$ir
  H = Y %*% solve(diag(as.vector(g)))
  Cst_r = itp$Cst_r
  inCst_r = solve(Cst_r)
  inCst_r = (inCst_r + t(inCst_r))/2
  Cst_g = itp$Cst_g
  inCst_g = solve(Cst_g)
  inCst_g = (inCst_g + t(inCst_g))/2
  # get distance matrix
  sqdist = dist(x = grid,y = grid)^2 #squared distance matrix for nodes
  class(sqdist) = "matrix"
  # set objects 
  C_memory = array(0,dim = c(K,N,Iter))
  g_memory = matrix(0,ncol = K,nrow = Iter)
  r_memory = matrix(0,ncol = K,nrow = Iter)
  Y_memory = array(0,dim = c(D,K,Iter))
  H_memory = array(0,dim = c(D,K,Iter))
  beta_memory = numeric(Iter)
  # set repeatedly used values
  rep_D = rep(0,D)
  rep_K = rep(0,K)
  diag_D = diag(D)
  lower=rep(0, length = K)
  upper=rep( Inf, length = K)
  parscale = rep(par_scale,K)
  # start iteration
  for(iter in 1:Iter){
    # with one-to-one constraints
    # calculation of C
    g = as.vector(g)
    Y = H %*% diag(g)
    Y_memory[,,iter] = Y
    tY = t(Y)
    # bipartite graph matching
    dis = dist(x = tY, y = data)
    lenrow = K
    lencol = N
    row.signs = rep ("<=", lenrow)
    row.rhs = rep(1,lenrow)
    col.signs = rep ("=", lencol)
    col.rhs = rep(1,lencol)
    C = lp.transport (dis, "min", row.signs, row.rhs, col.signs, col.rhs)$solution #calculation of C
    # calculation of beta
    sse = sum(C * (dis)^2)
    C_memory[,,iter] = C
    s = s0 + sse/2
    beta = (d-1)/s #calculation of beta
    beta_memory[iter] = beta
    # calculation of g
    vecG = rowSums(C)
    G = diag(vecG)
    diagh = diag(colSums(H * H))
    sigmag = ginv((beta * G %*% diagh) + inCst_g)
    sigmag = (sigmag + t(sigmag))/2
    vecCtXH = diag(C %*% data %*% H)
    myug = as.vector(beta * sigmag %*% vecCtXH)
    g = abs(myug) #calculation of g
    g_memory[iter,] = g
    diagg = diag(g)
    # calculation of H
    Ch = get_Ch(r,sqdist,epsilon = Ch_epsilon,L = L)
    inCh = solve(Ch)
    sigmah = ginv((beta * G %*% (diagg^2)) + inCh)
    sigmah = (sigmah + t(sigmah))/2
    myuhcol = beta * sigmah %*% diagg %*% C %*% data
    myuh = t(myuhcol)
    H = myuh
    H_memory[,,iter] = H
    # sampling of r
    r = as.vector(r)
    loglike_r = function(r){
      give_logDens(r = r,H = H,inCst = inCst_r,d = sqdist,epsilon = Ch_epsilon,L = L)
    }
    gr_loglike = function(r){
      give_dr(r = r,d = sqdist,H = H,inCst = inCst_r,epsilon = Ch_epsilon,L = L)
    }
    if(sum(is.na(gr_loglike(r))) != 0){
      r = rmvnorm(1,rep_K,Cst_r)
    }
    value = try(optim(par = r,fn = loglike_r,gr = gr_loglike,method = "BFGS",control = list(fnscale = -1,maxit = maxit_BFGS,parscale = parscale))$par)
    if(class(value) == "try-error"){
      r = optim(par = r,fn = loglike_r,method = "SANN",control = list(fnscale = -1,maxit = maxit_SANN))$par
    }else{
      r = value
    }
    r = as.vector(r)
    r_memory[iter,] = r
  }
  return(list(beta = beta_memory,g = g_memory,H = H_memory,r = r_memory,Y = Y_memory,C = C_memory,Cst_r = Cst_r,Cst_g = Cst_g))
}

## calculate log-likelihood through algorithm
get_likelihood = function(result,data){
  C = result$C
  Y = result$Y
  beta = result$beta
  Iter = length(beta)
  K = nrow(C)
  D = nrow(Y)
  N = ncol(C)
  likelihood = numeric(Iter)
  for(iter in 1:Iter){
    s =  t(Y[,,iter] %*% C[,,iter])
    likelihood[iter] = sum(dmvnorm(x = data - s,mean = rep(0,D),sigma = diag(D)/beta[iter],log = T)) + log(K) * (-N)
  }
  likelihood
}

## function for PTG (summarize all 3step algorithms:combine the three algorithms) 
PTG = function(grid_sum,data,n_sep = 1,lambda = 1,shape = 1,rate = 1,epsilon = 0.01,Iter = 10000,burn.in.samples = 10,maxit_BFGS = 10,maxit_SANN = 100,par_scale = 0.001,Cut_point = 5000,epsilon_Cf = 0.0001,shape_tun = 1,rate_tun = 1,Iter_tun = 10,maxit_BFGS_tun = 100,maxit_SANN_tun = 1000,par_scale_tun = 0.001,Ch_epsilon = 0.0001){
  ldlv = GTM_LDLV_ens(grid_sum = grid_sum,data = data,lambda = lambda,shape = shape,rate = rate,epsilon = epsilon,Iter = Iter,burn.in.samples = burn.in.samples,maxit_BFGS = maxit_BFGS,maxit_SANN = maxit_SANN,par_scale = par_scale,Cut_point = Cut_point,Ch_epsilon = Ch_epsilon)
  itp = interpolation(result_ens =  ldlv$result_ens,grid_sum = grid_sum,v_r = 1/lambda,l_r = lambda,v_g = 1/lambda,l_g = lambda,epsilon = epsilon,n_sep = n_sep,epsilon_Cf = epsilon_Cf)
  tun = itp_tuning(itp = itp,data = data,shape = shape_tun,rate = rate_tun,Iter = Iter_tun,maxit_BFGS = maxit_BFGS_tun,maxit_SANN = maxit_SANN_tun,par_scale = par_scale_tun,Ch_epsilon = Ch_epsilon)
  likelihood_tun = get_likelihood(result = tun,data = data)
  best_num = which(max(likelihood_tun) == likelihood_tun)
  return(list(ldlv = ldlv,itp = itp,tun = tun,best_num = best_num))
}

## get and summurize learned results of PTG (map, responsibility and others:map is result of dimension-reduction:proj is projection of map to data-space:resp is responsibility of each data-points for each grid-points:Y is projection of igrid to data-space)
get_resultof_PTG = function(result,data){
  num = result$best_num
  tun = result$tun
  Y = tun$Y[,,num]
  igrid = result$itp$igrid_sum$igrid
  C = tun$C[,,num]
  map = t(C) %*% igrid
  proj = t(C) %*% t(Y)
  
  tY = t(Y)
  dis = dist(x = tY, y = data)^2 #squared distance matrix between data and prototypes
  beta = tun$beta[num] 
  edis = exp(-(beta/2)*dis)
  N = nrow(data)
  K = ncol(Y)
  sums = matrix(colSums(edis),ncol = N,nrow = K,byrow = T)
  probs = edis/sums
  resp = t(probs)
  
  rownames(map) = rownames(data)
  colnames(proj) = colnames(data)
  rownames(proj) = rownames(data)
  tC = t(C)
  rownames(tC) = rownames(data)
  
  #add ldlv
  ldlv_map = result$ldlv$map
  
  return(list(map = map,C = tC,proj = proj,resp = resp,Y = tY,ldlv_map = ldlv_map))
}

## for creating PTG property landscapes
## calculate prediction surface of learned parameters by gaussian regression for interpretation of grid (lambda should be same to PTGs-lambda,n_sep is used for making grid for kriging)
kriging_calculation = function(result,n_sep = 3,lambda = 1,epsilon = 0.01,epsilon_Cf = 0.0001,round_num = 12){
  # get interpolated grid
  igrid_sum = result$itp$igrid_sum
  # object processing
  igrid_sum$grid = igrid_sum$igrid
  # increase grid for making prediction surface
  igrid_sum = itp_grid(grid_sum = igrid_sum,n_sep = n_sep,round_num = 12)
  # set values
  v_r = 1/lambda
  l_r = lambda
  v_g = 1/lambda
  l_g = lambda
  # get kriging grid
  grid = igrid_sum$grid
  igrid = igrid_sum$igrid
  obs = igrid_sum$obs
  # get values
  FULL = nrow(igrid) #full number
  N = nrow(grid) #number of grid-points(obs)
  K = FULL - N #number of positions to be interpolated
  L = ncol(grid) #dimension of grid
  # get matricies(r and g)
  sqdist = dist(x = igrid,y = igrid)^2 #squared distance matrix for nodes
  class(sqdist) = "matrix"
  pCst_g = gpcov(sqdist,v_g,l_g) #fixed
  Cst_g = (pCst_g + v_g * epsilon * diag(FULL)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_g = solve(Cst_g[obs,obs])
  inCst_g = (inCst_g + t(inCst_g))/2
  pCst_r = gpcov(sqdist,v_r,l_r) #fixed
  Cst_r = (pCst_r + v_r * epsilon * diag(FULL)) #epsilon is needed to guarantee nonsingularity of Cst matrix
  inCst_r = solve(Cst_r[obs,obs])
  inCst_r = (inCst_r + t(inCst_r))/2
  # get learned parameters by PTG
  num = result$best_num
  tun = result$tun
  g = tun$g[num,]
  r = tun$r[num,]
  Y = tun$Y[,,num]
  # get interpolated r and g
  new_r = numeric(FULL)
  new_r[obs] = r
  new_r[-obs] = Cst_r[-obs,obs] %*% inCst_r %*% r
  new_g = numeric(FULL)
  new_g[obs] = g
  new_g[-obs] = Cst_g[-obs,obs] %*% inCst_g %*% g
  # get interpolated Cf
  Ch = get_Ch(new_r,sqdist,epsilon = 0,L = L)
  G = new_g %*% t(new_g)
  Cf = G * Ch
  # get interpolated Y
  Y_obs = t(Y)
  inv_iCf = solve(Cf[obs,obs] + diag(length(obs)) * epsilon_Cf)
  inv_iCf = (inv_iCf + t(inv_iCf))/2
  Y_itp = t(Cf[-obs,obs] %*% inv_iCf %*% Y_obs)
  iY = matrix(0,ncol = ncol(Cf),nrow = ncol(Y_obs))
  iY[,obs] = t(Y_obs)
  iY[,-obs] = Y_itp
  return(list(Y = iY,g = new_g,r = new_r,krg_grid = igrid))
}

## show result of kriging (based on scatterplot3D)
show_kriging = function(krg,data,col_num,which_prm = "Y",pch = 15,cex = 2,phi = 0,theta = 30){
  grid = krg$krg_grid
  L = ncol(grid)
  K = nrow(grid)
  if(L == 1){
    if(which_prm == "Y"){
      Y = t(krg$Y)
      colnames(Y) = colnames(data) 
      colvar = Y[,col_num]
      scatter2D(x = grid[,1], y = rep(0,K), bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = colnames(Y)[col_num],xlim = c(-1,1), ylim = c(-1,1))
    }
    if(which_prm == "g"){
      colvar = krg$g
      scatter2D(x = grid[,1], y = rep(0,K), bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = "g(u)",xlim = c(-1,1), ylim = c(-1,1))
    }
    if(which_prm == "r"){
      colvar = exp(krg$r)
      scatter2D(x = grid[,1], y = rep(0,K), bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = "l(u)",xlim = c(-1,1), ylim = c(-1,1))
    }
  }
  if(L == 2){
    if(which_prm == "Y"){
      Y = t(krg$Y)
      colnames(Y) = colnames(data) 
      colvar = Y[,col_num]
      scatter2D(x = grid[,1], y = grid[,2], bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = colnames(Y)[col_num],xlim = c(-1,1), ylim = c(-1,1))
    }
    if(which_prm == "g"){
      colvar = krg$g
      scatter2D(x = grid[,1], y = grid[,2], bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = "g(u)",xlim = c(-1,1), ylim = c(-1,1))
    }
    if(which_prm == "r"){
      colvar = exp(krg$r)
      scatter2D(x = grid[,1], y = grid[,2], bty = "g",colvar = colvar
                , pch = pch, cex = cex, main = "l(u)",xlim = c(-1,1), ylim = c(-1,1))
    }
  }
  if(L == 3){
    if(which_prm == "Y"){
      Y = t(krg$Y)
      colnames(Y) = colnames(data) 
      colvar = Y[,col_num]
      scatter3D(x = grid[,1], y = grid[,2], z = grid[,3], phi = phi,theta = theta, bty = "g", colvar = colvar
                , pch = pch, cex = cex, ticktype = "detailed", main = colnames(Y)[col_num])
    }
    if(which_prm == "g"){
      colvar = krg$g
      scatter3D(x = grid[,1], y = grid[,2], z = grid[,3], phi = phi,theta = theta, bty = "g", colvar = colvar
                , pch = pch, cex = cex, ticktype = "detailed", main = "g(u)")
    }
    if(which_prm == "r"){
      colvar = exp(krg$r)
      scatter3D(x = grid[,1], y = grid[,2], z = grid[,3], phi = phi,theta = theta, bty = "g", colvar = colvar
                , pch = pch, cex = cex, ticktype = "detailed", main = "l(u)")
    }
  }
}
