###############################
##### Autor: Alexis Ayala #####
##### Fecha: Agosto 2025  #####
###############################

##### Funciones definidas

# ------------------------------------------
# Kernels & Anchos de Banda
# ------------------------------------------

# Kernel Triangular
trianKernel <- function(u) {
  return(ifelse(abs(u) <= 1, 1 - abs(u), 0))
}

# Kernel Gaussiano
gausKernel <- function(u) {
  return (apply(dnorm(u), 1, prod))
}

# Kernel Epanechnikov
epaneKernel <- function(u) {
  return(ifelse(abs(u) <= 1, 0.75*(1-u^2) , 0))
}

# Kernel Uniforme
uKernel <- function(u) {
  return(ifelse(abs(u) <= 1, 0.5, 0))
}

# Ancho de Banda de Silverman
silverman <- function(x) {
  IRQ <- quantile(x, c(0.25, 0.75))
  u1 <- sqrt(var(x))
  u2 <- (IRQ[2] - IRQ[1])/1.34
  return(1.06*min(u1,u2)*length(x)^(-1/5))
}

# Estimación de densidad de kernel
calculaKDE <- function(x, kernel, h, puntos) {
  U <- -sweep(data, 2, x, "-")/ h
  K <- gausKernel(U)
  kde <- 1/(n*h) * sum(K)
  
  return(kde)
}

# --------------------------------------------
# Pesos de Nadaraya y Estimación
# --------------------------------------------

# Weights Continuas (Pesos de Nadaraya)
W <- function(x0, X, h) {
  U <- -sweep(X, 2, x0, "-")/ h
  K <- gausKernel(U)
  W <- K / sum(K)
  return(W)
}

# Weights Categoricas (Pesos de Nadaraya)
WCat <- function(x0, X) {
  
  match_cat <- rowSums(X == matrix(x0, nrow(X), length(x0), byrow = TRUE)) == length(x0)
  K <- as.numeric(match_cat)
  
  W <- K / sum(K)
  return(W)
}

# Weights Extendido Cont + Cat (Pesos de Nadaraya)
WMix <- function(x0_cont, x0_cat, X_cont, X_cat, bw) {
  # Parte continua
  U <- -sweep(X_cont, 2, x0_cont, "-")/ bw
  K_cont <- gausKernel(U)
  
  # Parte categórica (peso 1 si coincide, 0 si no)
  match_cat <- rowSums(X_cat == matrix(x0_cat, nrow(X_cat), length(x0_cat), byrow = TRUE)) == length(x0_cat)
  K_cat <- as.numeric(match_cat)
  
  # Peso total
  K <- K_cont * K_cat
  W <- K / sum(K)
  return(W)
}

# Estimador NW
nadarayaWatson <- function(weights, y){
  est <- weights %*% y 
  return(est)
}

# --------------------------------------------
# Kaplan-Meier considerando IC
# --------------------------------------------

# Estimador de Kaplan-Meier
KM <- function(tiempo, censura) {
  # Creamos un data frame para sortear los datos
  df <- data.frame(tiempo = tiempo, censura = censura)
  df <- df[order(df$tiempo), ]
  
  # Encontramos los tiempos observados (no censurados)
  unicos <- unique(df$tiempo[df$censura == 1])
  
  # Hacemos un ciclo para encontrar la probabilidad en cada evento observado
  # Iniciamos con P(Surv) = 1
  KM <- 1
  prob <- 1
  err <- 1
  N <- length(unicos)
  suma <- 0
  
  if(N == 0) {
    res <- data.frame("tiempo" = c(1), "prob" = c(1))
  }else{
    for (i in 1:N) {
      # Nos paramos en un tiempo
      t <- unicos[i]
      
      # Vemos eventos observados
      d <- sum(df$tiempo == t & df$censura == 1)
      
      # Vemos eventos en riesgo
      n <- sum(df$tiempo >= t)
      
      ## Probabilidad de Supervivencia
      # Calculamos formula de KM
      KMForm <- 1 - d/n
      
      # Obtenemos probabilidad de supervivencia actualizada
      KM <- KM * KMForm
      
      # Guardamos la probabilidad
      prob[i] <- KM
      
      ## Intervalos de Confianza
      # Formula de Greenwood
      CI <- d/(n*(n-d))
      suma <- suma + CI
      
      # Guardamos el error
      err[i] <- suma
    }
    # Checamos el último caso donde se puede llegar a 0
    if(prob[N] == 0){
      #prob[N] <- prob[N-1]
      err[N] <- NaN
    }
    
    # Creamos data frame para regresar la curva de supervivencia
    res <- data.frame("tiempo" = unicos, "prob" = prob, "std" = sqrt((prob)^2 *err), "err" = err )
  }
  
  
  
  return(res)
}

# Funcion que determina la probabiliad de supervivencia basado en KM de un T0
KMT0 = function(KM, t0){
  # Caso especial, T0 se encuentra antes del primer tiempo (prob = 1)
  if(sum(KM$tiempo <= t0) == 0) {
    res <- 1
  }else{
    # Buscamos donde se encuentra T0 en la Curva de Supervivencia
    indice <- max(which(KM$tiempo <= t0))
    res <- KM$prob[indice]
  }
  return(res)
}

# Funcion que determina la probabiliad de supervivencia basado en KM de un T0 considerando 
# intervalos de confianza
KMT0CI = function(KM, t0, alpha){
  # Caso especial, T0 se encuentra antes del primer tiempo (prob = 1)
  if(sum(KM$tiempo <= t0) == 0) {
    prob <- 1
    err <- 0
  }else{
    # Buscamos donde se encuentra T0 en la Curva de Supervivencia
    indice <- max(which(KM$tiempo <= t0))
    prob <- KM$prob[indice]
    std <- KM$std[indice]
    err <- KM$err[indice]
  }
  
  #Intevalos Normales
  cDown <- prob - qnorm(1 - alpha/2) * sqrt( (prob)^2 *err )
  cUp <- prob + qnorm(1 - alpha/2) * sqrt( (prob)^2 *err )
  
  #Intervalos LogExp
  cDownLog <- log(-log(prob)) + qnorm(1 - alpha/2) * sqrt( 1/( log(prob) )^2 *err )
  cUpLog <- log(-log(prob)) - qnorm(1 - alpha/2) * sqrt( 1/( log(prob) )^2 *err )
  
  res <- data.frame(prob = prob,
                    down = cDown,
                    up = cUp,
                    downLog = exp(-exp(cDownLog)),
                    upLog = exp(-exp(cUpLog))
  )
  return(res)
}


# --------------------------------------------
# Estimador de Guessoum-Saïd
# --------------------------------------------

# Guessoum
guessoum = function(x0, X, t, d, h) {
  
  # Generamos los pesos de Nadaraya Watson
  W <- W(x0, X, h)
  
  # Generamos Kaplan Meier
  KMEst <- KM(t,1-d)
  
  # Aplicamos Kaplan Meier para toda t
  G <- sapply(t,function(T) {
    KMT0(KMEst, T)
  })
  
  # Juntamos los algoritmos
  r <- sum(W*(d * t)/G)
  return(r)
}

# Fn
fn = function(x0, X, hn){
  
  # Obtenemos el numero de observaciones
  n <- nrow(X)
  
  # Generamos el Kernel
  U <- -sweep(X, 2, x0, "-") / hn
  K <- gausKernel(U)
  
  # Juntamos resultados
  res <- sum(K) / (n*hn)
  return(res)
}

# R1
r1 = function(x0, X, t, d, hn){
  
  # Obtenemos el numero de observaciones
  n <- nrow(X)
  
  # Generamos el Kernel
  U <- -sweep(X, 2, x0, "-") / hn
  K <- gausKernel(U)
  
  # Generamos Kaplan Meier
  KMEst <- KM(t,1-d)
  
  # Aplicamos Kaplan Meier para toda t
  G <- sapply(t,function(T) {
    KMT0(KMEst, T)
  })
  
  # Juntamos los algoritmos
  r <- sum(K*(d * t)/G)
  res <- r / (n*hn)
  
  return(res)
}

# R2
r2 = function(x0, X, t, d, hn){
  
  # Obtenemos el numero de observaciones
  n <- nrow(X)
  
  # Generamos el Kernel
  U <- -sweep(X, 2, x0, "-") / hn
  K <- gausKernel(U)
  
  # Generamos Kaplan Meier
  KMEst <- KM(t,1-d)
  
  # Aplicamos Kaplan Meier para toda t
  G <- sapply(t,function(T) {
    KMT0(KMEst, T)
  })
  
  # Juntamos los algoritmos
  r <- sum(K*(d * (t^2) )/G)
  res <- r / (n*hn)
  
  return(res)
}

# Sigma Plug-in
sig2 = function(x0, X, t, d, hn){
  
  # Generamos los estimadores de fn, r1 y r2
  fnEst <- fn(x0, X, hn)
  r1Est <- r1(x0, X, t, d, hn)
  r2Est <- r2(x0, X, t, d, hn)
  
  # Obtenemos las dimensiones con las que estamos trabajando para kappa
  dim <- ncol(X)
  kappa <- (4 * pi)^(-dim / 2)
  
  # Juntamos los estimadores y kappa
  res <- kappa * ( (r2Est * fnEst^2 - r1Est^2 * fnEst) / fnEst^4 )
  
  return(res)
}

# Guessoum CI
guessoumCI = function(x0, X, t, d, hn, alpha){
  
  # Obtenemos el numero de observaciones
  n <- nrow(X)
  
  # Generamos el estimador de Guessoum
  r <- guessoum(x0, X, t, d, hn)
  
  # Generamos el estimador de sigma plug-in
  sigma2 <- sig2(x0, X, t, d, hn)
  
  # Generamos el quantil de la distribución normal
  z <- qnorm(1 - alpha/2) * sqrt(sigma2 /(n*hn))
  
  # Juntamos los estimadores
  resIzq <- r - z
  resDer <- r + z
  
  res <- data.frame(Est = r,
                    Izq = resIzq,
                    Der = resDer
  )
  
  return(res)
}

# Guessoum extendido
guessoumExt = function(x0Cont, x0Cat, XCont, XCat, t, d, h) {
  
  # Generamos los pesos de Nadaraya Watson
  # Preguntamos por casos: Primero si tenemos variables continuas y categoricas
  if( ( length(XCont) != 1 ) && ( length(XCat) != 1 ) ) {
    We <- WMix(x0Cont, x0Cat, XCont, XCat, h)
  }else{
    if( length(XCont) != 1 ){
      # Si solo tenemos continuas
      We <- W(x0Cont, XCont, h)
      
    }else{
      # Si solo tenemos categoricas
      We <- WCat(x0Cat, XCat)
      
    }
  }
  
  # Generamos Kaplan Meier
  KMEst <- KM(t,1-d)
  
  # Aplicamos Kaplan Meier para toda t
  G <- sapply(t,function(T) {
    KMT0(KMEst, T)
  })
  
  # Juntamos los algoritmos
  r <- sum(We*(d * t)/G)
  return(r)
}

# Error
miseIPCW = function(r, t, d) {
  
  # Tamaño de la muestra
  n <- length(t)
  
  # Calculamos las diferencias
  dif <- (r - t)^2
  
  # Generamos Kaplan Meier
  KMEst <- KM(t,1-d)
  
  # Aplicamos Kaplan Meier para toda t
  G <- sapply(t,function(T) {
    KMT0(KMEst, T)
  })
  
  # Corregimos el último valor de supervivencia para no dividir entre 0
  G[which(G == 0)] = G[which(G == 0)-1]
  
  # Juntamos los algoritmos
  r <- 1/n * sum((d * dif)/G)
  
  return(r)
}

# Función que genera estimación
aplicacionGuessoum = function(n, funcion, b){
  
  set.seed(1)
  #Generamos los datos
  X <- matrix(runif(n, min = 0, max = 1))
  e <- rnorm(n)
  Y <- funcion(X) + b * e
  C <- rnorm(n)
  T <- pmin(Y, C)
  delta <- as.numeric(Y <= C)
  h <- h.ucv(X[,1], kernel = "gaussian")$h
  
  # Definimos eje
  x0 <- seq(0, 1, length.out = 500)
  
  # Generamos el estimador
  est = sapply(x0, function(x0) {
    guessoum(x0, X, T, delta, h)
  })
  
  # Imprimimos el error
  cat("Error: ", miseIPCW(est, funcion(x0), delta))
  
  # Guardamos los datos
  res <- data.frame(x = x0,
                    y = funcion(x0),
                    est = est
  )
  return(res)
}

# Guessoum CI Bootstrap
bootstrapCI = function(x0, X, t, d, h, B){
  
  # Tamaño de los datos
  n <- nrow(X)
  
  # Create the vector where to save the estimators
  rB <- numeric(B)
  
  # Iteramos B veces
  for(i in 1:B){
    # Sampleamos con reemplazo
    ind <- sample(1:n, size = n, replace = TRUE)
    # Obtenemos los valores
    XB <- X[ind, ]
    tB <- t[ind, ]
    dB <- d[ind]
    # Calculamos guessoum considerando XB
    rB[i] <- guessoum(x0, matrix(XB), tB, dB, h)
    
  }
  
  # Quitamos posbiles NaNs generados por dividir entre 0
  rB <- na.omit(rB)
  
  # Obtenemos el error estandar "puro"
  err <- sd(rB)
  
  # Obtenemos guessoum para el punto x0
  r <- guessoum(x0, X, t, d, h)
  
  # Calculamos intervalos de confianza
  CI <- quantile(rB, c(0.025, 0.975))
  
  res <- data.frame("Low" = CI[1],
                    "Est" = r,
                    "Up" = CI[2]
  )
  
  return(res)
}

# Función que genera estimación considerando CI Bootstrap
aplicacionGuessoumBoot = function(n, funcion, b){
  
  set.seed(1)
  #Generamos los datos
  X <- matrix(runif(n, min = 0, max = 1))
  e <- rnorm(n)
  Y <- funcion(X) + b * e
  C <- rnorm(n)
  T <- pmin(Y, C)
  delta <- as.numeric(Y <= C)
  h <- h.ucv(X[,1], kernel = "gaussian")$h
  
  # Definimos eje
  x0 <- seq(0, 1, length.out = 500)
  
  # Generamos los intervalos
  IC <- lapply(x0, function(x0) {
    bootstrapCI(x0, X, T, delta, h, 30)
  })
  
  # Guardamos los datos
  res <- do.call(rbind, IC)
  res$x <- x0
  res$y <- funcion(x0)
  
  return(res)
}


