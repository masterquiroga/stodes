#! /usr/bin/Rscript
##' Proporciona una aproximación numérica de la trayectoria de una difusión que
##' soluciona una ecuación diferencial estocástica implementando el Método de
##' Milstein.
##' 
##' @author ＭＡＳＴＥＲキロガ <masterquiroga@protonmail.ch>
##' @license MIT License
# 
# Copyright (c) 2017- Víctor Gerardo González Quiroga
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##' 
##' @param b function el coeficiente de deriva del proceso
##' @param s function el coeficiente de difusión del proceso
##' @param s1 function el coeficiente  de tasa cambio de la difusión del proceso
##' @param X0 double la condición inicial del proceso
##' @param Tmax double el último tiempo a observar
##' @param n integer número de subintervalos
##' @param M integer número máximo de simulaciones
##' @param plotting boolean ¿se debe graficar o no?
##' @param col.mean character el color para el valor esperado
##' @return numeric matriz con las trayectorias simuladas (por columna)
##' @references Nualart D., "Lecture Notes on Stochastic Processes"
##' @usage 
##'   # Puente Browniano (T = 1)
##'   milstein(
##'     b = function(x, t){-x/(1-t)}, 
##'     s = function(x, t){1}, s1 = function(x,t){0}
##'     X0 = 0, Tmax = 1, 
##'     n = 1000, M = 500, plotting = T
##'   )
milstein <- function(
  b = function(x,t){1},
  s = function(x,t){1},
  s1 = function(x,t){0},
  X0 = 0,
  Tmax = 1,
  n = 1000,
  M = 500,
  plotting = T,
  col.mean = "black",
  col.path = rainbow(n = M, alpha = max(0.02, 1/M )),
  lwd.path = 0.5,
  prog.bar = T,
  ...
){
  i = 1:n            # Índices en partición
  t <- c(0,i)*Tmax/n # Partición del intervalo
  delta <- Tmax/n    # Longitud de cambio
  Xt <- NULL         # Procesos simulados
  ss1 <- function(x,t){s(x,t)*s1(x,t)} # Función auxiliar

  cat("Simulando",M,"trayectorias...\n")
  if (prog.bar) {
    pb <- progress::progress_bar$new(
      total = M*1,
      format = "(:spin) [:bar] :percent | ETA :eta | T :elapsed",
      clear = F)
  } else {
    pb <- c(tick = function(x){})
  }
  pb$tick(0)
  
  for (sim in 1:M) {
    DeltaB <- sqrt(delta)*rnorm(n, 0, 1)
    X <- X0
    for (k in i) {
      X <- c(X, 
        X[k] + b(X[k], t[k])*delta + s(X[k], t[k])*DeltaB[k] #Aproximación
        + (ss1(X[k], t[k]))*((DeltaB[k]^2 ) - delta)/2 # Corrección de Itô
      )
    }
    Xt <- cbind(Xt, X)
    pb$tick()
  }
  
  if (plotting) {
    cat("Graficando",M,"trayectorias...\n")
    if (prog.bar) {
      pb <- progress::progress_bar$new(
        total = M*1,
        format = "(:spin) [:bar] :percent | ETA :eta | T :elapsed",
        clear = F)
    } else {
      pb <- c(tick = function(x){})
    }
    pb$tick(0)
    
    plot(
      x = t,
      y = Xt[, 1],
      type = "l",
      ylim = range(Xt),
      col = col.path[1], 
      lwd = lwd.path,
      ...)
    pb$tick()
    
   for (sim in 2:M) {
      lines(
        x = t,
        y = Xt[, sim],
        col = col.path[sim],
        lwd = lwd.path)
      pb$tick()     
    }
    lines(
      x = t,
      y = rowMeans(Xt),
      col = col.mean)
  } 
  return(Xt)
}
