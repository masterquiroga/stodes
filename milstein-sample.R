#! /usr/bin/Rscript
##' Ejemplos de uso para el solucionador de ecuaciones diferenciales estocásticas
##' implementando el Método de Milstein.
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
source("milstein.R")

# Semilla inicial
set.seed(0)


# SCIENTIFIC AESTHETICS
#require(grDevices)
#require(fontcm)
#require(Cairo)
#X11(family = "Times")
#par(family = "ComputerModern")



##' -------------------------------------------------------------------------

# Puente Browniano
# dXt = -Xt/(1-t) + dBt

png(
  filename = "brownianbridge.png",
  width = 600,
  height = 600,
  type = "cairo-png")
milstein(
  b = function(x, t){-x/(1-t)}, 
  s = function(x, t){1}, s1 = function(x,t){0},
  X0 = 0, Tmax = 1, 
  n = 10000,
  M = 10,
  plotting = T,
  axes = T,
  bty = "n",
  font.main = 1,
  col.main = "black",
  col.lab = "gray",
  col.axis = "gray",
  ylab = "X(t)",
  #ylab = latex2exp::TeX("$X_t$"),
  xlab = "t")
dev.off()



##' -------------------------------------------------------------------------

## Ejemplo para una función periódica

# Número de simulaciones
M <- 10 # usar idealmente un número grande como 500...
# Número de periódos
k <- 5
# Escala de precisión
# 10000 ""revoluciones"" por periódo
n <- k*ceiling(10000/(2*pi))
# Graficación
png(
  filename = "periodic1.png",
  width = 600,
  height = 600, type = "cairo-png")
X <- milstein(
  X0 = 0, 
  b = function(x,t){t*sin(t*x)}, 
  s = function(x,t){sin(t)}, 
  s1 = function(x,t){0},
  Tmax = k*2*pi,
  M = M, n = n,
  col.mean = gray(level = 0, alpha = 0.38),
  lwd.path = 0.1,
  axes = T,
  bty = "n",
  main = latex2exp::TeX(
    string = "$dX_t = t \\sin(tX_t)dt + \\sin(t)dB_t , X_0 = 0$"), 
  font.main = 1,
  col.main = "black",
  col.lab = "gray",
  col.axis = "gray",
  ylab = "X(t)",
  #ylab = latex2exp::TeX("$X_t$"),
  xlab = "t")
dev.off()

##' -------------------------------------------------------------------------

## Ejemplo para otra función periódica

# Número de simulaciones
M <- 10
# Número de periódos
k <- 5
# Escala de precisión
# 10000 revoluciones por periódo
n <- k*ceiling(10000/(2*pi))
# Condición inicial
# Si no no es Lipschitz continua
X0 <- -2*pi/k
# Graficación
png(
  filename = "periodic2.png",
  width = 600,
  height = 600)
milstein(
  b = function(x,t){1*cos(t*x)*(log(t + 1) - sin(t*x))},
  s = function(x,t){-sin(t*x) - cos(min(t,x))},
  s1 = function(x,t){-t*cos(t*x) + sin(x)*(x >= t)}, Tmax = (2*k*pi),
  n = k*ceiling(10000/(2*pi)),
  M = M, 
  X0 = X0,
  col.mean = gray(0,0.38), 
  axes = T,
  bty = "n",
  main = "Periodic function 1",
  #main = latex2exp::TeX(
  #  string = "$dX_t = \\cos(tX_t) (\\log(t + 1) - sin (t X_t)) dt 
  #                   - (\\sin(t X_t) + \\cos(t \\wedge X_t))dB_t$"), 
  font.main = 1,
  #cex.main = 0.8,
  ylab = "X(t)",
  #ylab = latex2exp::TeX("$X_t$"),
  xlab = "t",
  col.axis = "gray",
  col.lab = "gray"
)
text(
  x = 0, 
  y = X0, 
  labels = "X(0) = -2 pi / 5",
  #labels = latex2exp::TeX(paste("$X_0 =", -2, "\\pi / ", k, "$")),
  cex = 0.5)
dev.off()
