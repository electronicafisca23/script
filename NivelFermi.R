#Descripción:  Cálculo del nivel de Fermi de una superred GaAs / AlAs.
#Año:          2023
#Bibliografía: Luscombe, J. H., Aggarwal, R., Reed, M. A., Frensley,
#              W. R., &amp; Luban, M. (1991). Theory of the fermi-level
#              energy in semiconductor superlattices. Physical Review B,
#              44(11), 5873–5876. https://doi.org/10.1103/physrevb.44.5873 

#Empezamos describiendo las constantes que necesitamos (en SI)
temp = 300
kB = 1.3806e-23
beta = 1/(kB*temp)
hbar = 1.0546e-34

#Y los parámetros que describen el sistema
n = 2e24 #en m^-3
m0 = 9.10938e-31
m_eff_GaAs = 0.067*m0
m_eff_AlAs = 0.21*m0
d_GaAs = 45.3e-10
d_AlAs = 17e-10
a = d_AlAs

#Calculamos la masa efectiva como se describe en el paper
delta_GaAs = d_GaAs / (d_GaAs + d_AlAs)
delta_AlAs = d_AlAs / (d_GaAs + d_AlAs)
m_eff = 1/(delta_AlAs/m_eff_AlAs + delta_GaAs/m_eff_GaAs)

#Introducimos manualmente los valores de posición (e0) y medio
#ancho de banda (W) para cada minibanda:
#e0 = c(0.121, 0.461, 0.950, 1.523, 2.229)*1.6022e-19
#W = (c(0.020, 0.099, 0.290, 0.630, 0.959)*1.6022e-19)/2
e0 = c(0.121, 0.461, 0.950)*1.6022e-19
W = (c(0.020, 0.099, 0.290)*1.6022e-19)/2

#Creamos un vector de valores para la E_F (afinar posteriormente)
EF = seq(0,0.3,len=2000)*1.6022e-19

#Definimos la función a integrar
intFermi = function(x,EF,e0,W) 1/(2*pi)*log(1 + exp(beta*(EF-e0+W*cos(x))))
n_est = function(EF) {
  ans = 0
  for(j in 1:length(e0)) {
    ans = ans + integrate(intFermi, lower=-pi, upper=+pi,
                          EF = EF, e0 = e0[j], W=W[j])$value
  }
  ans = m_eff * ans / (pi*a*beta*hbar^2)
  return(ans)
}

#Corremos la prueba para cada valor de EF, seleccionando el que menor
#error cometa
result = Inf
ref = 0
for(k in 1:length(EF)) {
 if(abs(n - n_est(EF[k]))<result) {
   result = abs(n - n_est(EF[k]))
   ref = EF[k]
   print(abs(n - n_est(EF[k])))
 }
}

#Nuestra estimación para la EF (expresada en eV sobre el pozo del potencial):
EF_estimada = ref/(1.602e-19)

#Calculamos el error relativo cometido por nuestra EF estimada: (n - n_estimado)/n*100
error = result/n*100
