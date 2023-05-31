# /home/jorge/anaconda3/bin/python
# -+- coding: utf-8 -+-

# Hacemos los cálculos y las gráficas de la práctica de Electronica Fisica

# Importamos el módulo de pylab para poder trabajar con libertad, también pondremos un print vacío con una función meramente estética.

from pylab import *

print(" ")

# Activamos el modo interactivo para visualizar las gráficas:

ion()

# Definimos la lista de datos de los voltajes:

volt = linspace(0.5,0.7,10) # [V]

# Definimos los parametros del pozo de potencial en Amstrong:

a = 10 # [A]
W = 70 # [A]

# Pasamos W a metros:

A = W*1E-10 # m

# Establecemos el resto de parámetros:

C = 50
q = -1.6E-19 # [C]
m0 = 9.109E-31 # [kg]
mef = 0.00119*m0 # [kg]
mefp = -0.00159*m0 # [kg]
k = 8.617332385E-5 # [eV/K]
T = 300 # [K]
Tsc = 10 # [K]

# Tomamos la energía mínima de la banda de conducción:

Ecmax = 0.191 # [eV]
Ecmin = 0.122 # [eV]

# Establecemos el nivel de Fermi y la distancia del campo

Ef = 0.1 # [eV]

d = 0.10 # [m]

# Calculamos el campo eléctrico:

E = volt/d # [V/m]

# Calculamos el rango de velocidades para n y p:

v = sqrt((-A*C*q*E)/mef) # [m/s]
vp = sqrt((A*C*q*E)/mefp) # [m/s]

# Calculamos el tiempo de relajación para n y p:

tn = C*A/v # [s]
tp = C*A/vp # [s]

# Calculamos la integral para obtener el número de electrones:

# Como en este caso n=p solo calculamos la integral una vez:

# Leemos los datos:

DOS, Egrid = loadtxt("1DDOSplot_d70.txt", unpack=True, skiprows=5, max_rows=2865, delimiter=",") 

ge = DOS*2.49424e+09 # #/eV/m

# Definimos la función de Fermi:

def f(x):
	return 1/(1+exp((x-Ef)/(k*T)))

# filtramos los datos para tomar los que están en la banda de conducción:
	
Egridint = [x for x in Egrid if Ecmin<= x <=Ecmax]
	
geint = ge[(Egrid >= Ecmin) & (Egrid<=Ecmax)] 
	
# Calculamos la integral:

s = 0.0

for i in range(len(Egridint)):
	if i==0:
		s += f(Egridint[i])*geint[i]*(Egridint[i]-Ecmin)
	else:
		s += f(Egridint[i])*geint[i]*(Egridint[i]-Egridint[i-1])

n = s

# Finalmente calculamos la j:

j = (tn/mef+tp/mefp)*(n*q**2*E)

# Guardamos los datos obtenidos:

# Creamos el contenido:

contenido = f"V [mV], "+"j [A]"+"\n"

for i in range(len(volt)):
	contenido += "%.8f"%(volt[i]*1000)+", "+"%s"%j[i]+"\n"

# Guardamos:

with open("valores.txt","w") as archivo:
	archivo.writelines(contenido)

# Graficamos los valores:

figure(1)

plot(volt,j,label="Simulado")

# Añadimos los ejes:

xlabel("Voltaje [mV]")
ylabel("Densidad de corriente [A]")
legend(loc="best")
grid()

# Guardamos la figura:

savefig("grafica.png")

###############################################################################################################################################################################################
###############################################################################################################################################################################################


# Hacemos lo mismo para la parte de superconductor:

def fsc(x):
	return 1/(1+exp((x-Ef)/(k*Tsc)))

s = 0.0

for i in range(len(Egridint)):
	if i==0:
		s += f(Egridint[i])*geint[i]*(Egridint[i]-Ecmin)
	else:
		s += f(Egridint[i])*geint[i]*(Egridint[i]-Egridint[i-1])

nsc = s

vsc = sqrt((-2*q*volt)/m0)

jsc = nsc*vsc*(-q)

# Guardamos los datos en otro archivo:

# Creamos el contenido:

contenido2 = f"V [mV], "+"j [A]"+"\n"

for i in range(len(volt)):
	contenido2 += "%.8f"%(volt[i]*1000)+", "+"%s"%jsc[i]+"\n"

# Guardamos

with open("valores_superconductor.txt","w") as archivo2:
	archivo2.writelines(contenido2)

# Graficamos los valores obtenidos:

figure(2)

plot(volt,jsc, label="Superconductor")

# Añadimos los ejes:

xlabel("Voltaje [mV]")
ylabel("Densidad de corriente [A]")
legend(loc="best")
grid()

# Guardamos la gráfica:

savefig("grafica_superconductor.png")
