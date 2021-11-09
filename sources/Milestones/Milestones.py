from Funciones_Milestones import *
from numpy import linspace, array, pi 
import matplotlib.pyplot as plt


#Kepler_orbit(6, 10000, EN_Runge_Kutta4)

Error = Error_extrapolacion_richardson( F_kepler, [1. , 0. , 0. , 1.], EN_Euler, 1, linspace(0,12*pi,10000) )
#print("Error de la funcion de Kepler =", Error)

#Mapa_convergencia_error( 0, 12*pi, F_kepler, [1. , 0. , 0. , 1.], EN_Euler, 1)

##### Oscilador Armonico 
#resuelvo mi problema de cauchy
N = 5000
time = linspace(0, 16*pi, N)
cond_inicial = array([1, 0])

Sol_osc_armonico(N, time, cond_inicial, EN_Euler)
