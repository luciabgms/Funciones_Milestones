import numpy as np
from numpy import array, concatenate, zeros, pi, linspace, matmul
from numpy.linalg import norm 
import matplotlib.pyplot as plt


##### Funcion F_kepler
def F_kepler(U,t):
    r=array(U[:2])
    drdt=array(U[2:])
    F = concatenate((drdt,-r/np.linalg.norm(r)**3))
    return F

##### Esquemas numericos = EN
def EN_Euler( F, t1, t2, U1):
    
    dt = t2 - t1
    t=t1
    U2 = U1 + dt * F( U1,t)
    return U2

#def EN_Inverse_Euler
#def EN_Crank_Nicolson

def EN_Runge_Kutta4 ( F, t1, t2, U1):
    dt = t2 - t1 
    t = t1

    k1 = F(U1, t)
    k2 = F( U1 + dt * k1/2, t + dt/2)
    k3 = F( U1 + dt * k2/2, t + dt/2)
    k4 = F( U1 + dt * k3,   t + dt )

    U2 = U1 + dt * ( k1 + 2*k2 + 2*k3 + k4 )/6

    return U2
#def EN_LeapFrog

#### Problema de Cauchy
def CauchyProblemS( time_domain, differential_operator, solution, scheme):
    N_steps = len(time_domain) - 1
    for i in range(N_steps):
        t1 = time_domain[i]
        t2 = time_domain[i + 1]
        solution[i+1,:] =  scheme( differential_operator, t1, t2, solution[i,:])
    return solution

#### Orbita de Kepler
def Kepler_orbit( n_vueltas, N, scheme):
    print("orbita de kepler")
    #quiero que mi orbita de 6 vueltas
    t0 = 0
    tf = n_vueltas*2*pi
   
    U = zeros( (N,4), dtype=float)
    #genero un array que empieza en t0 y da 6 vueltas en 10000 puntos
    Time = array( [t0+(tf-t0)*i/N for i in range(N)] )
    ### Time = np.arange(t0,tf,en incrementos de t) 

    #esto es la condicion inicial
    cond_kepler = [1. , 0. , 0. , 1.]
    U[0] = array(cond_kepler)
    CauchyProblemS(Time,F_kepler,U, scheme)
   

    plt.figure(0)
    plt.plot(Time, U[:,0])
    plt.plot(Time, U[:,1])
    plt.xlabel('Tiempo')
    plt.ylabel('Distancia')
    plt.title("Kepler orbit")

    plt.figure(1)
    plt.plot(U[:,0], U[:,1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Diagrama de fases")
    plt.axis('square')

    plt.show()

#### ERRORES 
def Error_Cauchy_Problem( time_domain, Differential_operator, initial_condition, Scheme, order ):

    t1 = time_domain
    N = len(time_domain)

    t2_len = 2*N - 1
    t2 = zeros( t2_len, dtype=float)
    for i in range (0, N-1):
        t2[2*i] = t1[i]
        t2[2*i+1] = ( t1[i] + t1[i+1])/2
    t2[t2_len - 1] = t1[N-1]

    U1 = zeros((N,4), dtype=float)
    U1[0] = array(initial_condition)
    CauchyProblemS(t1, Differential_operator, U1, Scheme)

    U2 = zeros((2*N-1,4), dtype=float)
    U2[0] = array(initial_condition)
    CauchyProblemS(t2, Differential_operator, U2, Scheme)

    Error = zeros((N,4), dtype=float)
    for i in range (0,N-1):
        Error[i,:] = ( U1[i,:] - U2[2*i,:] )/( 1 - 0.5**order)

    return Error

def Error_extrapolacion_richardson( differential_operator, cond_inicial, scheme, order, time_domain ):

    N_intervals = len(time_domain) - 1
    N_points = len(time_domain)

    t1 = time_domain
    t2 = zeros( [2*N_intervals + 1] )

    for i in range (0, N_points):
        t2[2*i] = t1[i]
    for i in range (0, N_points-1):
        t2[2*i + 1] = ( t1[i] + t1[i+1] )/2
       
    U1 = zeros( [N_intervals + 1,   4] )
    U1[0,:] = array( cond_inicial)
    CauchyProblemS( t1, differential_operator, U1, scheme)

    U2 = zeros( [2*N_intervals + 1, 4])
    U2[0,:] = array( cond_inicial)
    CauchyProblemS( t2, differential_operator, U2, scheme)

    Error = zeros( [N_intervals + 1, 4] )
    for i in range(0, N_intervals + 1):
       Error[i] = ( U1[i] - U2[2*i] )/(1 - 0.5**order) 
    
       error_final = norm(Error[N_intervals])
    return error_final

def Mapa_convergencia_error(t0, tf, differential_operator, cond_inicial, scheme, order):
    
    N = 10
    Error_final = zeros([N,1])

    for i in range(1, N, 1): 
        time = linspace(t0, tf, i+1)
        
        Error_final[i,1] = Error_extrapolacion_richardson( differential_operator, cond_inicial, scheme, order, time )
       
    return Error_final
       
    plt.figure(1)
    plt.plot(norm_Error[N], N)
    plt.xlabel('logN')
    plt.ylabel('logE')
    plt.title("Mapa de convergencia del error")
    plt.show()

#### Oscilador armonico 
def Osc_armonico(U,t):

    #la ecuacion es x''+ x = 0
    A = array([[0, 1], [-1, 0]])
    F = matmul(A,U)

    return F 

def Sol_osc_armonico(N, time, cond_inicial, scheme):

    U = zeros([N, 2])
    U[0] = array(cond_inicial) 

    Osc_armonico_Solution = CauchyProblemS( time, Osc_armonico, U, scheme)

    plt.figure(0)
    plt.plot(time, U[:,0])
    plt.xlabel('Tiempo')
    plt.ylabel('Posicion')
    plt.title("Posicion del oscilador")

    plt.figure(1)
    plt.plot(time, U[:,1])
    plt.xlabel('Tiempo')
    plt.ylabel('Velocidad')
    plt.title("Velocidad del oscilador")

    plt.show()