#import numpy as np
from numpy import zeros, array, abs, linspace
import matplotlib.pyplot as plt
### paradigma imperativo
def factorial_I(n):
    f = 1
    for i in range(1,n+1):
        f=f*i
    return f 

##print("factorial(4)=",factorial_I(4))

### paradigma funcional (se sustituye iteracion por recursion) 
def factorial_F(n):
    if n==0 : return 1
    else : return n*factorial_F(n-1)

##print("factorial(4)=",factorial_F(4))

### MAPING 
#def f(x):
#    return x**2

#x = np.array([1, 2, 3])
#y = f(x)
##print(y)

### FILTER 
#z = x[x>1]
##print(z)

### REDUCE 
##m = np.norm(x)
##print(m)

#s = np.sum(x)
#print(s)

#p = np.product(x)
#print(p)

#### JACOBIANO
def system_matrix( F, U0):
    N = len(U0)
    A = array( zeros([N,N]))
    delta = array(zeros([N]))
    eps = 1e-14 
    t=0
    
    for j in range(N):
            delta[:] = 0
            delta[j] = eps
            A[:,j] = (F(U0 + delta, t) - F(U0 - delta, t))/(2*eps)
    return A

def oscillator(U,t):
    return array([U[1], -U[0]])

def test_system_matrix():
    U0 = array([0., 0.])
    print("A = ", system_matrix(oscillator, U0))

#test_system_matrix() 

### REGION DE ESTABILIDAD PARA ESQUEMAS NUMERICOS DE 1 PASO 

def Euler(U1, t, dt, F):

    return U1 + dt * F( U1,t)

def stability_region( scheme, x, y):
    N = len(x)
    M = len(y)
    rho = array(zeros([N,M]))
    dt = 1.
    t = 0
    for i in range(N):
        for j in range(M):
            w = complex(x[i], y[j])
            U1 = complex(1., 0)
            r = scheme(U1, t, dt, lambda U, t : w*U)
            rho[i,j ] = abs(r)

    return rho

def test_stability_region():

    x = linspace(-2, +2, 100)
    y = linspace(-2, +2, 100)
    rho = stability_region(Euler,x, y)
    plt.contour(x, y, rho, linspace(0., 1., 11))
    plt.show() 

#test_stability_region()