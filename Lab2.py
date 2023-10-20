import numpy as np
import matplotlib.pyplot as pl
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")
# КОНСТАНТЫ -------------------------------------------------------------------
a0 = 1
b0 = -1
a1 = 1
b1 = -1
k0 = 1 
#------------------------------------------------------------------------------

#ФУНКЦИИ ДЛЯ РЕШЕНИЯ ЗАДАЧИ ---------------------------------------------------
def f(x, t):
    return pow(x,3) - 6*t*x + 2
    #return x ** 2 - 2 * t

def phi(x):
    return x*(1-x)
    #return x

def psi_0(t):
    return -1

def psi_1(t):
    return 1-2 * t

def k(u):
    return 1
    #return u ** 2 + u

def F(u):
    return 1
    #return u ** 2

def f_Cons(x, t, u):
    return F(u)*f(x,t)

def Solution(x, t):
    return t*pow(x,3)+x*(1-x)
    #return t * x ** 2 + x

def CalculateError(u, u_sol):
    maxdelta = -1
    for i in range(L+1):
        for j in range(N+1):
            delta = np.abs(u[i,j]-u_sol[i,j])
            if delta > maxdelta:
                maxdelta = delta
    return maxdelta

#------------------------------------------------------------------------------

#ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ------------------------------------------------------

def solveMatrix(a,b,c,f,N):
    x = np.empty(N+1)
    
    for i in range(1,N+1):
        m = a[i]/b[i-1]
        b[i] = b[i]-m*c[i-1]
        f[i] = f[i]-m*f[i-1]
        
    x[N] = f[N]/b[N]
    
    for i in range(N-1, 0,-1):
        x[i] = (f[i]-c[i]*x[i+1])/b[i]
        
    return x

def PlotError(x_step, t_step, error):
    fig = pl.figure()

    X, t = np.meshgrid(x_step, t_step)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, t, error)

    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')

    pl.show()
    
#------------------------------------------------------------------------------

T = 1.
N = 40
L = int(2 * k0 * N ** 2 * T)
h = 1. / N 
tau = T / L

x_step = np.empty(N+1)
t_step = np.empty(L+1)

for i in range(N+1):
    x_step[i] = i * h
for i in range(L+1):
    t_step[i] = i * tau

u_sol = np.empty((L+1,N+1), dtype="float")
for i in range(L+1):
    for j in range(N+1):
        u_sol[i,j] = Solution(j * h, i * tau)
        
print("N = ", N, ", L = ", L, ", tau = ", tau, ", h = ", h, "\n")


# =============================================================================
# Центрально-разностная явная схема -------------------------------------------
start_time = datetime.now()

#Инициализация массива для решения
u = np.empty((L+1,N+1), dtype="float")
   
for i in range(N+1):
    u[0,i] = phi(i*h)
    
    
#Алгоритм
for n in range(L):
    for i in range(1,N):
        u[n+1,i] = (k0*tau / pow(h, 2)) * u[n][i + 1] + (1 - (2 * k0 * tau / pow(h, 2))) * u[n][i] + (k0 * tau / pow(h, 2)) * u[n][i - 1] + tau * f(i * h, tau * (n + 1))
    
    u[n+1,0] = (h * psi_0(tau * (n + 1)) - b0 * u[n + 1,1]) / (a0 * h - b0)
    u[n+1,N] = (h * psi_1(tau * (n + 1)) + b1 * u[n + 1,N - 1]) / (a1 * h + b1)

print("Время явной схемы: ", datetime.now() - start_time, "\n")


#Вычисление максимальной погрешности и отрисовка погрешности
error = np.empty((L+1,N+1), dtype="float")
for i in range(L+1):
    for j in range(N+1):
        error[i,j] = np.abs(u_sol[i,j]-u[i,j])
        
PlotError(x_step, t_step, error)     
print("KRS Error: ", CalculateError(u, u_sol), "\n")
# =============================================================================




# =============================================================================
# НЕЯВНАЯ Центрально-разностная схема -----------------------------------------
start_time = datetime.now()

#Инициализация массива для решения
u = np.empty((L+1,N+1), dtype="float")
    
for i in range(N+1):
    u[0,i] = phi(i*h)


#Алгоритм
for n in range(L):
    A = np.empty(N+1)
    B = np.empty(N+1)
    C = np.empty(N+1)
    D = np.empty(N+1)
    
    B[0] = a0 * h - b0
    C[0] = b0
    D[0] = h * psi_0(tau * (n + 1))
    
    A[N] = -b1
    B[N] = a1 * h + b1
    D[N] = h * psi_1(tau * (n + 1))
    
    for i in range(1, N):
    	A[i] = k0 * tau / pow(h,2)
    	B[i] = -1 - (k0 * 2 * tau / pow(h, 2))
    	C[i] = k0 * tau / pow(h, 2)
    	D[i] = -u[n,i] - tau * f(i * h, tau * (n + 1))
    #u[n + 1] = solveMatrix(A, B, C, D, N)
    buf = solveMatrix(A, B, C, D, N)
    for i in range(len(buf)):
        u[n+1,i] = buf[i]
    
print("Время неявной схемы: ", datetime.now() - start_time, "\n")


#Вычисление максимальной погрешности и отрисовка погрешности
error = np.empty((L+1,N+1), dtype="float")
for i in range(L+1):
    for j in range(N+1):
        error[i,j] = np.abs(u_sol[i,j]-u[i,j])

PlotError(x_step, t_step, error)
print("Implicit KRS Error: ", CalculateError(u, u_sol), "\n")


# Схема Кранка – Николсона ----------------------------------------------------
start_time = datetime.now()

#Инициализация массива для решения
u = np.empty((L+1,N+1), dtype="float")
    
for i in range(N+1):
    u[0,i] = phi(i*h)


#Алгоритм
for n in range(L):
    A = np.empty(N+1)
    B = np.empty(N+1)
    C = np.empty(N+1)
    D = np.empty(N+1)
    
    B[0] = a0 * h - b0
    C[0] = b0
    D[0] = h * psi_0(tau * (n + 1))
    
    A[N] = -b1
    B[N] = a1 * h + b1
    D[N] = h * psi_1(tau * (n + 1))
    
    for i in range(1, N):
    	A[i] = k0 * tau / (2. * pow(h, 2))
    	B[i] = -1 - (k0 * tau / pow(h, 2))
    	C[i] = k0 * tau / (2.0* pow(h, 2))
    	D[i] = -u[n,i] - tau * f(i * h, tau * (n + 1)) - (k0 * tau / (2. * pow(h, 2)) * (u[n,i + 1] - 2. * u[n,i] + u[n,i - 1]))
    #u[n + 1] = solveMatrix(A, B, C, D, N)
    buf = solveMatrix(A, B, C, D, N)
    for i in range(len(buf)):
        u[n+1,i] = buf[i]
    
print("Время схемы Кранка-Николсона : ", datetime.now() - start_time, "\n")


#Вычисление максимальной погрешности и отрисовка погрешности
error = np.empty((L+1,N+1), dtype="float")
for i in range(L+1):
    for j in range(N+1):
        error[i,j] = np.abs(u_sol[i,j]-u[i,j])

PlotError(x_step, t_step, error)
print("Krank-Nicholson scheme Error: ", CalculateError(u, u_sol), "\n")
# =============================================================================


# КОНСЕРВАТИВНАЯ СХЕМА --------------------------------------------------------
start_time = datetime.now()

#Инициализация массива для решения
u = np.empty((L+1,N+1), dtype="float")
    
for i in range(N+1):
    u[0,i] = phi(i*h)


#Алгоритм
for n in range(L):
    A = np.empty(N+1)
    B = np.empty(N+1)
    C = np.empty(N+1)
    D = np.empty(N+1)
    
    B[0] = a0 * h - b0
    C[0] = b0
    D[0] = h * psi_0(tau * (n + 1))
    
    A[N] = -b1
    B[N] = a1 * h + b1
    D[N] = h * psi_1(tau * (n + 1))
    
    for i in range(1, N):
        v_summ = k((u[n,i] + u[n,i + 1]) / 2.)
        v_difference = k((u[n,i] + u[n,i - 1]) / 2.)
        A[i] = k0 * tau / (pow(h, 2)) * v_difference
       	B[i] = -1. - (k0 * tau / pow(h, 2)) * (v_summ + v_difference)
        C[i] = k0 * tau / (pow(h, 2)) * v_summ
        D[i] =  -u[n,i] - tau * (1 - k0) / h * (v_summ * (u[n,i + 1] - u[n,i]) / h - v_difference * (u[n,i] - u[n,i - 1]) / h) - tau * f_Cons(i * h, tau * (n + 1), u[n,i + 1])
    u[n + 1] = solveMatrix(A, B, C, D, N)

print("Время консервативная схема : ", datetime.now() - start_time, "\n")


#Вычисление максимальной погрешности и отрисовка погрешности
error = np.empty((L+1,N+1), dtype="float")
for i in range(L+1):
    for j in range(N+1):
        error[i,j] = np.abs(u_sol[i,j]-u[i,j])

PlotError(x_step, t_step, u_sol)
print("Сonservative scheme Error: ", CalculateError(u, u_sol), "\n")



    