import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

M=151  #dimensió matriu # M columnes
N=101                   # N files
t=500  #temps

v_p=3.43

dt=0.0001 #discretitzacions
dx=0.0003

C = 0.001       #constants
A_0 = 0.8
freq = 8*v_p/(5*(M+1)*dx)

alpha=(dt/dx)**2


T1=np.zeros((N+1,M+1))         #matrius inicials
T1o=np.zeros((N+1,M+1))      # T1n1 és a t n+1
T1n1=np.zeros((N+1,M+1))  
T1e=np.zeros((N+1,M+1))

U=np.zeros((M+1,1))

plotter = np.zeros((M+1,1))  

for j in range(1,N):
  T1[j,75]=A_0*np.sin(2*np.pi*(freq*dt))

for j in range(0,N):        #calculem el temps 1 imposant que t -1 = t0 -dt*C
  for i in range(0,M):
    T1o[j,i]=T1[j,i]
    T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(T1[j,i])-(dt*C)

for i in range(0,M):
  for j in range(0,N+1):
    T1[j,i]=T1n1[j,i]

for n in range(2,t):         #comencem a calcular els temps n+1

#-------------------------------------------------------------------Analítica
    
  for i in range(0,M-1):
      if i in range(0,int(M/2)):
            U[i,0]=A_0*np.sin(2*np.pi*i*dx*np.pi*freq/v_p)*np.cos(2*np.pi*freq*dt*n)
      else:
            U[i,0]=A_0*np.sin(2*np.pi*np.abs(i-M/2)*dx*np.pi*freq/v_p)*np.cos(2*np.pi*freq*dt*n)
            
#------------------------------------------------------------------Numèrica
 
  for j in range(0,N+1):            #condicions de contorn - columnes i=0 ; i=M
    T1[j,0]=0
    T1[j,M]=0


  if n < int(2*t/(4*freq)):
    for j in range(0,N+1):    
      T1[j,75]=A_0*np.sin(2*np.pi*(freq*n*dt))
  else:
    for j in range(0,N+1):
      T1[j,75]=A_0*np.sin(2*np.pi*(freq*n*dt))      
      
  for i in range(0,M):
    for j in range(0,N+1):
      T1e[j,i]=T1o[j,i]

  for i in range(0,M):
    for j in range(0,N+1):
      T1o[j,i]=T1[j,i]
      if j==N:
        T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[0,i]-4*T1[j,i]+T1[j-1,i])+(2*T1[j,i])-(T1e[j,i])          # j+1 -> 0
      if j==0:
        T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[N,i])+(2*T1[j,i])-(T1e[j,i])          # j-1 -> N
      if j in range (1,N):
          T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(2*T1[j,i])-(T1e[j,i])     
  
  for i in range(0,M):
    for j in range(0,N+1):
      T1[j,i]=T1n1[j,i]
    plotter[i,0]=T1[50,i] #guardar u 

  if n%5 == 0:
      figure(figsize=(10,4))      
      ax=plt.gca()
      ax.set_ylim([0,1])    
      plt.plot(np.abs(plotter), color='red', scaley = False)
      plt.plot(np.abs(U), color='gray', scaley = False,linestyle = 'dashed')
      plt.savefig('error_'+str(int(n)*0.00001)+'.png') 
      plt.show()

  