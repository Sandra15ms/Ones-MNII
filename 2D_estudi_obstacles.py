import numpy as np
import matplotlib.pyplot as plt


O1 = 150                # vores externes. Per no tenir-les en compte, posar valor de M i N
O2 = 150

M=150  #dimensió matriu # M columnes
N=150                   # N files
t=1000  
im_M = M
im_N = N
cx = 0

dt=0.0001             #discretitzacions                                        #dt_real = dt*t_0 = 0.0001*1s = 0.0001 s
dx = 0.0003                                                                    #dx_real = dx*x_0 = dx*t_0*v_p = 0.0003*1s*343m/s = 0.1029 m 

C = 0.00001       #constants
A_0 = 70         #amplitud constant
freq = 100       #freq ct


alpha=(dt/dx)**2


T1=np.zeros((O2+1,O1+1))         # matrius inicials
T1o=np.zeros((O2+1,O1+1))        # T1n1 és a t n+1
T1n1=np.zeros((O2+1,O1+1))  
T1e=np.zeros((O2+1,O1+1))  
show = np.zeros((O2+1, O1+1))     # parets
im_plot = np.zeros((im_N+1,im_M+1))
im_show = np.zeros((im_N+1,im_M+1))

tick=np.linspace(0,150,11)
real_pos_x=tick/3
real_pos_y=real_pos_x[::-1]

if M != O1:
    M = int(M + O1/4)
    N = int(N + O2/4)
    cx = int(O1/4)



# ---------------------------------- parets

def contorn():
    
# Parets
  for i in range(cx, M+1):           #condicions de contorn - paret superior j=0  
      T1n1[cx,i]=0
      show[cx,i]=70

  for i in range(cx, M+1):           #condicions de contorn - paret inferior j=N    
      T1n1[N,i]=0
      show[N,i]=70
  for j in range(cx, N+1):           #condicions de contorn - paret esquerra i=0  
      T1n1[j,cx]=0
      show[j,cx]=70

  for j in range(cx, N+1):           #condicions de contorn - paret dreta i=M    
      T1n1[j,M]=0
      show[j,M]=70  
      
      
#obstacles  
  for j in range(10,N-10):            
    T1n1[j,M-10]=0
    show[j,M-10]=70

  for j in range(50,N-40):            
    T1n1[j,M-50]=0
    show[j,M-50]=70     

  for j in range(70,N-40):            
    T1n1[j,M-80]=0
    show[j,M-80]=70   

  for j in range(50,71):            
    T1n1[j,M-65]=0
    show[j,M-65]=70  

    
    
  for j in range(100,N-20):            
    T1n1[j,M-30]=0
    show[j,M-30]=70     
    
  for j in range(100,N-20):            
    T1n1[j,30]=0
    show[j,30]=70     
        
    
  for i in range(30,M-30+1):
    T1n1[N-20,i]=0
    show[N-20,i]=70 

  for i in range(30,M-30+1):
    T1n1[N-20,i]=0
    show[N-20,i]=70       

  for i in range(M-80,M-65+1):
    T1n1[70,i]=0
    show[70,i]=70 

  for i in range(M-65,M-50+1):
    T1n1[50,i]=0
    show[50,i]=70 
    
  for i in range(M-80,M-50+1):
    T1n1[N-40,i]=0
    show[N-40,i]=70
    
# --------------------------------    

   

T1[50+cx,50+cx]=A_0*np.sin(2*np.pi*(freq*dt))           #fonts "converses" contínues


for j in range(0+cx,N):        #calculem el temps 1 imposant que t -1 = t0 -dt*C
  for i in range(0+cx,M):
    T1o[j,i]=T1[j,i]
    T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(T1[j,i])-(dt*C)

contorn() 

for i in range(0+cx,M):
  for j in range(0+cx,N+1):
    T1[j,i]=T1n1[j,i]

for n in range(2,t):         #comencem a calcular els temps n+1
 

  if n < int(1/(2*dt*freq)):
    for j in range(0,N+1):    
      T1[50+cx,50+cx]=A_0*np.sin(2*np.pi*(freq*n*dt)) 
  else:
    for j in range(0,N+1):
      T1[50+cx,50+cx]=0           #fonts "converses" contínues


  
  for i in range(0,O1):
    for j in range(0,O2+1):
      T1e[j,i]=T1o[j,i]

  for i in range(0,O1):
    for j in range(0,O2+1):
      T1o[j,i]=T1[j,i]
      if j in range (1,O2):
        T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(2*T1[j,i])-(T1e[j,i])      
        
  contorn()
  
  for i in range(0,O1):
    for j in range(0,O2+1):
          T1[j,i]=T1n1[j,i]

  t_past = n*dt             #temps real transcorregut
  
  for x in range(0,im_M+1):
      for y in range(0,im_N+1):
          im_plot[y,x]=np.abs(T1[y+cx,x+cx])
          im_show[y,x]=show[y+cx,x+cx]


  
      
  if n% 5==0:
      show = np.ma.masked_where(show == 0, show)
      fig = plt.subplots(1,1, figsize=(6,6))
           

      ax = plt.axes()
      ax.set_xlabel('x (m)')
      ax.set_ylabel('y (m)')
      ax.set_xticks(ticks=tick)
      ax.set_xticklabels(labels=real_pos_x)
      ax.set_yticks(ticks=tick)
      ax.set_yticklabels(labels=real_pos_y)
      
      plt.imshow(im_plot,cmap="gnuplot2", vmin = 0, vmax = 70)
      plt.imshow(im_show,cmap="gnuplot2", vmin = 0, vmax = 70, interpolation = 'none')
      plt.title("t = %f s" %t_past)
      plt.axes(ax)
      plt.colorbar()
#      plt.savefig('2D_obstacles_'+str(int(n)*0.00001)+'.png')      
      plt.show()