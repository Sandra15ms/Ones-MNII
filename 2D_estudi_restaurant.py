import numpy as np
import matplotlib.pyplot as plt
from numpy import random 

O1 = 300                # vores externes. Per no tenir-les en compte, posar valor de M i N
O2 = 300

M=150  #dimensió matriu # M columnes
N=150                   # N files
t=1000  
im_M = M
im_N = N
cx = 0

dt=0.0001             #discretitzacions                                        #dt_real = dt*t_0 = 0.0001*1s = 0.0001 s
dx = 0.0003                                                                    #dx_real = dx*x_0 = dx*t_0*v_p = 0.0003*1s*343m/s = 0.1029 m 

C = 0.00001       #constants
A_0 = random.normal(loc = 60, scale=3.87, size = (4,t))           #amplituds amb variabilitat
freqs = random.normal(loc = 170, scale=9, size = (4,t))           #frequencies amb variabiitat per les quatre fonts en el temps t

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
  for i in range(10+cx, M-10):           #condicions de contorn - paret superior j=20  
    if i not in range(70+cx,100+cx):
      T1n1[20+cx,i]=0
      show[20+cx,i]=70

  for i in range(10+cx, M-60):           #condicions de contorn - paret inferior j=N-20     
    if i not in range(50+cx,80+cx):
      T1n1[N-20,i]=0
      show[N-20,i]=70
      

      
  for j in range(20+cx,N+1-20):            #condicions de contorn - paret esquerra i=10
    T1n1[j,10+cx]=0
    show[j,10+cx]=70

  for j in range(20+cx,N+1-20):            #condicions de contorn - paret dreta i=M-10
    T1n1[j,M-10]=0
    show[j,M-10]=70
    
def obstacles():
    
  for i in range(0+cx,M):
    for j in range(0+cx,N+1):
        if (j in range(120+cx,N+1-20) and i in range(M-60,M-10)) or (j in range(90+cx,N+1-20) and i in range(M-25,M-10)):
            T1n1[j,i]=0
            show[j,i]=70
            
        if j==50+cx and i in range(10+cx,45+cx):
            T1n1[j,i]=0
            show[j,i]=70  
            
    
# --------------------------------    
    

T1[75+cx,40+cx]=A_0[0,0]*np.sin(2*np.pi*(freqs[0,0]*dt))           #fonts "converses" contínues
T1[75+cx,100+cx]=A_0[1,0]*np.sin(2*np.pi*(freqs[1,0]*dt))
T1[115+cx,40+cx]=A_0[2,0]*np.sin(2*np.pi*(freqs[2,0]*dt))           
T1[35+cx,100+cx]=A_0[3,0]*np.sin(2*np.pi*(freqs[3,0]*dt))



for j in range(0+cx,N):        #calculem el temps 1 imposant que t -1 = t0 -dt*C
  for i in range(0+cx,M):
    T1o[j,i]=T1[j,i]
    T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(T1[j,i])-(dt*C)

contorn()

for i in range(0+cx,M):
  for j in range(0+cx,N+1):
    T1[j,i]=T1n1[j,i]

for n in range(2,t):         #comencem a calcular els temps n+1
 

  T1[75+cx,40+cx]=A_0[0,n]*np.sin(2*np.pi*(freqs[0,n]*n*dt))           #fonts "converses" contínues
  T1[75+cx,100+cx]=A_0[1,n]*np.sin(2*np.pi*(freqs[1,n]*n*dt))
  T1[115+cx,40+cx]=A_0[2,n]*np.sin(2*np.pi*(freqs[2,n]*n*dt))
  T1[35+cx,100+cx]=A_0[3,n]*np.sin(2*np.pi*(freqs[3,n]*n*dt))

  
  for i in range(0,O1):
    for j in range(0,O2+1):
      T1e[j,i]=T1o[j,i]

  for i in range(0,O1):
    for j in range(0,O2+1):
      T1o[j,i]=T1[j,i]
      if j in range (1,O2):
        T1n1[j,i]=alpha*(T1[j,i+1]+T1[j,i-1]+T1[j+1,i]-4*T1[j,i]+T1[j-1,i])+(2*T1[j,i])-(T1e[j,i])      
        
  contorn()
  obstacles()      
    
  
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
#      plt.savefig('2D_restaurant_'+str(int(n)*0.00001)+'.png')      
      plt.show()



      

            
            
            
            
            
            
            
            
            
# https://courses.lumenlearning.com/suny-osuniversityphysics/chapter/16-2-mathematics-of-waves/
            #https://www.sltinfo.com/loudness/