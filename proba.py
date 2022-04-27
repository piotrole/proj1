
from math import radians
import numpy as np


import sys
sys.path.append('C:\\Users\Piotr\Desktop\INFA\proj1\projekt')

from projekt import *

geo = Transformacje(model = "grs80")

plik = 'wsp_inp.txt'
# odczyt z pliku: https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.genfromtxt.html
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
#np.savetxt("wsp_out.txt", tablica, delimiter=',', fmt = ['%10.2f', '%10.2f', '%10.3f'], header = 'konversja współrzednych geodezyjnych \\ anna radzikowska')


tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
rows,cols = np.shape(tablica)

xyz2flh = np.zeros((rows,cols))

flh2xyz = np.zeros((rows,cols))

neu=np.zeros((rows,cols))

u2000 = np.zeros((rows,2))

u1992 = np.zeros((rows,2))

xyz2elewacja = np.zeros((rows,2))



for i in range(rows):
    
    xyz2flh[i] = geo.xyz2flh(tablica[i,0],tablica[i,1],tablica[i,2])
    
    flh2xyz[i] = geo.flh2xyz(radians(xyz2flh[i,0]),radians(xyz2flh[i,1]),(xyz2flh[i,2]))

    u1992[i] = geo.u1992(radians(xyz2flh[i,0]), radians(xyz2flh[i,1]))
    
    u2000[i] = geo.u2000(radians(xyz2flh[i,0]), radians(xyz2flh[i,1]))
    
    xyz2elewacja[i] = geo.xyz2elewacja(tablica[i,0],tablica[i,1],tablica[i,2],radians(xyz2flh[i,0]), radians(xyz2flh[i,1]),xyz2flh[i,2])



