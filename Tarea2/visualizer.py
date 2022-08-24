#Visualizador de error relativo a lo largo de las iteraciones

import numpy as np
import itertools
import os
import matplotlib.pyplot as plt
from datetime import datetime
import argparse

with open("puntos.txt","r") as f:
    datos=f.read()
    datos=[float(d) for d in datos.split()]

parser=argparse.ArgumentParser()
parser.add_argument('--root',type=float,help='Raiz real de la funcion')



if __name__=='__main__':
    args=parser.parse_args()
    raiz_real=args.root
    roots=np.asarray([abs((raiz_real-datos[i]))/max(1,abs(raiz_real)) for i in range(len(datos))])
    image_name="Imagen_"+str(datetime.now())+'.jpg'
    plt.title("Error relativo a lo largo de las iteraciones")
    plt.plot(roots,linewidth=2,color="#04163A")
    plt.xlabel("Iteraciones")
    plt.ylabel("Error relativo")
    plt.grid((50,50),color='black',alpha=0.5)
    plt.savefig("./"+image_name)