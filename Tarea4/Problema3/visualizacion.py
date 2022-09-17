import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--file',type=str,help="Archivo con las soluciones")
if __name__=='__main__':
    args=parser.parse_args()
    with open(args.file) as f:
        sols=f.readlines()
    sols[0]='0'
    sols=[float(s) for s in sols]
    sols.append(2)
    X_axis=np.linspace(0,1,len(sols))
    Y_axis=np.asarray([x**2+x for x in list(X_axis)])
    plt.plot(X_axis,Y_axis,linewidth=2.5,color="#021648",alpha=0.7,label="Valor real")
    plt.plot(X_axis,np.asarray(sols),linewidth=3.0,color="#37060C",alpha=0.95,label="Valor estimado")
    plt.legend()
    plt.title("Aproximación con "+str(len(sols))+" nodos")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid((50,50),color='black',alpha=0.5)
    plt.savefig("Aproximacion_con_"+str(len(sols))+"_nodos.jpg")
