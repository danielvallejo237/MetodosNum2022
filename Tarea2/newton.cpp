/* Programación del método de newton, en este caso se asume que todas las funciones
evaluadas tiene expresiones analíticas en sus derivadas
Código escrito por Daniel Vallejo Aldana (danielvallejo237 on Github)
*/

#include<bits/stdc++.h>
#include<cmath>
#include"./fparser/fparser.hh"

using namespace std;

double newtonmethod(FunctionParser fp,FunctionParser dfp,double x0,int max_iters=100,double tol=1e-8)
{
  //Definimos el número máximo de iteraciones que se requieren en el algoritmo
  //De tolerancia definimos la raiz de epsilon para la tolerancia de nuestro algoritmo
  //En este caso debemos de recibir la derivada de la función
  double a1[1]; //En el auxiliar vamos a guardar el punto anterior que necesitaremos
  a1[0]=x0;
  double aux=x0;
  int auxiter=max_iters;
  do {
    /* Hacemos las primeras evaluaciones de x1 */
    aux=a1[0];
    if(dfp.Eval(a1)!=0) a1[0]=a1[0]-(fp.Eval(a1)/dfp.Eval(a1));
    else
      {
        cout<<"Division por cero encontrada"<<endl;
        exit(1);
      }
    }
  while(fabs(aux-a1[0])>tol && max_iters--);
  cout<<"Error |xt+1-xt|: "<<fabs(aux-a1[0])<<endl;
  cout<<"Iteraciones realizadas: "<<(auxiter-max_iters)<<endl;
  if (fabs(aux-a1[0])<tol) return a1[0];
  else
  {
    cout<<"Salida por máximo de iteraciones"<<endl;
    return a1[0]; //No hay solución dentro de los intervalos deseados y sale por máximo de iteraciones
  }
}

int main(int argc, char*argv[])
{
  double x1;
  x1=atof(argv[3]);
  FunctionParser fp;
  FunctionParser dfp;
  fp.AddConstant("pi",3.1415926535897932);
  fp.AddConstant("e",2.718281828459);
  dfp.AddConstant("pi",3.1415926535897932);
  dfp.AddConstant("e",2.718281828459);
  fp.Parse(argv[1],"x");
  dfp.Parse(argv[2],"x");
  double sol=newtonmethod(fp,dfp,x1);
  cout<<"Raiz de la funcion: "<<sol<<endl;
  return 0;
}
