/*
Implementación del método de bisección,
recibe una función de entrada y regresa la solución de la raiz de la función si es que existe
Código escrito por Daniel Vallejo Aldana
Código escrito por Daniel Vallejo Aldana (danielvallejo237 on Github)
*/


#include <bits/stdc++.h>
#include <cmath>
#include "./fparser/fparser.hh"


using namespace std;

double biseccion(FunctionParser fp, double x1, double x2,int max_iters=100,double tol=1e-8)
{
  //Definimos el número máximo de iteraciones que se requieren en el algoritmo
  //De tolerancia definimos la raiz de epsilon para la tolerancia de nuestro algoritmo
  double a1[1],a2[1],aux[1];
  a1[0]=x1;
  a2[0]=x2;
  double xm;
  int auxiter=max_iters;
  if (fp.Eval(a1)*fp.Eval(a2)>0)
  {
    cout<<"Los puntos deben de ser de diferente signo..."<<endl;
    exit(1);
  }
  while(fabs(fp.Eval(a1)-fp.Eval(a2))>tol && max_iters--)
  {
    xm=0.5*(a1[0]+a2[0]),aux[0]=xm;
    if(fp.Eval(a1)*fp.Eval(aux)>0) a1[0]=xm;
    else a2[0]=xm;
  }
  cout<<"Error |f(x1)-f(x2)|: "<<fabs(fp.Eval(a1)-fp.Eval(a2))<<endl;
  cout<<"Iteraciones realizadas: "<<(auxiter-max_iters)<<endl;
  if (fabs(fp.Eval(a1)-fp.Eval(a2))<tol) return 0.5*(a1[0]+a2[0]);
  else
  {
    cout<<"Salida por máximo de iteraciones"<<endl;
    return 0.5*(a1[0]+a2[0]); //No hay solución dentro de los intervalos deseados y sale por máximo de iteraciones
  }
}

int main(int argc, char*argv[])
{
  string function=argv[1];
  double x1,x2;
  x1=atof(argv[2]);
  x2=atof(argv[3]);
  FunctionParser fp;
  fp.AddConstant("pi",3.1415926535897932);
  fp.AddConstant("e",2.718281828459);
  fp.Parse(function,"x");
  double sol=biseccion(fp,x1,x2);
  cout<<"Raiz de la funcion: "<<sol<<endl;
  return 0;
}
