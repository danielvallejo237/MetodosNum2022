/*
Método de integración de montecarlo para integrales simples
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/
#include<bits/stdc++.h>

vector<double> Linspace(double a, double b, int N)
{
  double stepsize=(b-a)/N;
  vector<double> points(N+1);
  for(int i=0;i<=N;i++) points[i]=a+i*stepsize;
  return points;
}

pair<double,double> findMax1D(double(*f)(double),double a, double b,int type)
{
  /*
  Encontramos el máximo y el mínimo de una función f dentro de un intervalo a,b,
  el tipo nos determina si estamos considerando la parte positiva de la
  función o la parte negativa segun sea el caso
  */
  double max=0;
  vector<double> P=Linspace(a,b,1000);
  for(auto c: P) max=c>max?c:max;
  return max;
}
