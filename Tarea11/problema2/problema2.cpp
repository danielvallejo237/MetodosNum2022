
/*Probar extrapolación de richardson para acelerar la integración
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/
#include<bits/stdc++.h>
#include"../fparser/fparser.hh"
using namespace std;

pair<double,vector<double>> Linspace(double start, double end, int N)
{
  double stepsize=(end-start)/(double)(N-1);
  vector<double> Points(N);
  for(int i=0;i<N;i++) Points[i]=start+stepsize*i;
  return make_pair(stepsize,Points); //Regresamos el tamaño de paso y los puntos de evaluación
}

double NewtonCotesD2(FunctionParser fp, double start, double end, int N)
{
  pair<double,vector<double>> pts=Linspace(start,end,N); //Numero defaulr de puntos
  double sumaextremos=fp.Eval(&pts.second[0])+fp.Eval(&pts.second[N-1]);
  double sumaOdd, sumaEven;
  sumaOdd=0;
  sumaEven=0;
  double stS=(pts.second[2]-pts.second[0]);
  for(int i=1;i<N-1;i++)
  {
    if (i%2==1) sumaOdd+=fp.Eval(&pts.second[i]);
    else sumaEven+=fp.Eval(&pts.second[i]);
  }
  return (stS/6.0)*sumaextremos+(4.0*stS/6.0)*sumaOdd+(2.0*stS/6.0)*sumaEven;
}

double NewtonCotesD2wRE(FunctionParser fp, double start, double end)
{
  double ih=NewtonCotesD2(fp,start,end,10000);
  double i2h=NewtonCotesD2(fp,start,end,5000);
  return ih+(1.0/15.0)*(ih-i2h);
}

int main(int argc, char *argv[])
{
  string expresion; //Expresion que contiene la función que deberá ser evaluada
  expresion=argv[1];
  double a,b;
  FunctionParser fp;
  fp.AddConstant("pi",3.1415926535897932);
  fp.AddConstant("e",2.718281828459); //Valores más comunes encontrados en matemáticas
  fp.Parse(expresion,"x");
  a=atof(argv[2]);
  b=atof(argv[3]);
  cout<<"Valor de la integral NC3: "<<NewtonCotesD2(fp,a,b,10000)<<endl;
  cout<<"Valor de la integral NC3 con extrapolacion de Richardson: "<<NewtonCotesD2wRE(fp,a,b)<<endl;
  return 0;
}
