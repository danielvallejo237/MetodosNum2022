/*
Método de integración de montecarlo para integrales simples
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/
#include<bits/stdc++.h>
#include "./fparser/fparser.hh"
#include <omp.h>
#define MAX_NUM_THREADS 3
using namespace std;


vector<double> Linspace(double a, double b, int N)
{
  double stepsize=(b-a)/N;
  vector<double> points(N+1);
  for(int i=0;i<=N;i++) points[i]=a+i*stepsize;
  return points;
}

double random(double range_from, double range_to)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> unif(range_from,range_to);
    return unif(generator);
}

double findMax1D(FunctionParser fp,double a, double b,int type)
{
  /*
  Encontramos el máximo y el mínimo de una función f dentro de un intervalo a,b,
  el tipo nos determina si estamos considerando la parte positiva de la
  función o la parte negativa segun sea el caso
  */
  double max=0;
  vector<double> P=Linspace(a,b,1000);
  double e[1];
  if (type==1)
  {
    for(auto c: P)
    {
      e[0]=c;
      max=fp.Eval(e)>max?fp.Eval(e):max;
    }
  }
  else
  {
    for(auto c: P)
    {
      e[0]=c;
      max=fp.Eval(e)<max?fp.Eval(e):max;
    }
  }

  return max;
}

double biseccion(FunctionParser fp, double x1, double x2,vector<double> &Info,int max_iters=100,double tol=1e-8)
{
  //Definimos el número máximo de iteraciones que se requieren en el algoritmo
  //De tolerancia definimos la raiz de epsilon para la tolerancia de nuestro algoritmo
  double a1[1],a2[1],aux[1];
  a1[0]=x1;
  a2[0]=x2;
  double xm=(x1+x2)/2;
  aux[0]=xm;
  int auxiter=max_iters;
  if(fabs(fp.Eval(a1))<tol)
  {
    cout<<"Salida en la primera evaluacion"<<endl;
    return a1[0]; //Casos triviales que no requiren posteriore evaluaciones
  }
  if(fabs(fp.Eval(a2))<tol) {cout<<"Salida en la segunda evaluacion"<<fp.Eval(a2)<<endl;return a2[0];}
  if (fp.Eval(a1)*fp.Eval(a2)>0)
  {
    cout<<"Los puntos deben de ser de diferente signo..."<<endl;
    exit(1);
  }
  while(fabs(fp.Eval(aux))>tol && max_iters-- && fabs(0.5*(a2[0]-a1[0]))>tol)
  {
    xm=0.5*(a1[0]+a2[0]),aux[0]=xm;
    Info.push_back(xm);
    if(fp.Eval(a1)*fp.Eval(aux)>0) a1[0]=xm;
    else a2[0]=xm;
  }
  return aux[0];
}

int sign(double k)
{
  return k/fabs(k);
}

bool IsPositive(double k)
{
  return sign(k)>0;
}
void writePairstoTxt(string file,vector<pair<double,double>> P)
{
  ofstream os{file};
  for(auto p:P)
  {
    os<<p.first<<","<<p.second<<endl;
  }
}
vector<pair<double,double>> GenRandom1d(double startX,double endX,double startY,double endY,int N)
{
  vector<pair<double,double>> Pts(N);
  #pragma omp parallel
  {
    srand(int(time(NULL)) ^ omp_get_thread_num());
    #pragma omp for
    for(int i=0;i<N;i++)
    {
      if (IsPositive(endY))
      {
        double r1,r2;
        r1=random(startX,endX);
        r2=random(startY,endY);
        Pts[i]=make_pair(r1,r2);
      }
      else
      {
        double r1,r2;
        r1=random(startX,endX);
        r2=random(startY,endY);
        Pts[i]=make_pair(r1,r2);
      }
    }
  }
  writePairstoTxt("parejas.txt",Pts);
  return Pts;
}


double Proportion(FunctionParser fp,vector<pair<double,double>> Pairs,int type)
{
  int counter=0;
  #pragma omp parallel for reduction(+:counter)
  for(auto p:Pairs)
  {
    if(type==1)
    {
      counter+=p.second<=fp.Eval(&p.first);
    }
    if(type==2)
    {
      counter+=p.second>=fp.Eval(&p.first);
    }
  }
  return (double)counter/(double)Pairs.size();
}


double integralOverArea(FunctionParser fp, double a, double b)
{
    double mx=findMax1D(fp,a,b,1);
    double mn=findMax1D(fp,a,b,2);
    double area=(b-a)*(mx-mn); //El area del cuadrado
    vector<pair<double,double>> Pt=GenRandom1d(a,b,mn,mx,1000000);
    double eval[2];
   int counter=0;
   double integral,prop;
   for(auto p:Pt)
  {
    eval[0]=p.first;
    if(p.second>=0 &&fp.Eval(eval)>=p.second)counter+=1;
    if(p.second<0 &&fp.Eval(eval)<=p.second)counter-=1;
  }
	prop=(double)counter/(double)Pt.size();
    integral=prop*area;
  return integral;
}

int main(int argc, char *argv[])
{
  omp_set_num_threads(MAX_NUM_THREADS); //Se tiene que hacer en paralelo para mayor eficiencia
  string expresion; //Expresion que contiene la función que deberá ser evaluada
  expresion=argv[1];
  double a,b;
  FunctionParser fp;
  fp.AddConstant("pi",3.1415926535897932);
  fp.AddConstant("e",2.718281828459); //Valores más comunes encontrados en matemáticas
  fp.Parse(expresion,"x");
  a=atof(argv[2]);
  b=atof(argv[3]);
  double I=integralOverArea(fp,a,b);
  cout<<"Integral de "+string(argv[1])+" de "<<a<<" hasta "<<b<<": "<<I<<endl;
  return 0;
}
