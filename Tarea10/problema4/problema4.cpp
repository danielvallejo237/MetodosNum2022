/*
Método de integración de montecarlo para integrales dobles
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

double findMax2D(FunctionParser fp,double a,double b,double c,double d,int type)
{
  /*
  Encontramos el máximo y el mínimo de una función f dentro de un intervalo a,b,
  el tipo nos determina si estamos considerando la parte positiva de la
  función o la parte negativa segun sea el caso
  */
  double max=0;
  vector<double> P=Linspace(a,b,1000);
  vector<double> Q=Linspace(c,d,1000);
  double e[2];
  if (type==1)
  {
    for(auto c: P)
    {
      for(auto q:Q)
      {
        e[0]=c;
        e[1]=q;
        max=fp.Eval(e)>max?fp.Eval(e):max;
      }
    }
  }
  else
  {
    for(auto c: P)
    {
      for(auto q:Q)
      {
        e[0]=c;
        e[1]=q;
        max=fp.Eval(e)<max?fp.Eval(e):max;
      }
    }
  }

  return max;
}

int sign(double k)
{
  return k/fabs(k);
}

bool IsPositive(double k)
{
  return sign(k)>0;
}

void writePairstoTxt(string file,vector<pair<double,pair<double,double>>> P)
{
  ofstream os{file};
  for(auto p:P)
  {
    os<<p.second.first<<","<<p.second.second<<","<<p.first<<endl;
  }
}

vector<pair<double,pair<double,double>>> GenRandom2d(double startX,double endX,double startY,double endY,double startZ, double endZ,int N)
{
  vector<pair<double,pair<double,double>>> Pts(N);
  #pragma omp parallel
  {
    srand(int(time(NULL)) ^ omp_get_thread_num());
    #pragma omp for
    for(int i=0;i<N;i++)
    {
        double r1,r2,r3;
        r1=random(startX,endX);
        r2=random(startY,endY);
        r3=random(startZ,endZ);
        Pts[i]=make_pair(r3,make_pair(r1,r2));
    }
  }
  //writePairstoTxt("pairs.txt",Pts);
  return Pts;
}



double integralOverVolume(FunctionParser fp, double a, double b,double c,double d)
{
  double max=findMax2D(fp,a,b,c,d,1);
  double min=findMax2D(fp,a,b,c,d,2);
  vector<pair<double,pair<double,double>>> P=GenRandom2d(a,b,c,d,max,min,1000000);
  double eval[2];
  int counter=0;
  for(auto p:P)
  {
    eval[0]=p.second.first;
    eval[1]=p.second.second;
    if(p.first>=0 &&fp.Eval(eval)>=p.first)counter+=1;
    if(p.first<0 &&fp.Eval(eval)<=p.first)counter-=1;
  }
  double volume=(b-a)*(d-c)*(max-min);
  double prop=(double)counter/(double)P.size();
  return prop*volume;
}

int main(int argc, char *argv[])
{
  omp_set_num_threads(MAX_NUM_THREADS); //Se tiene que hacer en paralelo para mayor eficiencia
  string expresion; //Expresion que contiene la función que deberá ser evaluada
  expresion=argv[1];
  double a,b;
  double completeIntegral=0;
  double ay,by;
  ay=atof(argv[4]);
  by=atof(argv[5]);
  a=atof(argv[2]);
  b=atof(argv[3]);
  FunctionParser fp;
  fp.AddConstant("pi",3.1415926535897932);
  fp.AddConstant("e",2.718281828459); //Valores más comunes encontrados en matemáticas
  fp.Parse(expresion,"x,y");
  double k=integralOverVolume(fp,a,b,ay,by);
  cout<<"Integral de "+string(argv[1])+" : "<<k<<endl;
  return 0;
}
