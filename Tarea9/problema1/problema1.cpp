/* Interpolacion mediante splines lineales
Los splines lineales es simplemente ajustar una recta a los datos con base en nuestros puntos
de interpolación
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include<omp.h>
#include<bits/stdc++.h>

using namespace std;

void PrintMat(const vector<vector<double>> &vec)
{
  for(int i=0;i<vec.size();i++)
  {
    for(int j=0;j<vec[0].size();j++) cout<<vec[i][j]<<"\t";
    cout<<endl;
  }
}


double random(double range_from, double range_to)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> unif(range_from,range_to);
    return unif(generator);
}


double g(double x)
{
  return 3.45*sin(x)+2.5*sqrt(exp(x))-3*x;
}

double g2(double x)
{
  return x+ (x*sin(x/2.0))/3;
}

vector<double> Linspace(double Init, double End, int N)
{
  double step=(End-Init)/(double)(N-1);
  vector<double> steps(N);
  steps[0]=Init;
  steps[N-1]=End;
  for(int i=1;i<N-1;i++) {steps[i]=Init+i*step;}
  return steps;
}


vector<double> GenerateRandomPoints(int N, double Init , double End)
{
  vector<double> X(N,0);
  X[0]=Init;
  X[N-1]=End;
  for(int i=1;i<X.size()-1;i++)
  {
    X[i]=random(Init,End);
  }
  sort(X.begin(),X.end());
  return X;
}

pair<vector<double>,vector<double>> GeneratePoints(double(*f)(double),int N,double Init, double End)
{
  vector<double> X;
  X=GenerateRandomPoints(N+1,Init,End);
  cout<<X.size()<<endl;
  vector<double> Y(X.size(),0);
  double aux;
  for(int i=0;i<X.size();i++) Y[i]=f(X[i]);
  return make_pair(X,Y);
}


void SendtoText(string file,const vector<double> &v1,const vector<double> &v2)
{
  ofstream os{file};
  for(int i=0;i<v1.size();i++) os<<v1[i]<<","<<v2[i]<<endl;
}

void Vec2Text(string file,const vector<double> &v1)
{
  ofstream os{file};
  for(int i=0;i<v1.size();i++) os<<v1[i]<<endl;
}
void sendMat2Text(string file, const vector<vector<double>> &M)
{
  ofstream os{file};
  os<<M.size()<<" "<<M[0].size()<<endl;
  for(int i=0;i<M.size();i++)
  {
    for(int j=0;j<M[0].size();j++) os<<M[i][j]<<" ";
    os<<endl;
  }
}

vector<vector<double>> ComputeAB(vector<double> X, vector<double> Y)
{
  vector<vector<double>> M(X.size(),vector<double> (2,0));
  for(int i=0;i<X.size()-1;i++)
  {
    M[i][1]=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
    M[i][0]=Y[i];
    //cout<<M[i][1]<<endl;
  }
  M[X.size()-1][0]=Y[X.size()-1];
  return M;
}

double InterpolateValue(double val, vector<vector<double>> M,vector<double> X)
{
  int ind=0;
  for(int i=0;i<X.size();i++)
  {
    if(X[i]<=val)continue;
    else
    {
      ind=i-1;
      break;
    }
  }
  return M[ind][0]+M[ind][1]*(val-X[ind]); //El valor del punto en esa recta
}

double ComputeError(double(*f)(double),int points,double Init, double End, vector<vector<double>> Cf, vector<double> X)
{
  double suma=0;
  double a1,a2;
  double val;
  for(int i=0;i<points;i++)
  {
    val=random(Init,End);
    a1=InterpolateValue(val,Cf,X);
    a2=f(val);
    suma+=fabs(a2-a1);
  }
  return suma/points;
}



pair<vector<double>,vector<double>> GenerateInterpol(vector<double> Xs, vector<vector<double>> Cf, int N,double Init, double End)
{
  vector<double> pts=GenerateRandomPoints(N,Init,End);
  vector<double> ys(pts.size(),0);
  for(int i=0;i<ys.size()-1;i++) ys[i]=InterpolateValue(pts[i],Cf,Xs);
  return make_pair(pts,ys);
}

int main(int argc, char *argv[])
{
  pair<vector<double>,vector<double>> P=GeneratePoints(g2,atoi(argv[1]),atof(argv[2]),atof(argv[3]));
  vector<vector<double>> Smat=ComputeAB(P.first,P.second);
  pair<vector<double>,vector<double>> I=GenerateInterpol(P.first,Smat,1000,atof(argv[2]),atof(argv[3]));
  SendtoText("quadpuntos.txt",I.first,I.second);
  SendtoText("interpolpoints.txt",P.first,P.second);
  double error=ComputeError(g2,30,atof(argv[2]),atof(argv[3]),Smat,P.first);
  cout<<"Error |f(x)-fhat(h)| con 30 puntos aleatorios: "<<error<<endl;
  return 0;
}
