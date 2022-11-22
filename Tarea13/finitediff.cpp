#include <iostream>
#include <cstdlib>
#include <fstream>
#include "fparser/fparser.hh"

using namespace std;


double CenteredDifferenceC1(FunctionParser fp, double x, double step)
{
    double left=x-step;
    double right=x+step;
    return (fp.Eval(&right)-fp.Eval(&left))/(2*step);
}

double CenteredDifferenceC2(FunctionParser fp, double x, double step)
{
  double left=x-step;
  double right=x+step;
  double der=(CenteredDifferenceC1(fp,right,step)-CenteredDifferenceC1(fp,left,step))/(2*step);
  return der;
}

double CenteredDifferenceC3(FunctionParser fp, double x, double step)
{
  //Calcular la tercera derivada de una función con diferencias finitas centradas
  double right=x+step;
  double left=x-step;
  double der=(CenteredDifferenceC2(fp,right,step)-CenteredDifferenceC2(fp,left,step))/(2*step);
  return der;
}

double RightDifferenceC1(FunctionParser fp, double x, double step)
{
  double right=x+step;
  return (fp.Eval(&right)-fp.Eval(&x))/step;
}

double RightDifferenceC2(FunctionParser fp, double x, double step)
{
  double right=x+step;
  return (RightDifferenceC1(fp,right,step)-RightDifferenceC1(fp,x,step))/step;
}

double RightDifferenceC3(FunctionParser fp, double x, double step)
{
  double right=x+step;
  return (RightDifferenceC2(fp,right,step)-RightDifferenceC2(fp,x,step))/step;
}

double LeftDifferenceC1(FunctionParser fp, double x, double step)
{
  double left=x-step;
  return (fp.Eval(&x)-fp.Eval(&left))/step;
}

double LeftDifferenceC2(FunctionParser fp, double x, double step)
{
  double left=x-step;
  return (LeftDifferenceC1(fp,x,step)-LeftDifferenceC1(fp,left,step))/step;
}

double LeftDifferenceC3(FunctionParser fp, double x, double step)
{
  double left=x-step;
  return (LeftDifferenceC2(fp,x,step)-LeftDifferenceC2(fp,left,step))/step;
}

vector<double> linspace(double a, double b, int N)
{
  double stepsize=(b-a)/(double)(N-1);
  vector<double> points(N);
  points[0]=a;
  points[N-1]=b;
  for(int i=1;i<N-1;i++) points[i]=a+i*stepsize;
  return points;
}

void writeAll(const vector< vector<double> > &M)
{
  ofstream os("output.txt");
  for(int i=0;i<M.size();i++)
  {
    for(int j=0;j<M[i].size();j++) os<<M[i][j]<<",";
    os<<"\n";
  }
}

double MSE(const vector<vector<double> >&mat, int c1, int c2)
{
    double mse=0;
    for(int i=0;i<mat.size();i++)
    {
      mse+=(mat[i][c1]-mat[i][c2])*(mat[i][c1]-mat[i][c2]);
    }
    return mse/(double)mat.size();
}


int main(int argc, char *argv[])
{
  string expresion=argv[1];
  string derivative=argv[2];
  string secondder=argv[3];
  string thirder=argv[4];
  double a,b;
  a=atof(argv[5]); //Limits of derivation
  b=atof(argv[6]);
  int discretizacion=atoi(argv[8]);
  double step=atof(argv[7]);
  FunctionParser fp,dp,sdp,tdp;
  fp.Parse(expresion,"x");
  dp.Parse(derivative,"x");
  sdp.Parse(secondder,"x");
  tdp.Parse(thirder,"x");
  vector<double> toeval=linspace(a,b,discretizacion); //Evaluación en 1000 puntos para plotear la gráfica
  vector<vector<double> > Mat(toeval.size(),vector<double>(12,0));
  for(int i=0;i<toeval.size();i++)
  {
    Mat[i][0]=dp.Eval(&toeval[i]);
    Mat[i][1]=sdp.Eval(&toeval[i]);
    Mat[i][2]=tdp.Eval(&toeval[i]);
    Mat[i][3]=CenteredDifferenceC1(fp,toeval[i],step);
    Mat[i][4]=CenteredDifferenceC2(fp,toeval[i],step);
    Mat[i][5]=CenteredDifferenceC3(fp,toeval[i],step);
    Mat[i][6]=LeftDifferenceC1(fp,toeval[i],step);
    Mat[i][7]=LeftDifferenceC2(fp,toeval[i],step);
    Mat[i][8]=LeftDifferenceC3(fp,toeval[i],step);
    Mat[i][9]=RightDifferenceC1(fp,toeval[i],step);
    Mat[i][10]=RightDifferenceC2(fp,toeval[i],step);
    Mat[i][11]=RightDifferenceC3(fp,toeval[i],step);
  }

  cout<<"\n\n####### SUMMARY (DIFERENCIAS FINITAS CENTRADAS) #############"<<endl;
  cout<<"Error de diferencias centradas C1 con tamaño de paso "<<step<<" : "<<MSE(Mat,0,3)<<endl;
  cout<<"Error de diferencias centradas C2 con tamaño de paso "<<step<<" : "<<MSE(Mat,1,4)<<endl;
  cout<<"Error de diferencias centradas C3 con tamaño de paso "<<step<<" : "<<MSE(Mat,2,5)<<endl;
  cout<<"############ END SUMMARY ################################"<<endl;
  cout<<"\n\n";
  cout<<"####### SUMMARY (DIFERENCIAS FINITAS IZQUIERDAS) #############"<<endl;
  cout<<"Error de diferencias izquierdas C1 con tamaño de paso "<<step<<" : "<<MSE(Mat,0,6)<<endl;
  cout<<"Error de diferencias izquierdas C2 con tamaño de paso "<<step<<" : "<<MSE(Mat,1,7)<<endl;
  cout<<"Error de diferencias izquierdas C3 con tamaño de paso "<<step<<" : "<<MSE(Mat,2,8)<<endl;
  cout<<"############ END SUMMARY ################################"<<endl;
  cout<<"\n\n";
  cout<<"####### SUMMARY (DIFERENCIAS FINITAS DERECHAS) #############"<<endl;
  cout<<"Error de diferencias derechas C1 con tamaño de paso "<<step<<" : "<<MSE(Mat,0,9)<<endl;
  cout<<"Error de diferencias derechas C2 con tamaño de paso "<<step<<" : "<<MSE(Mat,1,10)<<endl;
  cout<<"Error de diferencias derechas C3 con tamaño de paso "<<step<<" : "<<MSE(Mat,2,11)<<endl;
  cout<<"############ END SUMMARY ################################"<<endl;
  cout<<"\n\n";
  writeAll(Mat);
  return 0;
}
