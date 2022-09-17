/*
Implementación de una clase de matriz de banda de tal forma que pueda almacenar y calcular
de forma eficiente una matriz enorme sin necesidad de hacer las N^{2} operaciones
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>

//Implementación de una estructura de una matriz de banda

using namespace std;

class MatBanda
{
public:
    int n,m,bandwidth;
    vector<double> matrix;
    MatBanda(){}
    MatBanda(int n,int m,int bandwidth)
    {
      //Necesitamos el ancho de banda para poder almacenar la matriz de forma eficiente
      this->n=n;
      this->m=m;
      this->matrix.resize(bandwidth*n-(0.5*(bandwidth-1)*bandwidth),0); //Nuevo tamaño de la matriz es de n + n-1,
      //en este caso vamos a suponer que la matriz es simétrica
      this->bandwidth=bandwidth;
    }
    MatBanda(string filename, int bandwidth)
    {
      this->bandwidth=bandwidth;
      ifstream ifs(filename);
      string line;
      getline(ifs, line);
      istringstream iss(line);
      int n,m;
      iss>>n>>m;
      this->n=n;
      this->m=m;
      int counter=0;
      double a;
      vector<double> tmp(n*m,0);
      while (ifs >> a) tmp[counter]=a,counter++;
      //Una vez guardada toda la matriz procederemos a guardar solo los elementos que nos importan
      this->matrix.resize(bandwidth*n-(0.5*(bandwidth-1)*bandwidth),0);
      int counter2=0;
      for(int i=0;i<bandwidth;i++) for (int j=0;j<(n-i);j++) this->matrix[i*n+j]=tmp[j*m+(j+i)],counter2++;
    }
    ~MatBanda(){}
    void printContent()
    {
      int counter2=0;
      for(int i=0;i<bandwidth;i++)
      {
        for (int j=0;j<(n-i);j++)
        {
          cout<<this->matrix[counter2]<<" ";
          counter2++;
        }
        cout<<endl;
      }
    }
    void CholeskyFactorize()
    {
      MatBanda cholesky(this->n,this->m,this->bandwidth);
      cholesky.matrix[0]=sqrt(this->matrix[0]); //Primer elemento de la factorización de cholesky
      double suma=0;
      for(int i=1;i<(this->n-1);i++)
      {
        suma=0;
        for(int j=1;j<bandwidth;j++) suma+=this->matrix[this->n*i+j];
        cholesky.matrix[i]=sqrt(this->matrix[i]-suma);
      }
      //Las element of the diagonal
      cholesky.printContent();
    }
};

int main(int argc, char *argv[])
{
  MatBanda M(argv[1],2); //Ancho de banda de dos es decir una matriz tridiagonal
  //M.printContent();
  M.CholeskyFactorize();
  return 0;
}
