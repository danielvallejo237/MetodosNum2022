/*Algoritmo de la potencia implementado para matrices simétricas */

//@Author Daniel Vallejo Aldana, (danielvallejo237) on Github

#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>
#include<omp.h>

#define MAX_NUM_THREADS 5

using namespace std;

//Incluimos las mismas librerías que en el problema anterior

class Matrix
{
    public:
    int n,m; //Dimensiones de la matriz
    vector<double> matrix;
    int type;
    Matrix(vector<double> &mat,int n, int m)
    {
        this->matrix=mat;
        this->n=n;
        this->m=m;
    }
    Matrix(int n, int m)
    {
      this->n=n;
      this->m=m;
      this->matrix.resize(n*m,0); //Inicializamos la matriz con puros 0's
    }
    Matrix(string filename)
    {
        //Podemos inicializar ahora una matriz con el nombre del archivo
        ifstream ifs(filename);
        string line;
        getline(ifs, line);
        istringstream iss(line);
        int n,m;
        iss>>n>>m;
        this->n=n;
        this->m=m;
        this->matrix.resize(this->n*this->m);
        int counter=0;
        double a;
        while (ifs >> a) this->matrix[counter]=a,counter++;
    }
    Matrix(){}
    ~Matrix(){}
    void set_dimensions(int n, int m)
    {
        this->n=n;
        this->m=m;
    }
    void set_vector(vector<double> &vec)
    {
        this->matrix=vec;
    }
    void print_content()
    {
        for(int i=0;i<this->n;i++)
        {
            for(int j=0;j<this->m;j++)
            {
                cout<<get(i,j)<<"\t";
            }
            cout<<endl;
        }
    }
    void print_to_text(string filename)
    {
        ofstream os{filename};
        os<<this->n<<" "<<this->m<<endl;
        for(int i=0;i<this->n;i++)
        {
            for(int j=0;j<this->m;j++) os<<this->matrix[i*this->m+j]<<" ";
            os<<endl;
        }
    }
    void swap_row(int a, int b)
    {
        //Funcion para cambiar los renglones de la matriz en caso de pivoteo
        for(int i=0;i<this->m;i++) swap(matrix[this->m*a+i],this->matrix[this->m*b+i]);
    }
    Matrix Transpose()
    {
        vector<double> T(this->n*this->m);
        for(int i=0;i<this->n;i++)for(int j=0;j<this->m;j++) T[j*this->n+i]=this->matrix[i*this->m+j];
        Matrix M(T,this->m,this->n);
        return M;
    }
    double get(int row,int col)
    {
        return this->matrix[row*this->m+col];
    }
    void put(int row, int col, double value)
    {
        this->matrix[row*this->m+col]=value;
    }
    pair<int,int> shape()
    {
        return make_pair(this->n,this->m);
    }
    //Sobrecarga de operadores para sumar matrices en lugar de sumarlos de una sola vez
    Matrix operator + (Matrix const &obj)
    {
      Matrix Sol(obj.n,obj.m);
      if(this->n==obj.n && this->m==obj.m)
      {
          //Tienen que tener las mismas dimensiones
          #pragma omp parallel for
          for(int i=0;i<this->n*this->m;i++) Sol.matrix[i]=this->matrix[i]+obj.matrix[i];
      }
      return Sol;
    }
    Matrix operator * (const double num)
    {
      //Necesitamos definir la multiplicación por un escalar
      Matrix Sol(this->n,this->m);
      for(int i=0;i<this->n*this->m;i++) Sol.matrix[i]=this->matrix[i]*num;
      return Sol;
    }

    vector<double> operator * (vector<double> const &obj)
    {
      vector<double> Sol(this->n,0);
      #pragma omp parallel for
      for(int j=0;j<this->n;j++)
      {
          for(int i=0;i<this->m;i++)
          {
            Sol[j]+=get(j,i)*obj[i];
          }
        }
      return Sol;
    }
    Matrix operator - (Matrix const &obj)
    {
        Matrix Sol(obj.n,obj.m);
        if(this->n==obj.n && this->m==obj.m)
        {
            //Tienen que tener las mismas dimensiones
            #pragma omp parallel for
            for(int i=0;i<this->n*this->m;i++) Sol.matrix[i]=this->matrix[i]-obj.matrix[i];
        }
        return Sol;
    }
    double norm()
    {
        double norma;
        if(this->n>1 && this->m>1)
        {
            norma=0;
            for(int i=0;i<this->n;i++)
            {
                double aux=0;
                for(int j=0;j<this->m;j++) aux+=fabs(get(i,j));
                if (aux>norma) norma=aux;
            }
        }
        else if (this->m==1 && this->n>1)
        {
            norma=0;
            for(int i=0;i<this->n;i++) norma+=get(i,0)*get(i,0);
            norma=sqrt(norma);
        }
        else if(this->m>1 && this->n==1)
        {
            norma=0;
            for(int i=0;i<this->m;i++) norma+=get(0,i)*get(0,i);
            norma=sqrt(norma);
        }
        else
        {
            norma=0;
        }
        return norma;
    }
};

double norm2(vector<double> v)
{
  //regresa la norma al cuadrado del
  double suma=0;
  #pragma omp parallel for reduction(+:suma)
  for(int i=0;i<v.size();i++) suma+=v[i]*v[i];
  return suma;
}

Matrix VectorMult(vector<double> v1, vector<double> v2)
{
  Matrix M((int)v1.size(),(int)v2.size());
  for(int i=0;i<v1.size();i++) for(int j=0;j<v2.size();j++) M.put(i,j,v1[i]*v2[j]);
  return M;
}

void normalizeVector(vector<double> &vec)
{
  double norma=0;
  #pragma omp parallel for reduction(+:norma)
  for(int i=0;i<vec.size();i++) norma+=vec[i]*vec[i];
  norma=sqrt(norma);
  #pragma omp parallel for
  for(int i=0;i<vec.size();i++) vec[i]/=norma;
}
double infiniteNorm(vector<double> vec)
{
  double norma=fabs(vec[0]);
  for(int i=0;i<vec.size();i++) if(fabs(vec[i])>norma) norma=fabs(vec[i]);
  return norma;
}
double getmax(vector<double> vec)
{
  double norma=fabs(vec[0]);
  for(int i=0;i<vec.size();i++) if(fabs(vec[i])>norma) norma=vec[i];
  return norma;
}
double multiplyVec(vector<double> a,vector<double> b)
{
  double c=0;
  #pragma omp parallel for reduction(+:c)
  for(int i=0;i<a.size();i++) c+=a[i]*b[i];
  return c;
}
double computeLambda(vector<double> v, vector<double> v2,Matrix A)
{
  double denom=multiplyVec(v,A*v);
  return denom;
}
vector<double> restaVecs(vector<double> v, vector<double> v2)
{
  vector<double> aux(v.size(),0);
  #pragma omp parallel for
  for(int i=0;i<v.size();i++) aux[i]=v[i]-v2[i];
  return aux;
}

vector<double> operator * (const vector<double>& v1, double v2)
{
    vector<double> v3(v1.size());
    #pragma omp parallel for
    for(int i=0;i<v1.size();i++) v3[i]=v1[i]*v2;
    return v3;
}

pair<vector<double>,double> LastEigenvalue(vector<double> v0,Matrix A,double TOL=1e-4,int maxiter=1000)
{
   //Esto calcula el valor propio y el vector propio más grande asociado
    vector<double> v1;
    double lambda=10000000;
    double old_lambda=100000000;
    double err=1000000;
    double mu;
    normalizeVector(v0);
    while(err>TOL)
    {
      v1=A*v0;
      lambda=computeLambda(v0,v1,A);
      err=fabs(lambda-old_lambda);
      old_lambda=lambda;
      normalizeVector(v1);
      v0=v1;
    }
    return make_pair(v1,lambda);
}

vector<double> Proyecta(vector<double> v1, vector<double> v2)
{
  vector<double> aux(v2.size());
  double norma=norm2(v2);
  double prod=multiplyVec(v1,v2);
  #pragma omp parallel for
  for(int i=0;i<v1.size();i++) aux[i]=(prod/norma)*v2[i];
  return aux;
}

vector<double> SumaVecs(vector<double> v1, vector<double> v2)
{
  vector<double> aux(v1.size(),0);
  for(int i=0;i<v1.size();i++) aux[i]=v1[i]+v2[i];
  return aux;
}

void printInfo(pair<vector<vector<double>>,vector<double>> Info)
{
  for (int i=0;i<Info.second.size();i++)
  {
    cout<<"Eigenvalue: "<<Info.second[i]<<endl;
    cout<<"Corresponding eigenvector: "<<endl;
    for(int j=0;j<Info.first[i].size();j++) cout<<Info.first[i][j]<<" ";
    cout<<endl;
  }
}

vector<double> get_column(Matrix A, int i)
{
  vector<double> Col(A.n,0);
  for (int j=0;j<A.n;j++) Col[j]=A.get(i,j);
  return Col;
}

pair<vector<double>,double> GeiIthEig(vector<double> v0,Matrix A,vector<vector<double>>M,int ind,double TOL=1e-4,int maxiter=1000)
{
  vector<double> v1;
  double lambda=10000000;
  double old_lambda=100000000;
  double err=1000000;
  vector<double> aux(v0.size(),0);
  normalizeVector(v0);
  while(err>TOL)
  {
    fill(aux.begin(),aux.end(),0);
    for(int i=0;i<ind;i++)
    {
      aux=SumaVecs(aux,Proyecta(v0,M[i]));
    }
    v1=A*restaVecs(v0,aux);
    lambda=computeLambda(v0,v1,A);
    err=fabs(lambda-old_lambda);
    old_lambda=lambda;
    normalizeVector(v1);
    v0=v1;
  }
  return make_pair(v1,lambda);
}

void printVec(vector<double> v)
{
  for(vector<double>::iterator it=v.begin();it!=v.end();++it) cout<<*it<<" ";
  cout<<endl;
}

vector<double> escalarDot(double l, vector<double> a)
{
  vector<double> aux(a.size(),0);
  for(int i=0;i<aux.size();++i) aux[i]=l*a[i];
  return aux;
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

void ComputeEigs(Matrix A, int Numvals)
{
  vector<vector<double>> P(Numvals,vector<double> (A.m));
  vector<double> eigs(Numvals);
  pair<vector<double>,double> out=LastEigenvalue(get_column(A,A.n-1),A);
  P[0]=out.first;
  eigs[0]=out.second;
  vector<double> vec;
  cout<<"Eigenvalue 1: "<<eigs[0]<<endl;
  cout<<"Error ||Ax-lx||: "<<norm2(restaVecs(A*out.first,out.first*out.second))<<endl;
  for(int i=1;i<Numvals;i++)
  {
    out=GeiIthEig(get_column(A,A.n-i-1),A,P,i);
    P[i]=out.first;
    eigs[i]=out.second;
    cout<<"Eigenvalue "<<(i+1)<<": "<<out.second<<endl;
    cout<<"Error ||Ax-lx||: "<<norm2(restaVecs(A*out.first,out.first*out.second))<<endl;
  }
  sendMat2Text("EvecsBiggest.txt",P);
}

int main(int argc, char* argv[])
{
  Matrix A(argv[1]);
  omp_set_num_threads(MAX_NUM_THREADS);
  ComputeEigs(A,atoi(argv[2]));
  return 0;
}
