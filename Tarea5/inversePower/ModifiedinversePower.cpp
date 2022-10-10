/*Algoritmo de la potencia Inversa implementado para matrices simétricas */

//@Author Daniel Vallejo Aldana, (danielvallejo237) on Github

#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>
#include<omp.h>

#define MAX_NUM_THREADS 4
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
  #pragma omp parallel for reduction (+:suma)
  for(int i=0;i<v.size();i++) suma+=v[i]*v[i];
  return suma;
}

Matrix VectorMult(vector<double> v1, vector<double> v2)
{
  Matrix M((int)v1.size(),(int)v2.size());
  #pragma omp parallel for
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
double getmin(vector<double> vec)
{
  double norma=fabs(vec[0]);
  for(int i=0;i<vec.size();i++) if(fabs(vec[i])<norma) norma=vec[i];
  return norma;
}
double multiplyVec(vector<double> a,vector<double> b)
{
  double c=0;
  #pragma omp parallel for reduction (+:c)
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
  #pragma omp parallel for
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
  #pragma omp parallel for
  for (int j=0;j<A.n;j++) Col[j]=A.get(i,j);
  return Col;
}

Matrix Forward_Substitution(Matrix A, Matrix b)
{
    vector<double> solutions(A.n);
    solutions[0]=b.matrix[0]/A.matrix[0]; //Caso base y asumiendo que jamás encontraremos una división por cero
    double value=0;
    for(int i=1;i<A.n;i++)
    {
        //Iteramos sobre todos los renglones que tenemos de la matriz
        for(int j=0;j<i;j++) value+=solutions[j]*A.matrix[i*A.m+j];
        solutions[i]=(b.matrix[i]-value)/A.matrix[i*(A.m+1)];
        value=0; //Regresamos a la inicialización del valor de valor
    }
    Matrix mat(solutions,A.n,1);
    return mat; //Regresamos la matriz de solución con sustitución hacia adelante
}
Matrix Backward_Substitution(Matrix A, Matrix b)
{
    vector<double> solutions(A.n);
    solutions[A.n-1]=b.matrix[A.n-1]/A.matrix[(A.n-1)*(A.m+1)]; //Caso base y asumiendo que jamás encontraremos una división por cero
    double value=0;
    for(int i=A.n-2;i>=0;i--)
    {
        //Iteramos sobre todos los renglones que tenemos de la matriz
        for(int j=i+1; j<A.m;j++) value+=solutions[j]*A.matrix[i*A.m+j];
        solutions[i]=(b.matrix[i]-value)/A.matrix[i*(A.m+1)];
        value=0; //Regresamos a la inicialización del valor de valor
    }
    Matrix mat(solutions,A.n,1);
    return mat; //Regresamos la matriz de solución con sustitución hacia adelante
}

pair<pair<Matrix,Matrix>,bool> LUDecomposition(Matrix &A)
{
    /*Recibimos una matriz y la descomponemos en dos matrices que satisfacen que una es triangular superior
     y la otra es triangular inferior*/
    vector<double> l1(A.n*A.m,0),u1(A.n*A.m,0);
    Matrix L(l1,A.n,A.m), U(u1,A.n,A.m);
    //Necesitamos definir que todos los elementos de la diagonal de L sean 1's, en este caso tenemos lo siguiente
    for(int i=0;i<A.n;i++) L.matrix[i*A.m+i]=1;
    U.matrix[0]=A.matrix[0]; //Inicialización de la descomposicion LU
    if(!U.matrix[0]) return make_pair(make_pair(L,U),false); //No es factorizable
    //Paso 2 para el primer renglón de U y la primera columna de L hacemos la siguiente inicialización
    for(int i=1;i<U.m;i++)U.matrix[i]=A.matrix[i]/L.matrix[0],L.matrix[i*A.m]=A.matrix[i*A.m]/U.matrix[0];
    for(int i=1;i<A.n;i++)
    {
        double value=0;
        for(int k=0;k<i;k++) value+=L.matrix[i*L.m+k]*U.matrix[k*U.m+i];
        U.matrix[i*U.m+i]=A.matrix[i*U.m+i]-value; //Para las diagonales de la matriz U
        if(!U.matrix[i*U.m+i]) return make_pair(make_pair(L,U),false);
        if (i<A.n-1)
        {
            for(int j=i+1;j<A.m;j++)
            {
                double value1,value2;
                value1=0,value2=0; //Inicialización con el neutro aditivo de los numeros
                for(int k=0;k<i;k++)
                {
                    value1+=L.matrix[i*L.m+k]*U.matrix[k*U.m+j];
                    value2+=L.matrix[j*L.m+k]*U.matrix[k*U.m+i];
                }
                U.matrix[i*U.m+j]=(1/L.matrix[i*L.m+i])*(A.matrix[i*A.m+j]-value1);
                L.matrix[j*L.m+i]=(1/U.matrix[i*U.m+i])*(A.matrix[j*A.m+i]-value2);
            }
        }
    }
    return make_pair(make_pair(L,U),true); //La factorización se hizo de forma correcta
}

vector<double> SolveLU(Matrix L, Matrix U, vector<double> x)
{
  Matrix tmp(x,x.size(),1);
  Matrix s1;
  Matrix s2;
  s1=Forward_Substitution(L,tmp);
  s2=Backward_Substitution(U,s1);
  vector<double> sal;
  sal=s2.matrix;
  return sal;
}
//pair<vector<double>,double>
pair<vector<double>,double> GetFirstInversePower(vector<double> v0,Matrix A,Matrix L, Matrix U,int maxiter=10000,double TOL=1e-4)
{
  vector<double> v1;
  double lambda=10000000;
  double old_lambda=100000000;
  double err=1000000;
  vector<double> aux(v0.size(),0);
  normalizeVector(v0);
  while(err>TOL)
  {
    v1=SolveLU(L,U,v0);
    lambda=computeLambda(v0,v1,A);
    err=fabs(lambda-old_lambda);
    old_lambda=lambda;
    normalizeVector(v1);
    v0=v1;
  }
  return make_pair(v1,lambda);
}

pair<vector<double>,double> IthInversePower(vector<double> v0,Matrix A,Matrix L, Matrix U,vector<vector<double>> M,int ind,int maxiter=10000,double TOL=1e-4)
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
    v1=SolveLU(L,U,restaVecs(v0,aux));
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

void ComputeEigs(Matrix A, int Numvals)
{
  vector<vector<double>> P(Numvals,vector<double> (A.m));
  vector<double> eigs(Numvals);
  pair<pair<Matrix,Matrix>,bool> LU=LUDecomposition(A);
  pair<vector<double>,double> out=GetFirstInversePower(get_column(A,0),A,LU.first.first,LU.first.second);
  P[0]=out.first;
  vector<double> vec;
  eigs[0]=out.second;
  cout<<"Eigenvalue 1: "<<eigs[0]<<endl;
  cout<<"Error ||Ax-lx||: "<<norm2(restaVecs(A*out.first,out.first*out.second))<<endl;
  for(int i=1;i<Numvals;i++)
  {
    out=IthInversePower(get_column(A,i),A,LU.first.first,LU.first.second,P,i);
    P[i]=out.first;
    eigs[i]=out.second;
    cout<<"Eigenvalue "<<(i+1)<<": "<<out.second<<endl;
    cout<<"Error ||Ax-lx||: "<<norm2(restaVecs(A*out.first,out.first*out.second))<<endl;
  }
}

int main(int argv, char* argc[])
{
  omp_set_num_threads(MAX_NUM_THREADS);
  Matrix A(argc[1]);
  ComputeEigs(A,atoi(argc[2]));
  return 0;
}
