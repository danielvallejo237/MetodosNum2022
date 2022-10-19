// Cálculo de la matriz inversa para una matrix dada
/*@Author: Daniel Vallejo Aldana (danielvallejo237 on Github) */

#include <omp.h>
#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>

#define MAX_NUM_THREADS 3
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


pair<pair<Matrix,Matrix>,bool> LUDecompositionBandedMatrix(Matrix &A, int bandwidth)
{
    /*Recibimos una matriz y la descomponemos en dos matrices que satisfacen que una es triangular superior
     y la otra es triangular inferior*/
    vector<double> l1(A.n*A.m,0),u1(A.n*A.m,0);
    Matrix L(l1,A.n,A.m), U(u1,A.n,A.m);
    int resta=bandwidth/2;
    //Necesitamos definir que todos los elementos de la diagonal de L sean 1's, en este caso tenemos lo siguiente
    for(int i=0;i<A.n;i++) L.matrix[i*A.m+i]=1;
    U.matrix[0]=A.matrix[0]; //Inicialización de la descomposicion LU
    if(!U.matrix[0]) return make_pair(make_pair(L,U),false); //No es factorizable
    //Paso 2 para el primer renglón de U y la primera columna de L hacemos la siguiente inicialización
    for(int i=1;i<=resta;i++)U.matrix[i]=A.matrix[i]/L.matrix[0],L.matrix[i*A.m]=A.matrix[i*A.m]/U.matrix[0];
    for(int i=1;i<A.n;i++)
    {
        double value=0;
        for(int k=max(0,i-resta);k<i;k++) value+=L.matrix[i*L.m+k]*U.matrix[k*U.m+i];
        U.matrix[i*U.m+i]=A.matrix[i*U.m+i]-value; //Para las diagonales de la matriz U
        if(!U.matrix[i*U.m+i]) return make_pair(make_pair(L,U),false);
        if (i<A.n-1)
        {
            for(int j=i+1;j<min(A.m,i+resta+1);j++)
            {
                double value1,value2;
                value1=0,value2=0; //Inicialización con el neutro aditivo de los numeros
                for(int k=max(0,i-resta);k<i;k++)
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

Matrix operator * (Matrix A, Matrix B)
{
  Matrix C(A.n,B.m);
  if(A.m!=B.n) return A; //No es posible hacer la multiplicación de las matrices
  for(int i=0;i<C.n;i++)
  {
    for(int j=0;j<C.m;j++)
    {
      double suma=0;
      #pragma omp parallel for reduction(+:suma)
      for(int k=0;k<A.m;k++)
      {
        suma+=A.get(i,k)*B.get(k,j);
      }
      C.put(i,j,suma);
    }
  }
  return C;
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

Matrix ExtractDiagonal(Matrix &U, int bandwidth)
{
    int resta=bandwidth/2;
    vector<double> Diagonal(U.n);
    for(int i=0;i<U.n;i++)
    {
        Diagonal[i]=U.matrix[i*U.m+i];
        for (int j=i;j<min(U.m,i+resta+1);j++) U.matrix[i*U.m+j]=U.matrix[i*U.m+j]/Diagonal[i];
    }
    Matrix M(U.n,U.n);
    for(int i=0;i<U.n;i++) M.put(i,i,Diagonal[i]);
    return M;
}

void InvertDiagonal(Matrix &M)
{
  for(int i=0;i<M.n;i++) M.matrix[i*M.m+i]=1/M.matrix[i*M.m+i];
}

void ForwardDiagonalMatrixMult(Matrix &D, Matrix &M,int bandwidth)
{
  int resta=bandwidth/2;
  for(int i=0;i<M.n;i++)
  {
    for(int j=0;j<=i;j++)
    {
      M.matrix[i*M.m+j]=D.matrix[i*M.m+i]*M.matrix[i*M.m+j];
    }
  }
}

void ChangeOffDiagonalSign(Matrix &BB,int type,int bandwidth)
{
  int resta=bandwidth/2;
  if (type==1)
  {
    //Lower triangular
    for(int i=1;i<BB.n;i++) for(int j=max(0,i-resta);j<i;j++) BB.matrix[i*BB.m+j]=-1*BB.matrix[i*BB.m+j];
  }
  else
  {
    for(int i=0;i<BB.n-1;i++) for(int j=(i+1);j<min(BB.m,i+resta+1);j++) {BB.matrix[i*BB.m+j]=-1*BB.matrix[i*BB.m+j];}
  }
}


Matrix MultipyBandedMat(Matrix A, Matrix B, int bandwidth)
{
  int resta=bandwidth/2;
  Matrix C(A.n,B.m);
  if(A.m!=B.n) return A; //No es posible hacer la multiplicación de las matrices
  for(int i=0;i<C.n;i++)
  {
    for(int j=max(0,i-resta);j<min(C.m,i+resta+1);j++)
    {
      double suma=0;
      #pragma omp parallel for reduction(+:suma)
      for(int k=0;k<A.m;k++)
      {
        suma+=A.get(i,k)*B.get(k,j);
      }
      C.put(i,j,suma);
    }
  }
  return C;
}

vector<double> CanonicalVector(int N, int ind)
{
  vector<double> e(N,0);
  e[ind]=1.0;
  return e;
}

void PutColumn(Matrix &M,vector<double> C, int ind)
{
    for(int i=0;i<M.n;i++)
    {
      M.matrix[ind*M.m+i]=C[i];
    }
}


vector<double> operator - (vector<double> v1, vector<double> v2)
{
  vector<double> v3(v1.size(),0);
  #pragma omp parallel for
  for(int i=0;i<v3.size();i++) v3[i]=v1[i]-v2[i];
  return v3;
}
double Norma2(vector<double> v1,vector<double> v2)
{
  double suma=0;
  #pragma omp parallel for reduction (+:suma)
  for(int i=0;i<v1.size();i++)
  {
    suma+=(v1[i]-v2[i])*(v1[i]-v2[i]);
  }
  return suma;
}
Matrix ComputeInverseBanded(Matrix &A, int bandwidth)
{
    pair<pair<Matrix,Matrix>,bool> P=LUDecompositionBandedMatrix(A,bandwidth);
    Matrix Inversa(A.n,A.m);
    #pragma omp parallel for
    for(int i=0;i<A.m;i++)
    {
    	vector<double> tmp,tmp2;
      cout<<"Solving system: "<<i<<endl;
      tmp2=CanonicalVector(A.n,i);
      tmp=SolveLU(P.first.first,P.first.second,tmp2);
      PutColumn(Inversa,tmp,i);
    }
    Inversa.print_to_text("PartialInverse.txt");
    return Inversa;
}
double maxOOD(Matrix A)
{
  double max=-1000000;
  for(int i=1;i<A.n-1;i++) for(int j=i+1;j<A.m;j++) if(fabs(A.get(i,j))>max) max=fabs(A.get(i,j));
  return max;
}


double MaxD(Matrix A)
{
  double max=-1000000;
  for(int i=0;i<A.n;i++) if (fabs(A.get(i,i))>max) max=fabs(A.get(i,i));
  return max;
}
int main(int argc, char *argv[])
{
  omp_set_num_threads(MAX_NUM_THREADS);
  Matrix A(argv[1]);
  Matrix Inversa;
  Inversa=ComputeInverseBanded(A,7);
  cout<<"Inversa calculada"<<endl;
  Matrix C=A*Inversa;
  cout<<"Pruebas AA^-1"<<endl;
  cout<<"Elemento maximo fuera de la diagonal AA^-1: "<<maxOOD(C)<<endl;
  cout<<"Elemento maximo en la diagonal AA^-1: "<<MaxD(C)<<endl;
  Matrix D=Inversa*A;
  cout<<"Pruebas A^-1A"<<endl;
  cout<<"Elemento maximo fuera de la diagonal A^-1A: "<<maxOOD(D)<<endl;
  cout<<"Elemento maximo en la diagonal A^-1A: "<<MaxD(D)<<endl;
  return 0;
}
