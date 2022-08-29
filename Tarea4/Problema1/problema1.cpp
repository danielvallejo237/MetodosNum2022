
/* Programación del método de Cholesky para matrices A definidas positivas
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>
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

pair<Matrix,bool> Cholesky(Matrix A)
{
    //Implementación del método de cholesky
    /* Recibe una matriz A y regresa una matriz L que posteriormente será transpuesta con la función transpose implementada*/
    Matrix L(A.n,A.m);
    double value=0;
    if (A.get(0,0)<0) return make_pair(L,false);
    L.put(0,0,sqrt(A.get(0,0)));
    for(int j=1;j<L.n;j++) L.put(j,0,A.get(j,0)/L.get(0,0));
    for(int i=1;i<L.n-1;i++)
    {
      value=0;
      for(int k=0;k<i;k++) value+=L.get(i,k)*L.get(i,k);
      if((A.get(i,i)-value)<0) return make_pair(L,false); //No se pudo hacer la factorización de Cholesky porque no es definida positiva
      L.put(i,i,sqrt(A.get(i,i)-value));
      for(int j=i+1;j<L.m;j++)
      {
        value=0;
        for(int k=0;k<i;k++) value+=L.get(j,k)*L.get(i,k);
        L.put(j,i,(A.get(j,i)-value)/L.get(i,i));
      }
    }
    value=0;
    for(int k=0;k<L.n;k++) value+=L.get(L.n-1,k)*L.get(L.n-1,k);
    L.put(L.n-1,L.n-1,sqrt(A.get(A.n-1,A.n-1)-value));
    return make_pair(L,true);
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

Matrix SolveUsingCholesky(Matrix A, Matrix b)
{
  Matrix S1,S2; //Solución del vector de matrices
  pair<Matrix,bool> CHOL=Cholesky(A);
  if(CHOL.second)
  {
    //La matriz se pudo factorizar por Cholesky
    S1=Forward_Substitution(CHOL.first,b);
    S2=Backward_Substitution(CHOL.first.Transpose(),S1);
  }
  return S2;
}

int main(int argv, char* argc[])
{
    Matrix A(argc[1]); //Leer las dos matices desde archivos de texto
    Matrix b(argc[2]);
    Matrix S=SolveUsingCholesky(A,b);
    S.print_to_text("SOLUCIONES_BIG.txt");
    return 0;
}
