/*Implementacion del metodo de iteracion en subespacio para encontrar los valores y vectores propios mas
grandes
Usamos en este caso la iteracion de Chevichev para acelerar la convergencia tan lenta del método
@Author Daniel Vallejo Aldana / danielvallejo237 on Github*/

#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>
using namespace std;

// Usamos la misma clase matriz que hemos usado a lo largo del curso para
//poder resolver el método de jacobi

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
        return matrix[row*m+col];
    }
    void put(int row, int col, double value)
    {
        this->matrix[row*m+col]=value;
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
    void ones(int n)
    {
      this->n=n;
      this->m=n;
      vector<double> aux(n*n,1.0);
      this->matrix=aux;
    }
    void eye(int n)
    {
      this->n=n;
      this->m=n;
      vector<double> aux(n*n,0.0);
      this->matrix=aux;
      for (int i=0;i<n;i++) this->matrix[i*m+i]=1.0;
    }
    void normalize(int dim)
    {
      //We normalize the vectors along a given dimension
      if(dim==1)
      {
        for(int i=0;i<m;i++)
        {
          double suma=0;
          for(int j=0;j<n;j++) suma+=get(j,i)*get(j,i);
          suma=sqrt(suma);
          for(int j=0;j<n;j++) put(j,i,get(j,i)/suma);
        }
      }
      else
      {
        for(int i=0;i<n;i++)
        {
          double suma=0;
          for(int j=0;j<m;j++) suma+=get(i,j)*get(i,j);
          suma=sqrt(suma);
          for(int j=0;j<m;j++) put(i,j,get(i,j)/suma);
        }
      }
    }
    vector<double> ToVector()
    {
      return this->matrix;
    }
    vector<vector<double>> GetColumns()
    {
      vector<vector<double>> C(this->m,vector<double> (this->n));
      for (int j=0;j<this->m;j++) for(int i=0;i<this->n;i++) C[j][i]=get(i,j);
      return C;
    }
    vector<vector<double>> GetMColumns(int numCols)
    {
      vector<vector<double>> C(numCols,vector<double> (this->n));
      for (int j=0;j<numCols;j++) for(int i=0;i<this->n;i++) C[j][i]=get(i,j);
      return C;
    }
    void fromCols(vector<vector<double>> M)
    {
      this->matrix.resize(M[0].size()*M.size());
      this->n=M[0].size();
      this->m=M.size();
      for(int i=0;i<M.size();i++) for (int j=0;j<M[0].size();j++) put(j,i,M[i][j]);
    }
    void fromRows(vector<vector<double>> M)
    {
      this->matrix.resize(M[0].size()*M.size());
      this->n=M.size();
      this->m=M[0].size();
      for(int i=0;i<M.size();i++) for (int j=0;j<M[0].size();j++) put(i,j,M[i][j]);
    }
    Matrix deepcopy()
    {
      Matrix CP(this->matrix,n,m);
      return CP;
    }
};

pair<Matrix,Matrix> GenExamMat()
{
  Matrix A(50,50);
  Matrix b(50,1);
  A.put(0,0,40);
  A.put(0,1,-8);
  A.put(0,2,-4);
  A.put(1,0,-8);
  A.put(1,1,40);
  A.put(1,2,-8);
  A.put(1,3,-4);
  b.put(0,0,20);
  b.put(1,0,50);
  b.put(48,0,50);
  b.put(49,0,20);
  A.put(48,46,-4);
  A.put(48,47,-8);
  A.put(48,48,40);
  A.put(48,49,-8);
  A.put(49,47,-4);
  A.put(49,48,-8);
  A.put(49,49,40);
  vector<double> c={-4,-8,40,-8,-4};
  int counter;
  for(int i=2;i<48;i++)
  {
    counter=0;
    for(int j=(i-2);j<=(i+2);j++)
    {
      A.put(i,j,c[counter]);
      counter++;
    }
    b.put(i,0,100);
  }
  return make_pair(A,b);
}

int main()
{
  pair<Matrix,Matrix> P=GenExamMat();
  P.first.print_to_text("Aexam.txt");
  P.second.print_to_text("bexam.txt");
  return 0;
}
