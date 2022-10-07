/* Programación de una clase Matriz de tipo banda en donde le daremos solamente el ancho
de banda a la matriz, en este caso si usaremos dobles apuntadores ya que resultan de forma más eficiente
que el hecho de usar solamente apuntadores simples

@Author Daniel Vallejo Aldana (danielvallejo237 on Github)
This code was written by danielvallejo237 for the Numerical Analysisi course (Winter 2022)
*/

#include<bits/stdc++.h>
#include<fstream>
#include<cmath>
#include<string>

using namespace std;

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
};

class BMat
{
public:
  vector<vector<double>> content;
  int n; //We assume we are working with square matrices
  int b; //Bandwidth
  BMat()
  {
    this->n=1;
    this->b=1;
    this->content.resize(1);
    this->content[0].resize(1);
    this->content[0][0]=1.0;
  }
  BMat(int n, int b)
  {
    vector<vector<double>> M(n,vector<double>(2*b-1));
    this->content=M; //Initializing from 0's
    this->n=n;
    this->b=b;
  }
  BMat(Matrix A,int b)
  {
    //We initialize the class using the exiting class Matrix
    vector<vector<double>> M(A.n,vector<double>(2*b-1,0.0));
    int counter;
    for(int i=0;i<A.n;i++)
    {
      counter=(2*b-1)/2-b/2;
      if(i+b/2>A.n-1) counter=i+(2*b-1)/2-A.n; //We start the counter some steps ahead
      else if(i-b/2<0) counter=(2*b-1)/2-i;
      for(int j=max(0,i-b/2);j<=min(A.n,i+b/2);j++)
      {
        M[i][counter]=A.get(i,j);
        counter++;
      }
    }
    this->n=A.n;
    this->b=b;
    this->content=M;
  }
  ~BMat(){}
  void printContent()
  {
    for(int i=0;i<this->n;i++){
      for(int j=0;j<(2*this->b-1);j++) cout<<this->content[i][j]<<"\t";
      cout<<endl;
    }
  }
  pair<int,int> convertCoordinates(int i,int j)
  {
    //Converting coordinates of a normal matrix to the newely created one
    return make_pair(i,min(max(0,b/2-(i-j)),this->n));
  }
  vector<double> ForwardSubsitute(vector<double> B)
  {
    vector<double> solutions(this->n,0);
    pair<int,int> coords;
    double suma;
    if (B.size()!= this->n) return solutions;
    coords=convertCoordinates(0,0);
    solutions[0]=B[0]/this->content[coords.first][coords.second];
    for(int i=1;i<this->n;i++)
    {
      suma=0;
      for(int j=max(0,i-this->b/2);j<i;j++)
      {
        coords=convertCoordinates(i,j);
        suma+=this->content[coords.first][coords.second]*solutions[j];
      }
      coords=convertCoordinates(i,i);
      solutions[i]=(B[i]-suma)/this->content[coords.first][coords.second];
    }
    return solutions;
  }
  vector<double> BackwardSubstitute(vector<double> B)
  {
    vector<double> solutions(this->n,0);
    pair<int,int> coords;
    double suma;
    int N=this->n;
    if(B.size()!=this->n) return solutions; //We have not found solutions for this problem
    coords=convertCoordinates(N-1,N-1);
    solutions[N-1]=B[N-1]/this->content[coords.first][coords.second];
    for(int i=N-2;i>=0;i--)
    {
      suma=0;
      for(int j=i+1;j<N;j++)
      {
        coords=convertCoordinates(i,j);
        suma+=solutions[j]*this->content[coords.first][coords.second];
      }
      coords=convertCoordinates(i,i);
      solutions[i]=(B[i]-suma)/this->content[coords.first][coords.second];
    }
    return solutions;
  }
  pair<BMat,BMat> LU()
  {
    /*
    Decomposing the matrix into an LU matrix where L is a lower triangular an U
    an upper triangular matrix
    */
    BMat L(this->n,this->b);
    BMat U(this->n,this->b); //Son matrices del mismo tamaño
    for(int i=0;i<this->n;i++) L.content[i][b/2]=1.0; //We have in the diagonal just ones
    U.content[0][b/2]=this->content[0][b/2];
    pair<int,int> coords;
    for(int i=1;i<=b/2;i++)
    {
      coords=convertCoordinates(i,0);
      cout<<coords.first<<" "<<coords.second<<endl;
      L.content[coords.first][coords.second]=this->content[coords.first][coords.second]/U.content[0][b/2];
    }
    //The LU factorization
    double suma;
    for(int i=1;i<n-1;i++)
    {
      suma=0;
      for(int j=0;j<b/2;j++)
      {
        suma+=U.content[j][i]*L.content[i][j];
      }
      coords=convertCoordinates(i,i);
      U.content[coords.first][coords.second]=this->content[coords.first][coords.second]-suma;
    }
    return make_pair(L,U);
  }
};

int main(int argv,char* argc[])
{
  Matrix A(argc[1]);
  BMat B(A,atoi(argc[2]));
  B.printContent();
  cout<<"\n\n";
  pair<BMat,BMat> LUDEC=B.LU();
  LUDEC.second.printContent();
  cout<<"\n\n";
  LUDEC.first.printContent();
  return 0;
}
