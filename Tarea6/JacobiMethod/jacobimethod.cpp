/*
Implementación del método de Jacobi para el cálculo de los valores y los vectores propios
@Author Daniel Vallejo Aldana / danielvallejo237 on Github
*/

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
};

double Todegrees(double radians)
{
  return radians*(180/M_PI);
}

void JacobiRotate(int p, int q,Matrix &A,Matrix &O)
{
  //We rotate a matrix in place
  double theta=0.5*atan2(2*A.get(p,q),(A.get(q,q)-A.get(p,p)));
  if(A.get(p,p)==A.get(q,q)) theta=M_PI/4;
  double s=sin(theta);
  double c=cos(theta);
  //Recives the matrix of ones and the original matrix
  //We first rotate the A matrix and then we compute the eigenvectors
  double a1,a2,a3,a4,a5;
  a1=A.get(p,p); //We will use them later
  a2=A.get(q,q);
  a3=A.get(p,q);
  for (int i=0;i<A.m;i++)
  {
    if(i!=p && i!=q)
    {
      a4=A.get(p,i); //Original values of the matrix
      a5=A.get(q,i);
      A.put(p,i,c*a4-s*a5);
      A.put(i,p,A.get(p,i)); //First change
      A.put(q,i,s*a4+c*a5);
      A.put(i,q,A.get(q,i));
    }
  }
  A.put(p,p,(c*c)*a1-2*s*c*a3+(s*s)*a2);
  A.put(q,q,(s*s)*a1+2*s*c*a3+(c*c)*a2);
  A.put(p,q,(s*c)*(a1-a2)+(c*c-s*s)*a3);
  A.put(q,p,A.get(p,q)); //This is the rotation bucle
  // We now compute the rest of the eigenvectors
  for(int i=0;i<A.n;i++)
  {
        a1=O.get(i,p);
        a2=O.get(i,q);
        O.put(i,p,c*a1-s*a2);
        O.put(i,q,s*a1+c*a2);
  }
}


pair<int,int> FindMaxOffDiagonal(Matrix &A)
{
  double maxi=fabs(A.get(0,1));
  pair<int,int> pareja(0,1);
  for(int i=0;i<A.n-1;i++)
  {
    for(int j=i+1;j<A.m;j++) if(fabs(A.get(i,j))>maxi) maxi=fabs(A.get(i,j)),pareja.first=i,pareja.second=j;
  }
  return pareja;
}

double ComputeError(Matrix &A)
{
  double sum=0;
  for(int i=0;i<A.n-1;i++)
  {
    for(int j=i+1;j<A.m;j++) sum+=fabs(A.get(i,j));
  }
  return sum;
}

pair<Matrix,Matrix> JacobiMethod(Matrix &A,double tol=1e-3, int maxiters=100)
{
  pair<int,int> toRotate=FindMaxOffDiagonal(A);
  Matrix O;
  O.eye(A.n);
  while(ComputeError(A)>tol)
  {
    JacobiRotate(toRotate.first,toRotate.second,A,O);
    toRotate=FindMaxOffDiagonal(A);
  }
  return make_pair(A,O);
}

void printRow(Matrix A,int row)
{
  for(int i=0;i<A.m;i++) cout<<A.get(i,row)<<"\t";
  cout<<endl;
}

Matrix getDiagonalOnly(Matrix A)
{
  Matrix D(A.n,1);
  for(int i=0;i<A.n;i++) D.put(i,0,A.get(i,i));
  return D;
}

double ComputeError(Matrix A,vector<double> x,double lambda)
{
  vector<double> Ax(A.n);
  vector<double> lx(A.n);
  double suma=0;
  for(int i=0;i<A.n;i++){lx[i]=lambda*x[i];for(int j=0;j<A.m;j++) Ax[i]+=A.get(i,j)*x[j];};
  for(int i=0;i<A.n;i++) lx[i]-=Ax[i],lx[i]*=lx[i],suma+=lx[i];
  return sqrt(suma);
}

vector<double> getCol(Matrix A,int i)
{
  vector<double> row(A.n);
  for(int j=0;j<A.n;j++) row[j]=A.get(j,i); 
  return row;
}
int main(int argc, char *argv[])
{
  Matrix A(argv[1]);
  Matrix ToCmp(argv[1]);
  pair<Matrix,Matrix> SOL=JacobiMethod(A);
  vector<double> errores(A.n);
  for(int i=0;i<A.n;i++)
  {
    cout<<"Valor propio: "<<endl;
    cout<<A.get(i,i)<<"\n\n";
    cout<<"Vector propio asociado"<<endl;
    printRow(SOL.second,i);
    cout<<"Error ||Ax-lx||: "<<endl;
    cout<<ComputeError(ToCmp,getCol(SOL.second,i),A.get(i,i))<<endl;
    errores[i]=ComputeError(ToCmp,getCol(SOL.second,i),A.get(i,i));
    cout<<endl;
  }
  Matrix Di=getDiagonalOnly(SOL.first);
  Matrix err(errores,A.n,1);
  Di.print_to_text("eigenvals2.txt");
  err.print_to_text("errores2.txt");
  return 0;
}
