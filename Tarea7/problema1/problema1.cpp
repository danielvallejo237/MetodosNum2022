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

void normalizeVector(vector<double> &v)
{
  double norm=0.0;
  for(int i=0;i<v.size();i++) norm+=v[i]*v[i];
  norm=sqrt(norm);
  for(int i=0;i<v.size();i++) v[i]/=norm; //La modificación de los vectore es in place
}

double DotProd(vector<double> v1,vector<double> v2)
{
  double dot;
  for(int i=0;i<v1.size();i++) dot+=v1[i]*v2[i];
  return dot;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++) v3[i]=v1[i]-v2[i];
    return v3;
}

vector<double> operator+(const vector<double>& v1, const vector<double>& v2)
{
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++) v3[i]=v1[i]+v2[i];
    return v3;
}

vector<double> operator * (const vector<double>& v1, double v2)
{
    vector<double> v3(v1.size());
    for(int i=0;i<v1.size();i++) v3[i]=v1[i]*v2;
    return v3;
}

double operator * (const vector<double>& v1,const vector<double>& v2)
{
  double dot=0;
  for(int i=0;i<v1.size();i++) dot+=v1[i]*v2[i];
  return dot;
}

vector<double> proyecta (const vector<double>& v1,const vector<double> &v2)
{
  //Calculamos la proyección de un vector sobre otro
  double num=v1*v2;
  double denom=v2*v2;
  vector<double> sal=v2*(num/denom);
  return sal;
}

double Norm(const vector<double>& v1)
{
  double norma=0;
  for(int i=0;i<v1.size();i++) norma+=v1[i]*v1[i];
  return sqrt(norma); //Calcular la norma del vector
}

void Orthonormalize(vector<vector<double>> &M)
{
  //Proceso de ortonormalización de gram schmidt para una matriz simétrica, deberían de
  // ser los vectores propios
  for (int i=1;i<M.size();i++)
  {
    for(int j=0;j<i;j++)
    {
      M[i]=M[i]-proyecta(M[i],M[j]);
    }
  }
  for(int i=0;i<M.size();i++)
  {
    normalizeVector(M[i]);
  }
}

vector<vector<double>> Copy(vector<vector<double>> M)
{
  vector<vector<double>> MC(M.size(),vector<double> (M[0].size()));
  for(int i=0;i<M.size();i++) for(int j=0;j<M[0].size();j++) MC[i][j]=M[i][j];
  return MC;
}

void GetJustDiagonal(Matrix &M)
{
  vector<double> D(M.n);
  for(int i=0;i<M.n;i++) D[i]=M.get(i,i);
  M.matrix=D;
  M.m=1;
}

Matrix MultiplyMat(Matrix A, Matrix B)
{
  Matrix R(A.n,B.m); //La dimensión de la nueva matriz
  double carry=0.0;
  if(A.m==B.n)
  {
    for(int i=0;i<R.n;i++)
    {
      for(int j=0;j<R.m;j++)
      {
        carry=0.0;
        for(int k=0;k<A.m;k++) carry+=A.get(i,k)*B.get(k,j);
        R.put(i,j,carry);
      }
    }
  }
  return R;
}

double sumOOD(Matrix A)
{
  double suma=0;
  for(int i=1;i<A.n-1;i++) for(int j=i+1;j<A.m;j++) suma+=fabs(A.get(i,j));
  return suma;
}

pair<Matrix,Matrix> SubspaceIteration(Matrix A,int nvals,double tol=1e-10, int Mit=10000)
{
  vector<vector<double>> V;
  V=A.GetMColumns(nvals);
  double error=10000000;
  Matrix Tmp1,Tmp2,Tmp3,evals;
  Orthonormalize(V);
  Tmp1.fromCols(V);
  while(fabs(1-error) > tol && Mit--)
  {
    Tmp2=MultiplyMat(A,Tmp1);
    V=Tmp2.GetColumns();
    Orthonormalize(V);
    Tmp2.fromCols(V);
    Tmp3=MultiplyMat(Tmp2.Transpose(),Tmp1);
    error=Tmp3.norm();
    Tmp1=Tmp2;
  }
  evals=MultiplyMat(Tmp1.Transpose(),MultiplyMat(A,Tmp1));
  GetJustDiagonal(evals);
  return make_pair(Tmp1.Transpose(),evals);
}

vector<double> GetRow(Matrix A, int row)
{
  vector<double> Row(A.m);
  for(int i=0;i<A.m;i++) Row[i]=A.get(row,i);
  return Row;
}

double ComputeError(Matrix A, Matrix Evecs, Matrix Evals)
{
  double suma=0;
  vector<double> aux;
  for(int i=0;i<Evals.n;i++)
  {
    aux=A*GetRow(Evecs,i)-GetRow(Evecs,i)*Evals.get(i,0);
    suma+=Norm(aux);
  }
  return suma/(double)Evals.n;
}

int main(int argc, char * argv[])
{
  Matrix A(argv[1]);
  pair<Matrix,Matrix> SOL;
  auto start = chrono::steady_clock::now();
  SOL=SubspaceIteration(A,atoi(argv[2]));
  auto end = chrono::steady_clock::now();
  double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
  cout << "Tiempo de ejecución en segundos: "<<elapsed_seconds<<" segs"<<endl;
  cout << "Error promedio ||Ax-lx||: "<<ComputeError(A,SOL.first,SOL.second)<<endl;
  string name(argv[1]);
  name.erase(name.end()-4,name.end());
  string n1=name+"_Eigenvals_SubspaceIt.txt";
  string n2=name+"_Eigenvecs_SubspaceIt.txt";
  SOL.first.print_to_text(n2);
  SOL.second.print_to_text(n1);
  return 0;
}
