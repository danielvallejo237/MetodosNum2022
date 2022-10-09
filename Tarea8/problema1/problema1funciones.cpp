#include <omp.h>
#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>
 // Podemos cambiar el numero de hilos de acuerdo a la computadora pero en este caso usaremos solo 2

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

vector<vector<double>> ComputeR(vector<vector<double>> Basis, vector<vector<double>> Obasis)
{
  vector<vector<double>> R(Obasis[0].size(),vector<double> (Obasis.size()));
  for(int i=0;i<R.size();i++) for(int j=0;j<=i;j++) R[j][i]=Obasis[j]*Basis[i];
  return R;
}

pair<Matrix,Matrix> QR(Matrix A)
{
  vector<vector<double>> B,OB,RM;
  OB=A.GetColumns();
  B=Copy(OB);
  Orthonormalize(OB);
  RM=ComputeR(B,OB);
  Matrix Q,R;
  Q.fromCols(OB);
  R.fromRows(RM);
  return make_pair(Q,R);
}

Matrix SolveQR(Matrix A, Matrix b)
{
  pair<Matrix,Matrix> QRDEC=QR(A);
  vector<double> tmp;
  tmp=QRDEC.first.Transpose()*b.ToVector();
  Matrix s(tmp,A.n,1);
  Matrix Sol;
  Sol=Backward_Substitution(QRDEC.second,s);
  return Sol;
}

double g1(double x)
{
    return cos(x); //La linea que queremos interpolar
}

double g2(double x)
{
  return sqrt(7*exp(x));
}

double g3(double x)
{
  return x*x;
}

double Combine(double(*g1)(double),double(*g2)(double),double(*g3)(double),double  x)
{
  return 3.45*g1(x)+2.5*g2(x)-3*g3(x);
}


vector<double> Linspace(double initRange, double endRange, int N)
{
  double Stepsize=(endRange-initRange)/(double)N;
  vector<double> points(N+1);
  for (int i=0;i<N+1;i++) points[i]=initRange+i*Stepsize;
  return points;
}

pair<vector<double>,vector<double>> GenerateFunctionPairsModified(double(*f)(double),double InitRange, double EndRange, int Npoints)
{
  vector<double> ls=Linspace(InitRange,EndRange,Npoints);
  vector<double> fs(ls.size());
  for(int i=0;i<Npoints+1;i++)fs[i]=f(ls[i]);
  return make_pair(ls,fs);
}

pair<pair<Matrix,Matrix>,bool> BuildLeastSquareMatrixModified(vector<double> X, vector<double> Y,double(*g1)(double),double(*g2)(double),double(*g3)(double))
{
  bool flag=false;
  Matrix A(X.size(),3);
  Matrix B(X.size(),1);
  for (int i=0;i<X.size();i++)
  {
    B.put(i,0,Y[i]);
    A.put(i,0,g1(X[i]));
    A.put(i,1,g2(X[i]));
    A.put(i,2,g3(X[i]));
  }
  return make_pair(make_pair(A,B),true);
}

Matrix operator * (Matrix A,Matrix B)
{
  Matrix C(A.n,B.m); //The shape of the new matrix
  double suma=0;
  for(int i=0;i<C.n;i++)
  {
    for(int j=0;j<C.m;j++)
    {
      suma=0;
      for(int k=0;k<A.m;k++) suma+=A.get(i,k)*B.get(k,j);
      C.put(i,j,suma);
    }
  }
  return C;
}
vector<double> SolveUsingLSModified(double(*f)(double) ,double InitRange, double EndRange, int Npoints)
{
  pair<vector<double>,vector<double>> P=GenerateFunctionPairsModified(f,InitRange,EndRange,Npoints);
  pair<pair<Matrix,Matrix>,bool> LS=BuildLeastSquareMatrixModified(P.first,P.second,&g1,&g2,&g3);
  Matrix A;
  Matrix B;
  Matrix C;
  if(LS.second)
  {
    //Si el Booleano no es falso
    A=LS.first.first.Transpose()*LS.first.first;
    B=LS.first.first.Transpose()*LS.first.second;
    C=SolveQR(A,B);
  }
  vector<double> Sol;
  Sol=C.ToVector();
  return Sol;
}

double ReconstructFun(double x, vector<double> A, double(*g1)(double), double(*g2)(double),double(*g3)(double))
{
  double res;
  res=A[0]*g1(x)+A[1]*g2(x)+A[2]*g3(x);
  return res;
}

Matrix ToInterpolate(double Init, double End, vector<double> coef,double(*g1)(double), double(*g2)(double),double(*g3)(double))
{
  double StepSize=(End-Init)/(double)(1000-1);
  Matrix S(1000,2);
  for(int i=0;i<1000;i++)
  {
    S.put(i,0,Init+i*StepSize);
    S.put(i,1,ReconstructFun(Init+i*StepSize,coef,g1,g2,g3));
  }
  return S;
}

vector<double> EvalFNCombine(double(*g1)(double),double(*g2)(double),double(*g3)(double),vector<double> X)
{
  vector<double> F(X.size());
  for (int i=0;i<X.size();i++)
  {
    F[i]=Combine(g1,g2,g3,X[i]);
  }
  return F;
}

double Constant(double x)
{
  return 3.45*sin(x)+2.5*sqrt(exp(x))-3*x;
}

int main(int argc, char * argv[])
{
  //definimos el numero de hilos que vamos a usar
  //Lo anterior es para paralelización de los métodos
  vector<double> Coef=SolveUsingLSModified(&Constant,atof(argv[1]),atof(argv[2]),atoi(argv[3]));
  cout<<"Coeficientes encontrados"<<endl;
  for (auto c: Coef) cout<<c<<" ";
  cout<<endl;
  Matrix P;
  P=ToInterpolate(atof(argv[1]),atof(argv[2]),Coef,&g1,&g2,&g3);
  vector<vector<double>> Cols=P.GetColumns();
  vector<double> FE=EvalFNCombine(&g1,&g2,&g3,Cols[0]);
  cout<<"||f-fhat||: "<<Norm(Cols[1]-FE)/atof(argv[3])<<endl;
  P.print_to_text("SalidaInterpolFunctions"+string(argv[3])+"_puntos.txt");
  return 0;
}
