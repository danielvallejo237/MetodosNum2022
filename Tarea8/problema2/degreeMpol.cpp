#include <omp.h>
#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>
#include<chrono>

#define MAX_NUM_THREADS 3 // Podemos cambiar el numero de hilos de acuerdo a la computadora pero en este caso usaremos solo 2

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
        #pragma omp parallel for
        for(int i=0;i<this->m;i++) swap(matrix[this->m*a+i],this->matrix[this->m*b+i]);
    }
    Matrix Transpose()
    {
        vector<double> T(this->n*this->m);
        #pragma omp parallel for
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
          #pragma omp parallel for
          //Tienen que tener las mismas dimensiones
          for(int i=0;i<this->n*this->m;i++) Sol.matrix[i]=this->matrix[i]+obj.matrix[i];
      }
      return Sol;
    }
    Matrix operator * (const double num)
    {
      //Necesitamos definir la multiplicación por un escalar
      Matrix Sol(this->n,this->m);
      #pragma omp parallel for
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
      #pragma omp parallel for
      for (int j=0;j<this->m;j++) for(int i=0;i<this->n;i++) C[j][i]=get(i,j);
      return C;
    }
    void fromCols(vector<vector<double>> M)
    {
      this->matrix.resize(M[0].size()*M.size());
      this->n=M[0].size();
      this->m=M.size();
      #pragma omp parallel for
      for(int i=0;i<M.size();i++) for (int j=0;j<M[0].size();j++) put(j,i,M[i][j]);
    }
    void fromRows(vector<vector<double>> M)
    {
      this->matrix.resize(M[0].size()*M.size());
      this->n=M.size();
      this->m=M[0].size();
      #pragma omp parallel for
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

double Line(double x)
{
    return 14*x+4; //La linea que queremos interpolar
}

double Cube(double x)
{
  return 3.23*(x*x*x)-0.25*(x*x)-10*x +12.47;
}

vector<double> Linspace(double initRange, double endRange, int N)
{
  double Stepsize=(endRange-initRange)/(double)N;
  vector<double> points(N+1);
  #pragma omp parallel for
  for (int i=0;i<N+1;i++) points[i]=initRange+i*Stepsize;
  return points;
}

pair<vector<double>,vector<double>> GenerateFunctionPairs(double(*f)(double),double InitRange, double EndRange, int Npoints)
{
  vector<double> ls=Linspace(InitRange,EndRange,Npoints);
  vector<double> fs(ls.size());
  #pragma omp parallel for
  for(int i=0;i<Npoints+1;i++)fs[i]=f(ls[i]);
  return make_pair(ls,fs);
}

pair<Matrix,Matrix> BuildPolinomialMatrix(vector<double> X, vector<double> Y,int degree)
{
  Matrix A(X.size(),degree+1);
  Matrix b(X.size(),1);
  for(int i=0;i<X.size();i++)
  {
    for(int j=degree;j>=0;j--)
    {
      A.put(i,(degree-j),pow(X[i],j));
    }
    b.put(i,0,Y[i]);
  }
  return make_pair(A,b);
}

double Interpolate(double x,vector<double> CF)
{
  double res=0;
  for(int i=0;i<CF.size();i++) res+=CF[i]*pow(x,i);
  return res;
}

Matrix InterpolateRange(double Init, double End, int Npts, vector<double> CF)
{
  vector<double> lsp=Linspace(Init,End,Npts);
  Matrix R(lsp.size(),2);
  for(int i=0;i<lsp.size();i++)
  {
    R.put(i,0,lsp[i]);
    R.put(i,1,Interpolate(lsp[i],CF));
  }
  return R;
}

vector<double> EvalFN(double(*f)(double),vector<double> X)
{
  vector<double> F(X.size());
  for (int i=0;i<X.size();i++)
  {
    F[i]=f(X[i]);
  }
  return F;
}

double Constant(double x)
{
  return 3.45*sin(x)+2.5*sqrt(exp(x))-3*x;
}

double func(double x)
{
  return 345*pow(x,4)+235;
}
int main(int argc, char *argv[])
{
  pair<vector<double>,vector<double>> P=GenerateFunctionPairs(&Constant,atof(argv[1]),atof(argv[2]),atoi(argv[3]));
  double degree=atoi(argv[3]);
  pair<Matrix,Matrix> MT=BuildPolinomialMatrix(P.first,P.second,degree);
  Matrix Sol=SolveQR(MT.first,MT.second);
  vector<double> Sl=Sol.ToVector();
  reverse(Sl.begin(),Sl.end());
  cout<<Sl.size()<<endl;
  Matrix O=InterpolateRange(atof(argv[1]),atof(argv[2]),atoi(argv[4]),Sl);
  O.print_to_text("Interpolacion_polinomios_"+string(argv[4])+"_puntos_de_"+string(argv[3])+"_grados.txt");
  vector<vector<double>> Cols=O.GetColumns();
  vector<double> FE=EvalFN(&Constant,Cols[0]);
  cout<<"||f-Pm(X)||: "<<Norm(Cols[1]-FE)/atof(argv[4])<<endl;
  return 0;
}
