/*
Programación del método de interpolación por elemento finito
*/

#include<bits/stdc++.h>
#include<omp.h>

#define MAX_NUM_THREADS 3

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


vector<double> Linspace(double a, double b, int N)
{
  vector<double> Points(N+1);
  double stepsize=(b-a)/(double)N;
  Points[0]=a;
  Points[N]=b;
  #pragma omp parallel for
  for(int i=1;i<N;i++)
  {
    Points[i]=a+i*stepsize;
  }
  return Points;
}

pair<vector<double>,vector<double>> GeneraPuntos(double(*f)(double),double a, double b, int N)
{
  vector<double> X=Linspace(a,b,N);
  vector<double> Y(X.size());
  #pragma omp parallel for
  for(int i=0;i<X.size();i++) Y[i]=f(X[i]);
  return make_pair(X,Y);
}

double N1(double x, double z1, double z2, double h)
{
  if(z1<=x && x <= z2)
        return 1.0-(1.0/h)*(x-z1);
    else
        return 0.0;
}

double N2(double x,double z1,double z2,double h)
{
    if(z1<=x && x <= z2)
        return (1.0/h)*(x-z1);
    else
        return 0.0;
}

vector<double> findPhi(vector<double> X,vector<double> Y,double start, double end, int N, double lambda)
{
  double h=(end-start)/(double)N;
  vector<double> z=Linspace(start,end,N);
  vector<double> a(N);
  vector<double> b(N);
  vector<double> c(N);
  vector<double> s(N);
  vector<double> t(N);
  vector<double> w(N+1);
  vector<double> phi;
  Matrix A(N+1,N+1);
  for(int i=0;i<N;i++)
  {
      double sum1=0,sum2=0,sum3=0,sum4=0,sum5=0;
      for(int j=0;j<X.size();j++)
      {
          sum1 += N1(X[j],z[i],z[i+1],h)*N1(X[j],z[i],z[i+1],h);
          sum2 += N1(X[j],z[i],z[i+1],h)*N2(X[j],z[i],z[i+1],h);
          sum3 += N2(X[j],z[i],z[i+1],h)*N2(X[j],z[i],z[i+1],h);
          sum4 += Y[j]*N1(X[j],z[i],z[i+1],h);
          sum5 += Y[j]*N2(X[j],z[i],z[i+1],h);
      }
      a[i] = sum1 + lambda/h;
      b[i] = sum2 - lambda/h;
      c[i] = sum3 + lambda/h;
      s[i] = sum4;
      t[i] = sum5;
    }
    A.put(0,0,a[0]);
    A.put(N,N,c[N-1]);
    for(int i=1;i<N;i++) A.put(i,i,c[i-1]+a[i]);
    for(int i=0;i<N;i++)A.put(i,i+1,b[i]);
    for(int i=1;i<N+1;i++) A.put(i,i-1,b[i-1]);
    for(int i=1;i<N;i++)w[i]=t[i-1] + s[i];
    w[0] = s[0];
    w[N] = t[N-1];
    pair<pair<Matrix,Matrix>,bool> LU=LUDecompositionBandedMatrix(A,3);
    phi=SolveLU(LU.first.first,LU.first.second,w);
    return phi;
}

pair<vector<double>,vector<double>> Interpolate(double a, double b,vector<double> phi)
{
  vector<double> z=Linspace(a,b,phi.size()-1);
  double h=(b-a)/(double)(phi.size()-1);
  vector<double> x=Linspace(a,b,1000);
  vector<double> y_aprox(x.size());
  for(int i=0;i<x.size();i++)
  {
        for(int j=0;j<z.size();j++){
            if(z[j]<= x[i] && x[i] <= z[j])
                y_aprox[i] = phi[j]*N1(x[i],z[j],z[j+1],h) + phi[j+1]*N2(x[i],z[j],z[j+1],h);
        }
    }
  return make_pair(x,y_aprox);
}

double g2(double x)
{
  return x+ (x*sin(x/2.0))/3;
}

void PintPairVec(string file,pair<vector<double>,vector<double>> pt)
{
  ofstream os{file};
  for(int i=0;i<pt.first.size();i++)
  {
    os<<pt.first[i]<<","<<pt.second[i]<<endl;
  }
}

double ComputeError(double(*f)(double),vector<double> X, vector<double> Y)
{
  double suma=0;
  for(int i=0;i<X.size();i++)
  {
    suma+=fabs(f(X[i])-Y[i])*fabs(f(X[i])-Y[i]);
  }
  return (1/(double)X.size())*sqrt(suma);
}


int main(int argc, char *argv[])
{
  if (argc<5) exit(1);
  omp_set_num_threads(MAX_NUM_THREADS);
  double a=atof(argv[1]);
  double b=atof(argv[2]);
  double M=atof(argv[3]);
  double N=atof(argv[4]);
  pair<vector<double>,vector<double>> P=GeneraPuntos(g2,a,b,M);
  vector<double> phi=findPhi(P.first,P.second,a,b,N,0.9);
  pair<vector<double>,vector<double>> I=Interpolate(a,b,phi);
  PintPairVec("ipol.txt",I);
  cout<<"MSELoss: "<<ComputeError(g2,I.first,I.second)<<endl;
  return 0;
}
