
/*
Algoritmo de creacion de un spline de continuidad 2 conocido como spline cúbico
@Author: Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include<bits/stdc++.h>
#include <omp.h>

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

pair<vector<vector<double>>,vector<double>> BuildCubicSplineMatrix(vector<double> X, vector<double> Y)
{
  vector<vector<double>> M(X.size(),vector<double> (X.size()));
  vector<double> hs(X.size()-1,0);
  for(int i=0;i<X.size()-1;i++) hs[i]=X[i+1]-X[i];
  vector<double> b(X.size());
  M[0][0]=M[X.size()-1][X.size()-1]=1.0;
  b[0]=0;
  b[X.size()-1]=0;
  for(int i=1;i<M.size()-1;i++)
  {
    M[i][i-1]=hs[i-1];
    M[i][i]=2.0*(hs[i-1]+hs[i]);
    M[i][i+1]=hs[i];
    b[i]=(3.0/hs[i])*(Y[i+1]-Y[i])-(3/hs[i-1])*(Y[i]-Y[i-1]);
  }
  return make_pair(M,b);
}


void PrintMat(const vector<vector<double>> &vec)
{
  for(int i=0;i<vec.size();i++)
  {
    for(int j=0;j<vec[0].size();j++) cout<<vec[i][j]<<"\t";
    cout<<endl;
  }
}


double random(double range_from, double range_to)
{
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> unif(range_from,range_to);
    return unif(generator);
}


double g(double x)
{
  return 3.45*sin(x)+2.5*sqrt(exp(x))-3*x;
}

double g2(double x)
{
  return x+ (x*sin(x/2.0))/3;
}

vector<double> Linspace(double Init, double End, int N)
{
  double step=(End-Init)/(double)(N-1);
  vector<double> steps(N);
  steps[0]=Init;
  steps[N-1]=End;
  for(int i=1;i<N-1;i++) {steps[i]=Init+i*step;}
  return steps;
}

vector<double> GenerateRandomPoints(int N, double Init , double End)
{
  vector<double> X(N,0);
  X[0]=Init;
  X[N-1]=End;
  for(int i=1;i<X.size()-1;i++)
  {
    X[i]=random(Init,End);
  }
  sort(X.begin(),X.end());
  return X;
}

pair<vector<double>,vector<double>> GeneratePoints(double(*f)(double),int N,double Init, double End)
{
  vector<double> X;
  X=GenerateRandomPoints(N+1,Init,End);
  cout<<X.size()<<endl;
  vector<double> Y(X.size(),0);
  double aux;
  for(int i=0;i<X.size();i++) Y[i]=f(X[i]);
  return make_pair(X,Y);
}


void SendtoText(string file,const vector<double> &v1,const vector<double> &v2)
{
  ofstream os{file};
  for(int i=0;i<v1.size();i++) os<<v1[i]<<","<<v2[i]<<endl;
}

void Vec2Text(string file,const vector<double> &v1)
{
  ofstream os{file};
  for(int i=0;i<v1.size();i++) os<<v1[i]<<endl;
}
void sendMat2Text(string file, const vector<vector<double>> &M)
{
  ofstream os{file};
  os<<M.size()<<" "<<M[0].size()<<endl;
  for(int i=0;i<M.size();i++)
  {
    for(int j=0;j<M[0].size();j++) os<<M[i][j]<<" ";
    os<<endl;
  }
}

Matrix VVtoMat(vector<vector<double>> M)
{
  Matrix A(M.size(),M[0].size());
  #pragma omp parallel for
  for(int i=0;i<M.size();i++)
  {
    for(int j=0;j<M[0].size();j++) A.put(i,j,M[i][j]);
  }
  return A;
}

vector<vector<double>> ComputeABCDMat(vector<double> C,vector<double> X,vector<double> Y)
{
  vector<vector<double>> Cf(C.size(),vector<double>(4));
  int n=(int)X.size()-1;
  Cf[n][2]=C[n];
  Cf[n][0]=Y[n];
  vector<double> hs(n);
  for(int i=0;i<hs.size();i++) hs[i]=X[i+1]-X[i];
  for(int i=n-1;i>=0;i--)
  {
    Cf[i][2]=C[i];
    Cf[i][1]=(Y[i+1]-Y[i])/hs[i]-hs[i]*(C[i+1]+2*C[i])/3.0;
    Cf[i][3]=(C[i+1]-C[i])/(3*hs[i]);
    Cf[i][0]=Y[i];
  }
  return Cf;
}

double InterpolateValue(double val, vector<vector<double>> Cf,vector<double> X)
{
  int ind=0;
  for(int i=0;i<X.size();i++)
  {
    if (X[i]<=val) continue;
    else
    {
      ind=i-1;
      break;}
  }
  return Cf[ind][0]+Cf[ind][1]*(val-X[ind])+Cf[ind][2]*(val-X[ind])*(val-X[ind])+Cf[ind][3]*(val-X[ind])*(val-X[ind])*(val-X[ind]);
}

pair<vector<double>,vector<double>> GenerateInterpol(vector<double> Xs, vector<vector<double>> Cf, int N,double Init, double End)
{
  vector<double> pts=Linspace(Init,End,N);
  vector<double> ys(pts.size(),0);
  for(int i=0;i<ys.size()-1;i++) ys[i]=InterpolateValue(pts[i],Cf,Xs);
  return make_pair(pts,ys);
}


double ComputeError(double(*f)(double),int points,double Init, double End, vector<vector<double>> Cf, vector<double> X)
{
  double suma=0;
  double a1,a2;
  double val;
  for(int i=0;i<points;i++)
  {
    val=random(Init,End);
    a1=InterpolateValue(val,Cf,X);
    a2=f(val);
    suma+=fabs(a2-a1);
  }
  return suma/points;
}

int main(int argc, char *argv[])
{
  pair<vector<double>,vector<double>> P=GeneratePoints(g2,atoi(argv[1]),atof(argv[2]),atof(argv[3]));
  pair<vector<vector<double>>,vector<double>> SplineMat=BuildCubicSplineMatrix(P.first,P.second);
  Matrix A;
  A=VVtoMat(SplineMat.first);
  pair<pair<Matrix,Matrix>,bool> LU=LUDecomposition(A);
  vector<double> C=SolveLU(LU.first.first,LU.first.second,SplineMat.second);
  vector<vector<double>> Coefs=ComputeABCDMat(C,P.first,P.second);
  pair<vector<double>,vector<double>> I=GenerateInterpol(P.first,Coefs,10000,atof(argv[2]),atof(argv[3]));
  SendtoText("quadpuntos.txt",I.first,I.second);
  SendtoText("interpolpoints.txt",P.first,P.second);
  double error=ComputeError(g2,30,atof(argv[2]),atof(argv[3]),Coefs,P.first);
  cout<<"Error |f(x)-fhat(h)| con 30 puntos aleatorios: "<<error<<endl;
  return 0;
}
