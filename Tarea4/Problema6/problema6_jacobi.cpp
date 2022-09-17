/*
Solución del problema elíptico unidimensional usando diferencias finitas para un determinado numero de nodos
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
    Matrix operator * (const double num)
    {
      //Necesitamos definir la multiplicación por un escalar
      Matrix Sol(this->n,this->m);
      for(int i=0;i<this->n*this->m;i++) Sol.matrix[i]=this->matrix[i]*num;
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

pair<Matrix,Matrix> BuildMatrixEliptical(double stepsize,int NNodes)
{
  //Construimos una matriz cuadrada de (NNodes-2 x NNodes-2) para estimar ese número de variables
  Matrix S(NNodes-2,NNodes-2); //Inicializamos una matriz con 0's
  Matrix b(NNodes-2,1); //Estimamos el vector b  para poder calcular la solucion
  for(int i=0;i<(NNodes-2);i++)
  {
    if(i==0)
    {
      S.put(i,i,-2.0);
      S.put(i,i+1,1.0);
      b.put(i,0,2.0);
    }
    else if(i==NNodes-3)
    {
      S.put(i,i-1,1.0);
      S.put(i,i,-2.0);
      b.put(i,0,2-(2/(stepsize*stepsize)));
    }
    else
    {
      S.put(i,i,-2.0);
      S.put(i,i-1,1.0);
      S.put(i,i+1,1.0);
      b.put(i,0,2.0);
    }
  }
  S=S*(1/(stepsize*stepsize));
  return make_pair(S,b);
}


pair<Matrix,bool> JacobiMethod(Matrix A, Matrix b, Matrix &XO, int MaxIter=1000, double TOL=1e-8)
{
    Matrix aux(XO.n,XO.m); //Inicializamos un vector auxiliar
    Matrix aux2(XO.n,XO.m);
    double value=0;
    bool Flag=true;
    int tmp=MaxIter;
    while(MaxIter--)
    {
      //Hacemos el proceso de acuerdo al número de iteraciones
      for(int i=0;i<A.n;i++)
      {
        value=0;
        for(int j=0;j<A.n;j++) if(j!=i) value+=A.get(i,j)*XO.get(j,0);
        aux.put(i,0,(1/A.get(i,i))*(-value+b.get(i,0)));
      }
      aux2=aux-XO;
      if(aux2.norm()/aux.norm()<TOL) break; //Encontramos la solucion
      else XO=aux;
    }
    if(MaxIter<=0) cout<<"Salida por maximo de iteraciones"<<endl,Flag=false;
    else cout<<"Iteraciones necesarias: "<<tmp-MaxIter<<endl;
    return make_pair(aux,Flag);
}

vector<double> linspace(double low,double high, int N)
{
  //Función auxiliar de linspace para equiespaciado de los puntos
  vector<double> espaciado;
  double step=(high-low)/(double)N;
  for(int i=1;i<N-1;i++) espaciado.push_back(low+i*step);
  return espaciado;
}

Matrix EvaluateFunc(double(*fun)(double),double low, double high, int N)
{
  vector<double> linsp=linspace(low,high,N);
  Matrix M(N-2,1);
  for(int i=0;i<linsp.size();i++) M.put(i,0,fun(linsp[i]));
  return M;
}

double funcion(double x){return x*x+x;}

int main(int argc, char const *argv[])
{
  //El unico argumento de entrada en este caso es el numero de nodos en el sistema
  int Nnodos=atoi(argv[1]);
  pair<Matrix,Matrix> Sistema=BuildMatrixEliptical(1.0/(double)Nnodos,Nnodos);
  Matrix XO=Matrix(Nnodos-2,1);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  pair<Matrix,bool> Soluciones=JacobiMethod(Sistema.first,Sistema.second,XO);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
  Matrix Real=EvaluateFunc(funcion,0,1,Nnodos);
  Matrix Error=Real-Soluciones.first;
  cout<<"Error de estimacion: "<<Error.norm()/Real.norm()<<endl;
  return 0;
}
