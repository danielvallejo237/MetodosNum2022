/*
Implementación del método de Gauss-Seidel para solución de sistemas de ecuaciones
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

pair<Matrix,bool> GaussSeidelMethod(Matrix A, Matrix b, Matrix &XO, int MaxIter=200, double TOL=1e-8)
{
    Matrix aux(XO.n,XO.m); //Inicializamos un vector auxiliar
    Matrix aux2(XO.n,XO.m);
    double value1=0;
    double value2=0;
    bool Flag=true;
    int tmp=MaxIter;
    while(MaxIter--)
    {
      //Hacemos el proceso de acuerdo al número de iteraciones
      for(int i=0;i<A.n;i++)
      {
        value1=0,value2=0;
        for(int j=0;j<i;j++) value1+=A.get(i,j)*aux.get(j,0);
        for(int j=i+1;j<A.n;j++) value2+=A.get(i,j)*XO.get(j,0);
        aux.put(i,0,(1/A.get(i,i))*(-(value1+value2)+b.get(i,0)));
      }
      aux2=aux-XO;
      if(aux2.norm()/aux.norm()<TOL) break; //Encontramos la solucion
      else XO=aux;
    }
    if(MaxIter<=0) cout<<"Salida por maximo de iteraciones"<<endl,Flag=false;
    else cout<<"Iteraciones necesarias: "<<tmp-MaxIter<<endl;
    return make_pair(aux,Flag);
}

int main(int argc, char *argv[])
{
  Matrix A(argv[1]);
  Matrix b(argv[2]);
  Matrix XO(b.n,1);
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  pair<Matrix,bool> Sol=GaussSeidelMethod(A,b,XO);
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
  if (Sol.second) Sol.first.print_to_text("SalidaGS.txt");
  if (argc>3)
  {
    Matrix Comp(argv[3]);
    Matrix S=Comp-Sol.first;
    cout<<"Error de estimacion: "<<S.norm()/Comp.norm()<<endl;
  }
  return 0;
}
