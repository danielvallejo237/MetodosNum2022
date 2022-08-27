/* Implementar un algoritmo para resolver un sistema de ecuaciones
mediante eliminacion gaussiana 
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/ 

#include<bits/stdc++.h>
#include<cmath> 
#include<fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <chrono>
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
        this->type=1;
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
    void set_type(int type)
    {
        this->type=type;
    }
    void print_content()
    {
        for(vector<double>::iterator it=this->matrix.begin();it!=this->matrix.end();++it) cout<<*it<<" ",
        cout<<endl;
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
};

Matrix Read_from_file(string file, Matrix &Mat)
{
    ifstream ifs(file);
    string line;
    getline(ifs, line);
    istringstream iss(line);
    int n,m;
    iss>>n>>m;
    Mat.set_dimensions(n,m);
    double a;
    while (ifs >> a) Mat.matrix.push_back(a);
    return Mat;
}

bool toRREF(Matrix &Mat, Matrix &b,Matrix &order) 
{
    //Suponemos que orden ya esta definido desde antes 
    //Vamos a mandar una matriz a su forma escalonada mediante eliminación gaussiana
    for(int k=0;k<Mat.n;k++)
    {
        //Pivoteo parcial para las matrices
        int i_max=k;
        double v_max=Mat.matrix[i_max*Mat.m+k]; //El primer elemento de la matriz
        for (int i=k+1;i<Mat.n;i++) 
        {
            if(abs(Mat.matrix[i*Mat.m+k])>v_max) v_max = Mat.matrix[i*Mat.m+k], i_max = i;
        }
        if(!Mat.matrix[k*Mat.m+i_max]) return false; //La matriz es singular
        if(i_max!=k) Mat.swap_row(i_max,k),b.swap_row(i_max,k),order.swap_row(i_max,k);
        //Hasta aquí termina la parte del pivoteo parcial
        for(int i=k+1;i<Mat.n;i++)
        {
            double f=Mat.matrix[i*Mat.m+k]/Mat.matrix[k*Mat.m+k];
            for (int j=k+1;j<Mat.m;j++) Mat.matrix[i*Mat.m+j]-=Mat.matrix[k*Mat.m+j]*f;
            b.matrix[i]-=b.matrix[k]*f;
            Mat.matrix[i*Mat.m+k]=0;
        }
    }
    return true; //Si se pudo completar llevar a forma escalonada la matriz en cuestión.
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

Matrix SolverGaussianElimination(Matrix &A, Matrix &b)
{
    Matrix solucion;
    vector<double> orden(b.n);
    for(int i=0;i<orden.size();i++) orden[i]=i;
    Matrix Orden(orden,orden.size(),1);
    bool indic=toRREF(A,b,Orden);
    if(!indic) {cout<<"La matriz es singular"<<endl;return solucion;} //Regresamos una matriz vacía
    else
    {
        solucion=Backward_Substitution(A,b);
        return solucion;
    }
}

int main(int argv, char *argc[])
{
    Matrix A;
    Matrix b;
    Read_from_file(argc[1],A);
    Read_from_file(argc[2],b);
    //Creamos el orden de las variables
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    Matrix Sol=SolverGaussianElimination(A,b);
    end = std::chrono::system_clock::now();
    Sol.print_to_text("solucionesBIG.txt");
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
    return 0;
}