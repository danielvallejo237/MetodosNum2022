/*
Implementación del método de Doolittle para resolver sistemas de ecuaciones, 
el desgloce de las cuentas viene en el reporte entregado en conjunto con esta tarea

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

pair<pair<Matrix,Matrix>,bool> LUDoolittle(Matrix &A)
{
    /*Recibimos una matriz y la descomponemos en dos matrices que satisfacen que una es triangular superior
     y la otra es triangular inferior*/
    vector<double> l1(A.n*A.m,0),u1(A.n*A.m,0);
    Matrix L(l1,A.n,A.m), U(u1,A.n,A.m);
    for(int i=0;i<A.n;i++)
    {
        for(int k=i;k<A.m;k++)
        {
            //Esto es para la triangular superior
            double value=0;
            for(int j=0;j<i;j++) value+=L.matrix[i*L.m+j]*U.matrix[j*U.m+k];
            U.matrix[i*L.m+k]=A.matrix[i*A.m+k]-value;
        }
        //A continuación tenemos el ciclo para la triangular inferior
        for(int k=i;k<A.m;k++)
        {
            if(i==k) L.matrix[i*L.m+k]=1; //1s en las diagonales
            else 
            {
                double value=0;
                for (int j=0;j<i;j++) value+=L.matrix[k*A.m+j] * U.matrix[j*A.m+i];
                L.matrix[k*A.m+i]=(A.matrix[k*A.m+i]-value)/U.matrix[i*A.m+i];
            }
        }
    }
    return make_pair(make_pair(L,U),true);
}

Matrix SolveUsingDoolitlle(Matrix &A,Matrix &b)
{
    Matrix S2;
    pair<pair<Matrix,Matrix>,bool> LUDEC=LUDoolittle(A);
    if(LUDEC.second)
    {
        Matrix S1=Forward_Substitution(LUDEC.first.first,b);
        S2=Backward_Substitution(LUDEC.first.second,S1);
    }
    return S2;
}
int main(int argv, char* argc[])
{
    Matrix A;
    Matrix b;
    Read_from_file(argc[1],A);
    Read_from_file(argc[2],b);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    Matrix Sol=SolveUsingDoolitlle(A,b);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
    Sol.print_to_text("solucionesBIG.txt");
    return 0;
}