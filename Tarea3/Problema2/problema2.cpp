/*
Codigo del problema 2, para resolver un sistema de ecuaciones de la forma Ax=b donde A es una matriz diagonal
calcular el determinante y la inversa de la función

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
};

void read_only_diagonal(Matrix &mat)
{
    vector<double> diagonal;
    //Función que solamente lee la diagonal de una matriz general
    for (int i=0;i<mat.n;i++) diagonal.push_back(mat.matrix[i*mat.m+i]);
    mat.set_vector(diagonal);
    mat.set_dimensions(mat.n,1);
    mat.set_type(2);
}

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

Matrix Solve_For_Diagonal(Matrix A, Matrix b)
{
    //Como se piden más cosas las vamos a guardar en matrices diferentes aunque se podrían guardar en una sola
    vector<double> solucion(A.n);
    for(int i=0;i<A.n;++i) solucion[i]=b.matrix[i]/A.matrix[i]; //Suponemos que la matriz no es 0
    Matrix X(solucion,A.n,1);
    return X;
}

Matrix FindInverseForDiagonal(Matrix A)
{
    //Como se piden más cosas las vamos a guardar en matrices diferentes aunque se podrían guardar en una sola
    vector<double> solucion(A.n);
    for(int i=0;i<A.n;++i) solucion[i]=1/A.matrix[i]; //Suponemos que la matriz no es 0
    Matrix X(solucion,A.n,1);
    return X;
}

double Compute_Determinant(Matrix A)
{
    double det=1.0;
    for (int i=0;i<A.n;i++) det*=A.matrix[i];
    return det;
}

int main(int argv, char *argc[])
{
    Matrix A;
    Matrix b;
    cout<<"Reading Matrix A"<<endl;
    Read_from_file(argc[1],A);
    cout<<"Reading Matrix B"<<endl;
    Read_from_file(argc[2],b);
    read_only_diagonal(A);
    Matrix p,q;
    double det;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    p=Solve_For_Diagonal(A,b);
    p.print_to_text("soluciones.txt");
    q=FindInverseForDiagonal(A);
    q.print_to_text("solucionesInversa.txt");
    det=Compute_Determinant(A);
    cout<<endl<<"Determinante de la matriz: "<<det<<endl;
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
    return 0;
}