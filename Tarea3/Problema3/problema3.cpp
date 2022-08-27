/* 
Codigo de sustitucion hacia adelante para una matriz triangular inferior

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

/*
Podemos almacenar las matrices triangulares superiores de tal forma que no considereremos las entradas con 0,
no obstante lo anterior consume tiempo O(n^2) que es lo que principalmente nos interesa
*/

void Create_Triangular_Matrix(Matrix &A, int type)
{
    //Convertir una matriz completa a una matriz triangular, el tipo nos indica si queremos una del tipo
    //Triangular superior o del tipo triangular inferior
    if(type == 1) for(int i=0;i<A.n-1;i++) for(int j=i+1;j<A.m;j++) A.matrix[i*A.m+j]=0;
    else for(int i=1;i<A.n;i++) for (int j=0;j<i;j++) A.matrix[i*A.m+j]=0; //Triangular superior

}

Matrix Forward_Substitution(Matrix A, Matrix b)
{
    vector<double> solutions(A.n);
    solutions[0]=b.matrix[0]/A.matrix[0]; //Caso base y asumiendo que jamás encontraremos una división por cero
    double value=0;
    for(int i=1;i<A.n;i++)
    {
        //Iteramos sobre todos los renglones que tenemos de la matriz
        for(int j=0; j<i;j++) value+=solutions[j]*A.matrix[i*A.m+j];
        solutions[i]=(b.matrix[i]-value)/A.matrix[i*(A.m+1)];
        value=0; //Regresamos a la inicialización del valor de valor
    }
    Matrix mat(solutions,A.n,1);
    return mat; //Regresamos la matriz de solución con sustitución hacia adelante
}

int main(int argv, char *argc[])
{
    Matrix A;
    Matrix b;
    cout<<"Reading Matrix A"<<endl;
    Read_from_file(argc[1],A);
    cout<<"Reading Matrix B"<<endl;
    Read_from_file(argc[2],b);
    Create_Triangular_Matrix(A,1);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    Matrix Sol=Forward_Substitution(A,b);
    Sol.print_to_text("soluciones.txt");
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
    return 0;
}