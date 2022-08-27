/* Implementación de la factorización LDL, en este caso ya que la LDU pues no tiene mucho sentido
, notemos que este es una caso especial de la LDL y esperamos que la matriz sea positiva definda 
@Author Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include<bits/stdc++.h>
#include<cmath> 
#include<fstream>
#include<sstream>
#include<string>
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

Matrix Solve_For_Diagonal(Matrix A, Matrix b)
{
    //Como se piden más cosas las vamos a guardar en matrices diferentes aunque se podrían guardar en una sola
    vector<double> solucion(A.n);
    for(int i=0;i<A.n;++i) solucion[i]=b.matrix[i]/A.matrix[i]; //Suponemos que la matriz no es 0
    Matrix X(solucion,A.n,1);
    return X;
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

Matrix ExtractDiagonal(Matrix &U)
{
    vector<double> Diagonal(U.n);
    for(int i=0;i<U.n;i++)
    {
        Diagonal[i]=U.matrix[i*U.m+i];
        for (int j=i;j<U.m;j++) U.matrix[i*U.m+j]/=Diagonal[i];
    }
    Matrix M(Diagonal,U.n,1);
    return M;
}

pair<Matrix,pair<Matrix,Matrix>> LDUFactorization(Matrix &A)
{
    //Es basicamente hacer factorización LU con pasos extra y mínima ganancia al hacerlo pero fine
    pair<pair<Matrix,Matrix>,bool> LUDEC=LUDecomposition(A);
    Matrix D=ExtractDiagonal(LUDEC.first.second);
    return make_pair(LUDEC.first.first,make_pair(D,LUDEC.first.second));
}

Matrix SolveUsingLDU(Matrix A, Matrix b)
{
    pair<Matrix,pair<Matrix,Matrix>> LDUDEC=LDUFactorization(A);
    Matrix S2,S3,S1;
    S1=Forward_Substitution(LDUDEC.first,b);
    S2=Solve_For_Diagonal(LDUDEC.second.first,S1);
    S3=Backward_Substitution(LDUDEC.second.second,S2);
    return S3;
}

int main(int argv, char* argc[])
{
    Matrix A;
    Matrix b;
    Read_from_file(argc[1],A);
    Read_from_file(argc[2],b);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    Matrix Sol=SolveUsingLDU(A,b);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Tiempo transcurrido: " << elapsed_seconds.count() << "s\n";
    Sol.print_to_text("solucionesBIG.txt");
    return 0;
}
