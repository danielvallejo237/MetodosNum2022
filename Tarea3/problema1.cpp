/*
Código de encontrar los máximos y los mínimos de una función mediante diferencias finitas

@Author Daniel Vallejo Aldana (danielvallejo237) on Github 

*/

#include<bits/stdc++.h>
#include<cmath>
#include <ostream>

using namespace std;

double finite_derivates(double x, double(*fun)(double), int degree)
{
    //Recibe un punto, un apuntador a la función de la cual queremos encontrar su derivada y el grado de la derivada
    double delta=x*0.0001+0.005; //Para evitar las divisiones por 0 cuando x es cero
    switch (degree)
    {
        case 1:
        return (fun(x+delta)-fun(x))/delta;
        break;
        case 2:
        return (fun(x+delta)-2*fun(x)+fun(x-delta))/(delta*delta);
        break;
        default:
        break;
    }
    return delta;
}

vector<pair<double,bool>> FindMaximalAndMinimal(vector<double> rango,  vector<double> fd,vector<double> sd)
{
    double tolerance=0.0001;
    vector<pair<double,bool>> points;
    for(vector<double>::iterator it=rango.begin(),it2=fd.begin(),it3=sd.begin();it!=rango.end();++it,++it2,++it3)
    {
        if(abs(*it2)<tolerance && *it3>0) points.push_back(make_pair(*it,0));
        else if (abs(*it2)<tolerance && *it3<0) points.push_back(make_pair(*it,1));
        else continue;
    }
    return points;
}
double Cubo(double x){return x*x*x*x+1;} //x^{3}+1

int main()
{
    vector<double> range(1000);
    vector<double> first_derivate(1000);
    vector<double> second_derivate(1000);
    double a=-10,b=10;
    for (int i=0;i<1000;i++)
    {
        range[i]=a+i*((b-a)/1000);
        first_derivate[i]=finite_derivates(range[i],&Cubo,1);
        second_derivate[i]=finite_derivates(range[i],&Cubo,2);
    }

    vector<pair<double,bool>> puntos;
    puntos=FindMaximalAndMinimal(range,first_derivate,second_derivate);
    cout<<"Punto"<<" "<<"Es maximo (1) o minimo (0)?"<<endl;
    for(vector<pair<double,bool>>::iterator it=puntos.begin();it!=puntos.end();++it) cout<<(*it).first<<"\t"<<(*it).second<<endl;
    //Estas salidas son para visualización
    /*ofstream os{"salida_der.txt"};
    for(vector<double>::iterator it=first_derivate.begin();it!=first_derivate.end();++it) os<<*it<<" ";
    ofstream os2{"salida_2der.txt"};
    for(vector<double>::iterator it=second_derivate.begin();it!=second_derivate.end();++it) os2<<*it<<" ";*/
    return 0;
}
