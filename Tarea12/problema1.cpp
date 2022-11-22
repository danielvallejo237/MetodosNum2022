/*
Programación del algoritmo de cuadratura gaussiana para integración numérica
@Author: Daniel Vallejo Aldana (danielvallejo237) on Github
*/

#include <vector>
#include <iostream>
#include <cstdlib>
#include "./fparser/fparser.hh"

using namespace std;

class LegendrePolynomial 
{
public:
    LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations)
        : mLowerBound(lowerBound), mUpperBound(upperBound), mNumberOfIterations(numberOfIterations), mWeight(numberOfIterations+1), mRoot(numberOfIterations+1) {
        calculateWeightAndRoot();
    }

    const std::vector<double> & getWeight() const {
        return mWeight;
    }

    const std::vector<double> & getRoot() const {
        return mRoot;
    }

private:
    const static double EPSILON;

    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}
    };

    void calculateWeightAndRoot() {
        for(int step = 0; step <= mNumberOfIterations; step++) {
            double root = cos(M_PI * (step-0.25)/(mNumberOfIterations+0.5));
            Result result = calculatePolynomialValueAndDerivative(root);

            double newtonRaphsonRatio;
            do {
                newtonRaphsonRatio = result.value/result.derivative;
                root -= newtonRaphsonRatio;
                result = calculatePolynomialValueAndDerivative(root);
            } while (fabs(newtonRaphsonRatio) > EPSILON);

            mRoot[step] = root;
            mWeight[step] = 2.0/((1-root*root)*result.derivative*result.derivative);
        }
    }

    Result calculatePolynomialValueAndDerivative(double x) {
        Result result(x, 0);

        double value_minus_1 = 1;
        const double f = 1/(x*x-1);
        for(int step = 2; step <= mNumberOfIterations; step++) {
            const double value = ((2*step-1)*x*result.value-(step-1)*value_minus_1)/step;
            result.derivative = step*f*(x*value - result.value);

            value_minus_1 = result.value;
            result.value = value;
        }

        return result;
    }

    const double mLowerBound;
    const double mUpperBound;
    const int mNumberOfIterations;
    std::vector<double> mWeight;
    std::vector<double> mRoot;
};

const double LegendrePolynomial::EPSILON = 1e-15;

double gaussLegendreIntegral(double a, double b, int n, FunctionParser fp) {
    const LegendrePolynomial legendrePolynomial(a, b, n);
    const std::vector<double> & weight = legendrePolynomial.getWeight();
    const std::vector<double> & root = legendrePolynomial.getRoot();

    const double width = 0.5*(b-a);
    const double mean = 0.5*(a+b);

    double gaussLegendre = 0;
    for(int step = 1; step <= n; step++) {
        double numero=width * root[step] + mean;
        gaussLegendre += weight[step]*fp.Eval(&numero);
    }
    return gaussLegendre * width;
}

int main(int argc, char *argv[])
{
    string expresion; //Expresion que contiene la función que deberá ser evaluada
    expresion=argv[1];
    double a,b;
    FunctionParser fp;  
    fp.AddConstant("pi",3.1415926535897932);
    fp.AddConstant("e",2.718281828459); //Valores más comunes encontrados en matemáticas
    fp.Parse(expresion,"x");
    a=atof(argv[2]);
    b=atof(argv[3]);
    double integral;
    integral=gaussLegendreIntegral(a,b,1000,fp);
    cout<<"Valor de la integral con Gauss-Legendre: "<<integral<<endl;
    return 0;
}