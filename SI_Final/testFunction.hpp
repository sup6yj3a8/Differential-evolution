//
//  testFunction.hpp
//  SI_Final
//
//  Created by Aaron on 2021/6/9.
//

#ifndef testFunction_hpp
#define testFunction_hpp

#include <iostream>
#include <vector>
#include <cmath>

class Constraints{
public:
    Constraints(const double lower = 0.0, const double upper = 1.0, const bool isConstrained = false);
    bool check(const double candidate);
    
    inline double getLower(){return lower;};
    inline double getUpper(){return upper;};
    inline bool getIsConstrained(){return isConstrained;};
    
private:
    double lower;
    double upper;
    bool isConstrained;
};

class ObjectFunction{
public:
    ObjectFunction(const unsigned int dims = 30,
                   const double low = 0.0,
                   const double up = 1.0)
    :dimension(dims), lower(low), upper(up)
    {}
    
    std::vector<Constraints> getConstraints() const{
        std::vector<Constraints> constr(getDimension());
        for (auto& c : constr){
            c = Constraints(lower, upper, true);
        }
        return constr;
    };
    
    inline unsigned int getDimension() const {return dimension;}
    inline double getUpper()const {return upper;}
    inline double getLower()const {return lower;}
    
    virtual double evaluteCost(const std::vector<double> inputs) const{return 0;}
    virtual ~ObjectFunction(){}; // virtual destructor
    
private:
    const unsigned int dimension;
    const double lower;
    const double upper;
};

const unsigned int DIMS = 30;

/* Ackley function*/
class Ackley : public ObjectFunction{
public:
    explicit Ackley(const int dims = DIMS)
    :ObjectFunction(dims, -32.768, 32.768)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        const double a = 20;
        const double b = 0.2;
        const double c = 2 * M_PI;
        double squareSum = 0.0;
        double cosSum = 0.0;
        
        for (const auto &input : inputs){
            squareSum += input * input;
            cosSum += cos(c * input);
        }
        
        return - exp(cosSum / dimension) + exp(1) -a * exp(-b * sqrt(squareSum / dimension)) + a ;
    };
};


/* Griewank Function */
class Griewank : public ObjectFunction{
public:
    explicit Griewank(const int dims = DIMS)
    :ObjectFunction(dims, -600, 600)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double sum = 0.0;
        double product = 1;
        
        for (int i=0; i<dimension; ++i) {
            sum += inputs[i] * inputs[i] / 4000;
            product *= cos(inputs[i] / sqrt(i + 1));
        }
        
        // When the value of input is much small, sum is also much small.
        // Then, product ~= 1 -> sum - product + 1 -> sum is ignored.
        // However, expression is 1 - product + sum.
        return 1 - product + sum;
    }
    
};

/* Michalewicz function*/
class Michalewicz : public ObjectFunction{
public:
    explicit Michalewicz(const int dims = DIMS)
    :ObjectFunction(dims, 0, M_PI)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        const double m = 10;
        double sum = 0.0;
        
        for (int i=0; i<dimension; ++i){
            sum += sin(inputs[i]) * pow(sin(i * inputs[i] * inputs[i] / M_PI), 2 * m);
        }
        return -sum;
    }
};

/* Powell Function */
class Powell : public ObjectFunction{
public:
    explicit Powell(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -4, 5)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double sum = 0.0;
        
        for (int i=1; i<=dimension/4; ++i) {
            const double &x1 = inputs[4 * i - 3 - 1];
            const double &x2 = inputs[4 * i - 2 - 1];
            const double &x3 = inputs[4 * i - 1 - 1];
            const double &x4 = inputs[4 * i - 1];
                        
            sum += pow(x1 + 10 * x2, 2) + 5 * pow(x3 - x4, 2) +
                   pow(x2 - 2 * x3, 4) + 10 * pow(x1 - x4, 4);
        }
        return sum;
    }
};

/* Rastrigin function*/
class Rastrigin : public ObjectFunction{
public:
    explicit Rastrigin(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -5.12, 5.12)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        const double A = 10;
        double squareSum = 0.0;
        double sum = 0.0;
        
        for (const auto &input : inputs){
            squareSum += input * input;
            sum += A * cos(2 * M_PI * input);
        }
        return A * dimension - sum + squareSum;
    }
};

/* Rosenbrock function */
class Rosenbrock : public ObjectFunction{
public:
    explicit Rosenbrock(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -2.048, 2.048)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double sum = 0.0;
        
        for (int i=0; i<dimension-1; ++i) {
            sum += 100 * pow(inputs[i+1]-inputs[i] , 2) + pow(inputs[i] -1, 2);
        }
        
        return sum;
    }
};

/* Schwefel Function */
class Schwefel : public ObjectFunction{
public:
    explicit Schwefel(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -500, 500)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double sum = 0.0;
        
        for (int i=0; i<dimension; ++i) {
            sum += inputs[i] * sin(sqrt(abs(inputs[i])));
        }
        return 418.9829 * dimension - sum;
    }
};


/* Sphere function */
class Sphere : public ObjectFunction{
public:
    explicit Sphere(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -5.12, 5.12)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        double squareSum = 0.0;
        
        for (const auto &input : inputs){
            squareSum += input * input;
        }
        
        return squareSum;
    }
};

/* Sum Squares Function */
class SumSquares : public ObjectFunction{
public:
    explicit SumSquares(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -10, 10)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double squareSum = 0.0;
        
        for (int i=0; i<dimension; ++i) {
            squareSum += i * inputs[i] * inputs[i];
        }
        
        return squareSum;
    }
};

/* Zakharov Function */
class Zakharov : public ObjectFunction{
public:
    explicit Zakharov(const unsigned int dims = DIMS)
    :ObjectFunction(dims, -5, 10)
    {}
    
    double evaluteCost(const std::vector<double> inputs) const override{
        const int dimension = getDimension();
        double squareSum = 0.0;
        double sum1 = 0.0;
        double sum2 = 0.0;
        
        for (int i=0; i<dimension; ++i) {
            squareSum += inputs[i] * inputs[i];
            sum1 += pow(0.5 * (i + 1) * inputs[i], 2);
            sum2 += pow(0.5 * (i + 1) * inputs[i], 4);
        }
        return squareSum + sum1 + sum2;
    }
};



#endif /* testFunction_hpp */
