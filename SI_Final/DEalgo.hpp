//
//  DEalgo.hpp
//  SI_Final
//
//  Created by Aaron on 2021/6/9.
//

#ifndef DEalgo_hpp
#define DEalgo_hpp

#include <vector>
#include <random>
#include <limits>
#include <algorithm>
#include "testFunction.hpp"

class Point{
public:
    // Default constructor
    Point(const int dim)
    :factor(0.8), crossingRate(0.9),vec(dim),
    minCost(std::numeric_limits<double>::max()), rank(0)
    {}
    
    // Overloading copy assignment operator (=)
    Point operator= (const Point &b){
        vec = b.vec;
        minCost = b.minCost;
        rank = b.rank;
        return *this;
    }
    
    // Overloading operator[]
    inline double& operator[] (const int &i){return vec[i];}
    
    // Overloading relationship operators
    inline int operator< (const Point &b)const{return this->minCost < b.minCost;}
    inline int operator> (const Point &b)const{return this->minCost > b.minCost;}
    inline int operator== (const Point &b)const{return this->minCost == b.minCost;}
    
    // Data member
    double factor;
    double crossingRate;
    std::vector<double> vec;
    double minCost;
    int rank;
};

// APTDE status monitor
class Status{
public:
    // Data member
    int RM; // Redundance Monitor
    int NM; // Stagnation Monitor
    int UM; // Upper bound monitor
    int LM; // Lower bound monitor
    
    // Default constructor
    Status():RM(0), UM(0), NM(0), LM(0){}
};

// EPSDE, MPEDE
class Strategy{
public:
    // Default constructor
    Strategy(const double value = 0.5)
    :meanF(value), meanCR(value),
    succF(0), succCR(0),
    improve(0), fes(0),
    mutation(0)
    {}
    
    // Data member
    double meanF;
    double meanCR;
    std::vector<double> succF;
    std::vector<double> succCR;
    double improve; // The improvement of fitness value
    int fes; // Function evaluations
    int mutation; // The mutation strategy
};


// Basic DE (DE/rand/bin/1)
class DifferentialEvolution{
public:
    explicit DifferentialEvolution(const ObjectFunction &costFunction,
                                   const unsigned int &argPopulationSize,
                                   const bool &isCheckConstraints);
    void initPopulation();
    virtual void selectionAndCrossing(const int &maxGen,const int &currGen);
    void crossover(std::vector<double> &trial, const std::vector<double> &mutant, const int &i);
    
    Point getBestPoint(){return population[bestIndexPointer];}
    double getBestCost(){return bestIndexPointer;}
    
    void printPopulation();
    void optimize(int generation, int detail = 0);
    
    // Sort rule
    static int cmp(const Point &x,const Point &y){return x < y;}
    static int cmp2(const int &x,const int &y){return x < y;}
    
    virtual ~DifferentialEvolution(){}; // virtual destructor
    
protected:
    /* data member */
    // Basic DE coefficient
    const ObjectFunction& cost;
    unsigned int populationSize;
    const unsigned int DIMENSION;
    
    std::default_random_engine generator;
    std::vector<Point> population;
    std::vector<Constraints> constraints;
    std::vector<int> popIndex;
    
    // Global optimal
    int bestIndexPointer;
    double minCost;
    
    // Max and min number
    const double defaultLower = -std::numeric_limits<double>::infinity();
    const double defaultUpper = std::numeric_limits<double>::infinity();
    
    // Constraints
    bool isCheckConstraints;
    bool checkConstraints(std::vector<double> point);
};

// APTDE
class APTDE : public DifferentialEvolution {
public:
    APTDE(const ObjectFunction &costFunction,
          const unsigned int &argPopulationSize,
          const bool &isCheckConstraints);
    
    void selectionAndCrossing(const int &maxGen,const int &currGen) override;
    
    // APTDE functions
    void statusMonitor(const double currMinCost, const int currBestIndex);
    void cutStrategy();
    void incremental_Strategy(const int &maxGen,const int &currGen);
    
    int calDeltaSize(const int &maxGen,const int &currGen);
    double calDelta();
    
private:
    // APTDE variables and coefficient
    Status status;
    const int uPopSize;
    const int lPopSize;
    int k1;
    int k2;
    int p;
    std::vector<int> dimIndex;
};

// EPSDE
class EPSDE : public DifferentialEvolution {
public:
    EPSDE(const ObjectFunction &costFunction,
          const unsigned int &argPopulationSize,
          const bool &isCheckConstraints);
    
    void selectionAndCrossing(const int &maxGen,const int &currGen) override;
    
    // EPSDE member funciotns
    void initCoefficientAndMutation();
    
private:
    double probSuceess;
    std::vector<Strategy> strategies;
};

// MPEDE
class MPEDE : public DifferentialEvolution {
public:
    MPEDE(const ObjectFunction &costFunction,
          const unsigned int &argPopulationSize,
          const bool &isCheckConstraints);
    
    void selectionAndCrossing(const int &maxGen,const int &currGen) override;
    
    // MPEDE member functions
    void calCoef();
    int isStrategy(const int &i);
    
private:
    // MPEDE coefficients and variables
    const int ng; // The number of generation
    const double subPopRate;
    const double c; // The control parameter of mean
    int k; // The index of reward subpopulation
    std::vector<Strategy> strategies;
};


#endif /* DEalgo_hpp */
