//
//  DEalgo.cpp
//  SI_Final
//
//  Created by Aaron on 2021/6/9.
//

#include "DEalgo.hpp"
#include <iostream>
#include <vector>
#include <limits>
#include <ctime>
#include <random>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace std;

/* Basic DE (DE/rand/nin/1) */

// Defult constructor
DifferentialEvolution::DifferentialEvolution(
    const ObjectFunction &costFunction,
    const unsigned int &argPopulationSize,
    const bool &isCheckConstraints)
:cost(costFunction),
 populationSize(argPopulationSize),
 isCheckConstraints(isCheckConstraints),
 bestIndexPointer(0),
 minCost(numeric_limits<double>::infinity()),
 DIMENSION(costFunction.getDimension()),
 constraints(costFunction.getConstraints()),
 popIndex(populationSize)
{
    generator.seed(static_cast<unsigned int>(time(nullptr)));
    
    // The size of population is populationSize * DIMENSION
    for (int i=0; i<populationSize; ++i) {
        population.push_back(Point(DIMENSION));
    }
    
};

// Initialization of population
void DifferentialEvolution::initPopulation(){
    
    // Initialize population based on a uniform probability distribution of cost function
    uniform_real_distribution<double> *distribution;
    
    for ( auto &agent : population) {
        for (int i=0; i<DIMENSION; ++i) {
            if (constraints[i].getIsConstrained()) {
                distribution = new uniform_real_distribution<double>(constraints[i].getLower(), constraints[i].getUpper());
            }else{
                distribution = new uniform_real_distribution<double>(defaultLower, defaultUpper);
            }
            agent[i] = (*distribution)(generator);
            free(distribution);
            distribution = NULL;
        }
    }
    
    // Initialize min cost per point and best pointer
    for (int i=0; i<populationSize; ++i){
        population[i].minCost = cost.evaluteCost(population[i].vec);
        
        if (population[i].minCost<minCost) {
            minCost = population[i].minCost;
            bestIndexPointer = i;
        }
    }
    
    // Initialize popIndex
    for (int i=0; i<populationSize; ++i) {
        popIndex[i] = i;
    }
}

// Mutation -> crossing -> selection
void DifferentialEvolution::selectionAndCrossing(const int &maxGen,const int &currGen){
    uniform_real_distribution<double> distribution(0, populationSize);
    
    double currMinCost = population[0].minCost;
    int currBestIndexPointer = 0;
    
    for (int i=0; i<populationSize; ++i) {
        
        // Randomly select 3 different indexes are different from i pointer
        int randNum = 3;
        shuffle(popIndex.begin(), popIndex.end(), generator);
        vector<int> r(popIndex.begin(), popIndex.begin()+randNum);
        
        // r1, r2, r3 must be different from each other and from i
        for (int j=0; j<randNum; ++j) {
            if (r[j] == i) {r[j] = popIndex[j+randNum];}
        }
        
        // 1. Mutation: Form mutant solution of i pointer
        vector<double> mutant(DIMENSION);
        for (int j=0; j<DIMENSION; ++j) {
            mutant[j] = population[r[0]][j] + population[i].factor * (population[r[1]][j] - population[r[2]][j]);
        }
        
        // 2. Crossover
        vector<double> trial(DIMENSION, i);
        crossover(trial, mutant, i);
        
        // Check if trial satisfies constraints and skip it if not.
        if ( isCheckConstraints && !checkConstraints(trial) ) {
            --i;
            continue;
        }
        
        // 3. Selection: Calculate new cost and decide trial should be kept.
        double newCost = cost.evaluteCost(trial);
        if (newCost < population[i].minCost) {
            population[i].vec = trial;
            population[i].minCost = newCost;
        }
        
        // Track the global best point
        if (population[i].minCost < currMinCost) {
            currMinCost = population[i].minCost;
            currBestIndexPointer = i;
        }
    }
    
    // Update min cost and best index
    minCost = currMinCost;
    bestIndexPointer = currBestIndexPointer;
    
}

// The binomial crossover
void DifferentialEvolution::crossover(vector<double> &trial, const vector<double> &mutant, const int &i){
    // Randomly choose an index (R)
    shuffle(popIndex.begin(), popIndex.end(), generator);
    int R = popIndex[0];
    
    // Randomly choose probability (pr) for each dimension
    vector<double> pr(DIMENSION);
    uniform_real_distribution<double> probDistribution(0, 1);
    for (auto &var : pr){
        var = probDistribution(generator);
    }
    
    // Crossing
    for (int j=0; j<DIMENSION; ++j) {
        if (pr[j] < population[i].crossingRate || pr[j] == R) {
            trial[j] = mutant[j];
        }else{
            trial[j] = population[i][j];
        }
    }
    
    // Avoid the element of trial out-of-bounds
    for (auto &t : trial){
        if (t < cost.getLower()) {
            t = cost.getLower();
        }else if (t > cost.getUpper()){
            t = cost.getUpper();
        }
    }
}

// Print population
void DifferentialEvolution::printPopulation(){
    for ( const auto& pointer : population) {
        for ( const auto& var : pointer.vec) {
            cout << var << " ";
        }
        cout << endl;
    }
}

// Execute differential evolution
void DifferentialEvolution::optimize(int generation, int detail){
    int evaluation = 0;
    initPopulation();
    for (int i=0; i<generation; ++i) {
        selectionAndCrossing(generation, i);
        evaluation += populationSize;
        switch (detail) {
            case 0:
                cout << "## Generation " << i+1 << " :" << endl;
                cout << "Current min cost: " << minCost << endl;
                cout << "Best pointer:";
                for (int j=0; j<DIMENSION; ++j) {
                    cout << population[bestIndexPointer][j] << " ";
                }
                cout << endl;
                break;
                
            case 1:
                cout << i+1 << " " ;
                cout << minCost << endl;
                break;
            
            case 2:
                cout << i+1 << " ";
                cout << evaluation << " ";
                cout << minCost << endl;
                break;
                
            case 3:
                cout << evaluation << " ";
                cout << minCost << endl;
                break;
                
            case 4:
                cout << evaluation << " , ";
                cout << minCost << " , "<< endl;
                break;
        }
    }
    
    cout << defaultfloat << minCost << " " << endl;
    
    if (!detail) {
        cout << "Terminated." << endl;
    }
};

// private membert function
// Check if a point each dimension is valid
bool DifferentialEvolution::checkConstraints(std::vector<double> point){
    for (int i=0; i<DIMENSION; ++i) {
        if ( !constraints[i].check(point[i]) ) {
            return false;
        }
    }
    return true;
}

/* APTDE */

APTDE::APTDE(const ObjectFunction &costFunction,
             const unsigned int &argPopulationSize,
             const bool &isCheckConstraints)
:DifferentialEvolution(costFunction, argPopulationSize, isCheckConstraints),
 uPopSize(1.5 * argPopulationSize),
 lPopSize(0.5 * argPopulationSize),
 k1(4),
 k2(4),
 p(3),
 dimIndex(DIMENSION)
 {
     // Initialize dimIndex
     for (int i=0; i<DIMENSION; ++i) {
         dimIndex[i] = i;
     }
 }

void APTDE::selectionAndCrossing(const int &maxGen,const int &currGen){
    uniform_real_distribution<double> distribution(0, populationSize);
    
    double currMinCost = population[0].minCost;
    int currBestIndexPointer = 0;
    
    for (int i=0; i<populationSize; ++i) {
        
        // Randomly select 3 different indexes are different from i pointer
        int randNum = 3;
        shuffle(popIndex.begin(), popIndex.end(), generator);
        vector<int> r(popIndex.begin(), popIndex.begin()+randNum);
        
        // r1, r2, r3 must be different from each other and from i
        for (int j=0; j<randNum; ++j) {
            if (r[j] == i) {r[j] = popIndex[j+randNum];}
        }
        
        // 1. Mutation: Form mutant solution of i pointer
        vector<double> mutant(DIMENSION);
        for (int j=0; j<DIMENSION; ++j) {
            mutant[j] = population[r[0]][j] + population[i].factor * (population[r[1]][j] - population[r[2]][j]);
        }
        
        // 2. Crossover
        vector<double> trial(DIMENSION);
        crossover(trial, mutant, i);
        
        // Check if trial satisfies constraints and skip it if not.
        if ( isCheckConstraints && !checkConstraints(trial) ) {
            --i;
            continue;
        }
        
        // 3. Selection: Calculate new cost and decide trial should be kept.
        double newCost = cost.evaluteCost(trial);
        if (newCost < population[i].minCost) {
            population[i].vec = trial;
            population[i].minCost = newCost;
        }
        
        // Track the global best point
        if (population[i].minCost < currMinCost) {
            currMinCost = population[i].minCost;
            currBestIndexPointer = i;
        }
    }
    
    // Update status monitor each generation
    statusMonitor(currMinCost, currBestIndexPointer);
    
    // Decrease population size
    if ( (status.RM > k1 && populationSize >= lPopSize) || status.UM > k2 ) {
        cutStrategy();
        status.RM = 0;
        status.UM = 0;
    }
    
    // Increase population size
    if ( (status.NM > k1 && populationSize <= uPopSize) || status.LM > k2) {
        incremental_Strategy(maxGen, currGen);
        status.NM = 0;
        status.LM = 0;
    }
    
}

// Upadta status monitors, min cost and best index
void APTDE::statusMonitor(const double currMinCost, const int currBestIndex){
    // Update redundance (RM) and stagnation Monitor (NM)
    if (currMinCost < minCost) {
        ++ status.RM;
        status.NM = 0;
        
        // Update min cost and best index
        minCost = currMinCost;
        bestIndexPointer = currBestIndex;
    }else{
        status.RM = 0;
        ++ status.NM;
    }
    
    // Update upper (UM) and lower (LM) bound monitor
    if (populationSize >= uPopSize) {
        ++ status.UM;
        status.LM = 0;
    }else if (populationSize <= lPopSize){
        status.UM = 0;
        ++ status.LM;
    }
    
}


//
void APTDE::cutStrategy(){
    // Rearrange the population with its fitness function values from small to large
    sort(population.begin(), population.end(), cmp);
    sort(popIndex.begin(), popIndex.end(), cmp2);
    
    // 1. Calculate rank and select inferior points
    int gridSize = 0.8 * populationSize;
    double currMinCost = population[0].minCost;
    double currMaxCost = population[populationSize - 1].minCost;
    vector<int> IA; // Inferior-Archive
    for (int i=0; i<populationSize; ++i) {
        // Calculate rank each pointer
        if (i < populationSize - 1 ) {
            population[i].rank = (population[i].minCost - currMinCost) * gridSize / (currMaxCost - currMinCost);
        }else{
            population[i].rank = population[populationSize - 1 - 2].rank;
        }
        
        // Probabiltively select inferior points
        uniform_real_distribution<double> probDistribution(0, 1);
        double prob = probDistribution(generator);
        double t = 1 - 1 / (population[i].rank + 1); // Normalization
        if (prob < t) {
            IA.push_back(i);
        }
    }
    
    // 2. Delete the chosen individuals from population
    int delSize = std::min<int>(populationSize * p / 100, IA.size() );
    if (delSize == IA.size() ) {
        // Delete the chosen individuals from population
        auto popIter = population.begin();
        for (auto iter=IA.rbegin(); iter!=IA.rend(); ++iter) {
            population.erase(popIter + *iter);
            popIndex.pop_back(); // Pop popIndex
        }
    }else{
        // Randomly select Del worst individuals from IA
        shuffle(IA.begin(), IA.end(), generator);
        while (IA.size() != delSize) {IA.pop_back();}
        sort(IA.begin(), IA.end(), cmp2);
        
        // Delete the chosen individuals from population
        auto popIter = population.begin();
        for (auto iter=IA.rbegin(); iter!=IA.rend(); ++iter) {
            population.erase(popIter + *iter);
            popIndex.pop_back(); // Pop popIndex
        }
    }
    
    // Update populationSize
    populationSize = static_cast<int>(population.size());
    
    // Clear
    IA.clear();
}

//
void APTDE::incremental_Strategy(const int &maxGen,const int &currGen){
    // Rearrange the population with its fitness function values from small to large
    sort(population.begin(), population.end(), cmp);
    
    int incrSize = ceil(populationSize * p / 100); // The number of increment populations
    int deltaSize = calDeltaSize(maxGen, currGen); // The number of mutation dimensions
    for (int i=0; i<incrSize; ++i) {
        // Avoid deltaSize is greater than DIMENSION
        if (deltaSize > DIMENSION) {deltaSize = DIMENSION;}
        
        // Mutate deltaSize dimensions each elite copy
        shuffle(dimIndex.begin(), dimIndex.end(), generator);
        population.push_back(population[i]); // Copy elite into end
        for (int j=0; j<deltaSize; ++j) {
            double delta = calDelta();
            if (cost.getUpper() >= population[i][dimIndex[j]] + delta && population[i][dimIndex[j]] + delta >= cost.getLower() ) {
                population[population.size() - 1][dimIndex[j]] = population[i][dimIndex[j]] + delta;
            }else{
                --j;
            }
            
            // Calculate the cost value of mutated elite copy
            population[population.size() - 1].minCost = cost.evaluteCost(population[population.size() - 1].vec);
        }
        
        // Increase popIndex
        popIndex.push_back(populationSize + i - 1);
    }
    
    // Update populationSize
    populationSize = static_cast<int>(population.size());
}

// Calculate the size of delta
int APTDE::calDeltaSize(const int &maxGen,const int &currGen){
    const int UDelta = 5;
    const int LDelta = 1;
    const double a = UDelta - LDelta;
    const double b = LDelta / a;
    const double lamda = 2 / maxGen;
    return a * ( sqrt(1 - lamda * (currGen - 1 / lamda) ) / (sqrt(1 - lamda * (currGen - 1 / lamda) ) + sqrt(1 + lamda * (currGen - 1 / lamda) ) ) + b);
}

// Calculate the value of delta
double APTDE::calDelta(){
    normal_distribution<double> NDstr(0, 1/9); // mean value 0 and standard deviation 1/9.
    double R = abs(NDstr(generator));
    const double c1 = 0.5; // Seggest c1 = 0.02 ~ 0.7
    const double c2 = 0.01; // Seggest c2 = 0 ~ 0.02
    const double Ud = c1 * cost.getUpper(); // Upper boundaries of delta
    const double Ld = c2 * cost.getLower(); // Lower boundaries of delta
    
    while (R >= 1) {R = abs(NDstr(generator));} // Ensure the range of R is [0, 1)
    
    if (0 < R && R <= 0.5) {
        return (Ud - Ld) * (2 * sqrt(R) + Ld / (Ud - Ld) );
    }else if (0.5 <= R && R < 1){
        return (Ud - Ld) * (1 - 2 * sqrt(R - 1) + Ld / (Ud - Ld) );
    }
    return 0;
}

/* EPSDE */

EPSDE::EPSDE(const ObjectFunction &costFunction,
             const unsigned int &argPopulationSize,
             const bool &isCheckConstraints)
:DifferentialEvolution(costFunction, argPopulationSize, isCheckConstraints),
 probSuceess(0.9)
{
    // The size of strategies is populationSize
    for (int i=0; i<populationSize; ++i) {
        strategies.push_back(Strategy());
    }
    
    initCoefficientAndMutation();
}

// Initialization of coefficient (factors and crossing rates) and mutation strategies.
void EPSDE::initCoefficientAndMutation(){
    // The range of factor value is 0.4 ~ 0.9 in steps of 0.1
    uniform_real_distribution<double> fDistribution(0.4, 1);
    // The range of crossing rate value is 0.1 ~ 0.9 in steps of 0.1
    uniform_real_distribution<double> crDistribution(0.1, 1);

    // There are three different mutation strategie.
    uniform_int_distribution<int> mDistribution(0, 2);
    
    // Initialize factor, crossing rates and mutation strategie based on a uniform probability distribution
    for (int i=0; i<populationSize; ++i) {
        strategies[i].improve = 1;
        strategies[i].mutation = mDistribution(generator);
        population[i].factor = floor(fDistribution(generator) * 10) / 10;
        population[i].crossingRate = floor(crDistribution(generator) * 10) / 10;
    }
}

// Mutation -> crossing -> selection
void EPSDE::selectionAndCrossing(const int &maxGen,const int &currGen){
    uniform_int_distribution<int> distribution(0, populationSize - 1);
    
    double currMinCost = population[0].minCost;
    int currBestIndexPointer = 0;
    
    for (int i=0; i<populationSize; ++i) {
        
        // Randomly select 4 different indexes are different from i pointer
        int randNum = 4;
        shuffle(popIndex.begin(), popIndex.end(), generator);
        vector<int> r(popIndex.begin(), popIndex.begin()+randNum);
        
        // r1, r2, r3m r4 must be different from each other and from i
        for (int j=0; j<randNum; ++j) {
            if (r[j] == i) {r[j] = popIndex[j+randNum];}
        }
        
        // 1. Mutation: Form mutant solution of i pointer based on different mutation strategies
        vector<double> mutant(DIMENSION);
        for (int j=0; j<DIMENSION; ++j) {
            switch (strategies[i].mutation) {
                case 0: // DE/rand/1/bin
                    mutant[j] = population[r[0]][j] +
                    population[i].factor * (population[r[1]][j] - population[r[2]][j]);
                    break;
                    
                case 1: // DE/best/2/bin
                    mutant[j] = population[bestIndexPointer][j] +
                    population[i].factor * (population[r[0]][j] - population[r[1]][j]) +
                    population[i].factor * (population[r[2]][j] - population[r[3]][j]);
                    break;
                    
                case 2: // DE/current-to-rand/1/bin
                    uniform_real_distribution<double> kDistribution(0, 1);
                    double k = kDistribution(generator);
                    mutant[j] = population[i][j] +
                    k * (population[r[0]][j] - population[i][j]) +
                    population[i].factor * (population[r[1]][j] - population[r[2]][j]);
                    break;
            }
        }
        
        // 2. Crossover
        vector<double> trial(DIMENSION, i);
        crossover(trial, mutant, i);
        
        // Check if trial satisfies constraints and skip it if not.
        if ( isCheckConstraints && !checkConstraints(trial) ) {
            --i;
            continue;
        }
        
        // 3. Selection: Calculate new cost and decide trial should be kept.
        double newCost = cost.evaluteCost(trial);
        if (newCost < population[i].minCost) {
            population[i].vec = trial;
            population[i].minCost = newCost;
            strategies[i].improve = 1;
        }else{
            // Upadating factors, crossing rates and mutation strategies.
            strategies[i].improve = 0;
            uniform_real_distribution<double> probDistribution(0, 1);
            double pr = probDistribution(generator);
            if ( pr < probSuceess && probSuceess >= 2 / populationSize ) { // Reinitialization based on successful combination
                shuffle(popIndex.begin(), popIndex.end(), generator);
                int counter = 0;
                int pool = popIndex[counter];
                while (pool == i && strategies[i].improve == 0 ) {pool = popIndex[++counter];}
                
                population[i].factor = population[pool].factor;
                population[i].crossingRate = population[pool].crossingRate;
                strategies[i].mutation = strategies[pool].mutation;
            }else{ // Random reinitialization
                // The range of factor value is 0.4 ~ 0.9 in steps of 0.1
                static uniform_real_distribution<double> fDistribution(0.4, 1);
                // The range of crossing rate value is 0.1 ~ 0.9 in steps of 0.1
                static uniform_real_distribution<double> crDistribution(0.1, 1);
                // There are three different mutation strategie.
                static uniform_int_distribution<int> mDistribution(0, 2);
                
                population[i].factor = floor(fDistribution(generator) * 10) / 10;
                population[i].crossingRate = floor(crDistribution(generator) * 10) / 10;
                strategies[i].mutation = mDistribution(generator);
            }
        }
        
        // Updating the probability of successful combination
        int sum = 0;
        for (int j=0; j<populationSize; ++j) {
            if (strategies[j].improve == 1) {++sum;}
        }
        probSuceess = sum / populationSize;
        
        // Track the global best point
        if (population[i].minCost < currMinCost) {
            currMinCost = population[i].minCost;
            currBestIndexPointer = i;
        }
        
    }
    
    // Update min cost and best index
    minCost = currMinCost;
    bestIndexPointer = currBestIndexPointer;
    
}

/* MPEDE */

MPEDE::MPEDE(const ObjectFunction &costFunction,
             const unsigned int &argPopulationSize,
             const bool &isCheckConstraints)
:DifferentialEvolution(costFunction, argPopulationSize, isCheckConstraints),
 ng(20),
 subPopRate(0.2),
 c(0.5)
{
    // There are 3 diifferent Strategies
    for (int i=0; i<3; ++i) {
        strategies.push_back(Strategy(0.5));
    }
    
    // Randomly initialize the index (k) of reward subpopulation
    uniform_int_distribution<int> kDst(0,2);
    k = kDst(generator);
}

void MPEDE::selectionAndCrossing(const int &maxGen,const int &currGen){
    uniform_real_distribution<double> distribution(0, populationSize);
    
    double currMinCost = population[0].minCost;
    int currBestIndexPointer = 0;
    
    calCoef(); // Calculate CR and F
    
    // Mutation -> crossing -> selection each Point
    for (int i=0; i<populationSize; ++i) {
        
        // Randomly select 4 different indexes are different from i pointer
        int randNum = 4;
        shuffle(popIndex.begin(), popIndex.end(), generator);
        vector<int> r(popIndex.begin(), popIndex.begin()+randNum);
        
        // r1, r2, r3, r4 must be different from each other and from i
        for (int j=0; j<randNum; ++j) {
            if (r[j] == i) {r[j] = popIndex[j+randNum];}
        }
                
        // Assign a strategy each pointer
        const int strI = isStrategy(i);
        
        // 1. Mutation: Form mutant solution of i pointer
        vector<double> mutant(DIMENSION);
        switch (strI) {
            case 0: // DE/rand/1/bin
                for (int j=0; j<DIMENSION; ++j) {
                    mutant[j] = population[r[0]][j] +
                    population[i].factor * (population[r[1]][j] - population[r[2]][j]);
                }
                break;
            
            case 1: // DE/current-to-best/1/bin
                for (int j=0; j<DIMENSION; ++j) {
                    mutant[j] = population[i][j] +
                    population[i].factor* (population[bestIndexPointer][j] - population[i][j]) +
                    population[i].factor * (population[r[0]][j] - population[r[1]][j]);
                }
                break;
                
            case 2: // DE/current-to-rand/1/bin
                uniform_real_distribution<double> kDistribution(0, 1);
                double k;
                
                for (int j=0; j<DIMENSION; ++j) {
                    k = kDistribution(generator);
                    mutant[j] = population[i][j] +
                    k * (population[r[0]][j] - population[i][j]) +
                    population[i].factor * (population[r[1]][j] - population[r[2]][j]);
                }
                break;
        }
        
        // 2. Crossover
        vector<double> trial(DIMENSION, i);
        crossover(trial, mutant, i);
        
        // Check if trial satisfies constraints and skip it if not.
        if ( isCheckConstraints && !checkConstraints(trial) ) {
            --i;
            continue;
        }
        
        // 3. Selection: Calculate new cost and decide trial should be kept.
        double newCost = cost.evaluteCost(trial);
        if (newCost < population[i].minCost) {
            // Record successful improvment, fes,  F and CR
            strategies[strI].improve += population[i].minCost - newCost;
            ++strategies[strI].fes;
            strategies[strI].succF.push_back(population[i].factor);
            strategies[strI].succCR.push_back(population[i].crossingRate);
            
            // Update global best
            population[i].vec = trial;
            population[i].minCost = newCost;
        }
        
        // Track the global best point
        if (population[i].minCost < currMinCost) {
            currMinCost = population[i].minCost;
            currBestIndexPointer = i;
        }
    }
    
    // Update min cost and best index
    minCost = currMinCost;
    bestIndexPointer = currBestIndexPointer;
    
    // Reward the highest performance of a strategy every ng generation
    if (currGen % ng == 0) {
        double max = 0;
        int maxI = 0;
        for (int i=0; i<3; ++i) {
            if (max < strategies[i].improve / strategies[i].fes) {
                max = strategies[i].improve / strategies[i].fes;
                maxI = i;
            }
            // Reset improvment and fes
            strategies[i].improve = 0;
            strategies[i].fes = 0;
        }
        // Upadte k
        k = maxI;
    }
    
    // Shuffle populations
    shuffle(population.begin(), population.end(), generator);
}

// Calculate F and CR
void MPEDE::calCoef(){
    for (int i=0; i<3; ++i) {
        // Update the mean of factor
        if (strategies[i].succF.size() > 0) {
            double sum1 = 0;
            double sum2 = 0;
            for (auto &succ : strategies[i].succF ) {
                sum1 += succ * succ;
                sum2 += succ;
            }
            strategies[i].succF.clear();
            // Lehmer mean
            strategies[i].meanF = (1 - c) * strategies[i].meanF + c * sum1 / sum2;
        }
        
        // Update the mean of crossing rate
        if (strategies[i].succCR.size() >0) {
            double sum1 = 0;
            double sum2 = 0;
            for (auto &succ : strategies[i].succCR ) {
                sum1 += succ;
                ++sum2;
            }
            strategies[i].succCR.clear();
            strategies[i].meanCR = (1 - c) * strategies[i].meanCR + c * sum1 / sum2;
        }
    }
    
    // Update factor and crossing rate
    for (int i=0; i<DIMENSION; ++i) {
        cauchy_distribution<double> FDst(strategies[isStrategy(i)].meanF, 0.1);
        normal_distribution<double> CRDst(strategies[isStrategy(i)].meanCR, 0.1);
        population[i].factor = FDst(generator);
        population[i].crossingRate = CRDst(generator);
        
        if (population[i].factor < 0.4 ) {
            population[i].factor = 0.4;
        }else if (population[i].factor > 0.9){
            population[i].factor = 0.9;
        }

        if (population[i].crossingRate < 0.1 ) {
            population[i].crossingRate = 0.1;
        }else if (population[i].crossingRate > 0.9){
            population[i].crossingRate = 0.9;
        }
        
    }
    
}

// Determine a strategy into pointer
int MPEDE::isStrategy(const int &i){
    // Detetmine subpopulation and reward subpopulation size
    int subPopSize = ceil(subPopRate * populationSize);
    int rewardPopSize = populationSize - 3 * subPopSize;
    
    if ( (0 <= i && i <= subPopSize - 1) ||
         ( 3 * subPopSize <= i && i <= rewardPopSize - 1 && k == 0 ) ) {
        return 0;
    }else if ((subPopSize <= i && i <= 2 * subPopSize - 1) ||
              ( 3 * subPopSize <= i && i <= rewardPopSize - 1 && k == 1 )) {
        return 1;
    }else if ((2 * subPopSize <= i && i <= 3 * subPopSize - 1) ||
              ( 3 * subPopSize <= i && i <= rewardPopSize - 1 && k == 2 )){
        return 2;
    }
    
    return 0; // error
}
