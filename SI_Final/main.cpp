//
//  main.cpp
//  SI_Final
//
//  Created by Aaron on 2021/6/9.
//

#include <iostream>
#include "testFunction.hpp"
#include "DEalgo.hpp"
using namespace std;

void test(const int DIMS = 5,
          const int POPULATIONSIZE = 50,
          const int GENERATION = 20,
          int detail = 0);

int main() {
    const int dims = 30;
    const int populationSize = 50;
    const int generation = 10000;
    // 0:population, 1:行數+min, 2:行數+eva+mim, 3:eva+min, 4:eva+min(csv), 5: min
    const int detail = 3;
    const int times = 1;
    
    for (int i=0; i<times; ++i) {
        test(dims, populationSize, generation, detail);
    }
    
}

void test(const int DIMS, const int POPULATIONSIZE, const int GENERATION,
          int detail){
    
    vector<string> funcName = {
        "Ackley", "Griewank", "Michalewicz", "Powell", "Rastrigin",
        "Rosenbrock", "Schwefel", "Sphere", "Sum Squares", "Zakharov"};
    
    vector<ObjectFunction*> cost(funcName.size());
    cost[0] = new Ackley      (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    cost[1] = new Griewank    (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    cost[2] = new Michalewicz (DIMS); // d=5 -> f(x*) = -4.687658
    cost[3] = new Powell      (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    cost[4] = new Rastrigin   (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    
    cost[5] = new Rosenbrock  (DIMS); // f(x*)=0, at x* = (1, 1, ..., 1)
    cost[6] = new Schwefel    (DIMS); // f(x*)=0, at x* = (420.9687,., 420.9687)
    cost[7] = new Sphere      (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    cost[8] = new SumSquares  (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    cost[9] = new Zakharov    (DIMS); // f(x*)=0, at x* = (0, 0, ...., 0)
    
    // 4 different DE
    vector<string> deName = {"Basic DE", "APTDE", "EPSDE", "MPEDE"};
    vector<DifferentialEvolution> Basic_DE;
    vector<APTDE> APT_DE;
    vector<EPSDE> EPS_DE;
    vector<MPEDE> MPE_DE;
    for (int i=0; i<cost.size(); ++i) {
        Basic_DE.push_back(DifferentialEvolution(*cost[i], POPULATIONSIZE, true));
        APT_DE.push_back(APTDE(*cost[i], POPULATIONSIZE, true));
        EPS_DE.push_back(EPSDE(*cost[i], POPULATIONSIZE, true));
        MPE_DE.push_back(MPEDE(*cost[i], POPULATIONSIZE, true));
    }
    
    // Optimize
    for (int i=0; i<10; ++i) {
        cout << deName[0] << endl;
        cout << funcName[i] << endl;
        Basic_DE[i].optimize(GENERATION, detail);

        cout << deName[1] << endl;
        cout << funcName[i] << endl;
        APT_DE[i].optimize(GENERATION, detail);

        cout << deName[2] << endl;
        cout << funcName[i] << endl;
        EPS_DE[i].optimize(GENERATION, detail);

        cout << deName[3] << endl;
        cout << funcName[i] << endl;
        MPE_DE[i].optimize(GENERATION, detail);
    }
    
    // Free pointer
    for (auto &c : cost) {
        free(c);
        c = nullptr;
    }
    
}
