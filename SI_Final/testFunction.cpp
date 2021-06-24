//
//  testFunction.cpp
//  SI_Final
//
//  Created by Aaron on 2021/6/9.
//

#include "testFunction.hpp"
using namespace std;

/*  Constraints  */

// Default constructor
Constraints::Constraints(const double lower, const double upper, const bool isConstrained)
:lower(lower), upper(upper), isConstrained(isConstrained)
{}

// Check if a value is valid
bool Constraints::check(const double candidate){
    if (isConstrained){
        if (candidate <= upper && candidate >= lower){
            return true;
        }else{
            return false;
        }
    }else{
        return true;
    }
}
