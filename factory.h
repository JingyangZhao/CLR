#ifndef TTRPCMAKE_FACTORY_H
#define TTRPCMAKE_FACTORY_H

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

class factory {
    struct star_struct
    {
        double starcost;
        int index;
    };
    static bool myfunction(star_struct& a, star_struct& b);
public:
    static void openStatus(int m, int n, int capacity, double** connectionCosts, vector<double>& openCost,vector<bool>& openFlag,vector<int>&demand);
};


#endif //TTRPCMAKE_FACTORY_H
