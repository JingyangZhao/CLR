#ifndef TTRPCMAKE_PS_H
#define TTRPCMAKE_PS_H

#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include <vector>
#include "tsp.h"
#include "ts.h"
#include "tour.h"
#include "constant.h"
#include "factory.h"

using namespace std;

class ps {
    static void findPath(int from, int to, int **edge, vector<int> &path, int depNum, int n);
public:
    void solve(int& n,int& cliNum,int&depNum,int&capacity,
               vector<int>& demand,vector<double>& openCost, double** dist);
};


#endif //TTRPCMAKE_PS_H
