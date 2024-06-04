#ifndef TTRPCMAKE_GRAPH_H
#define TTRPCMAKE_GRAPH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "constant.h"

using namespace std;

class graph {
public:
    static double **buildGraphFromInstance(int &n, int &cliNum, int &depNum, int &capacity,
                                     vector<int> &demand, vector<double> &openCost);
};


#endif //TTRPCMAKE_GRAPH_H
