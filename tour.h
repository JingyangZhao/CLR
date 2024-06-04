#ifndef TTRPCMAKE_TOUR_H
#define TTRPCMAKE_TOUR_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <climits>
#include <map>
#include "constant.h"

using namespace std;

struct treeNode {
    treeNode *parent;
    vector<treeNode *> child;
    int index;
    int demand;
    int totalDemand;
};
class tour {
    struct stack {
        int top, node[M];
    } sEuler;
    void dfs(int x, int num, vector<vector<int> > &GG);
    void fleury(int num, vector<vector<int> > &info, vector<int> &euler);
    double
    cal_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, int delta, vector<int> &ans, bool containRoot,
            int rootToIndex, vector<int>& demand, map<int, int> &indexToTree, int depNum);
    double
    seek_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, map<int, int> &indexToTree, bool containRoot,
             int rootToIndex, vector<vector<int>> &tourAns,vector<int>& demand, int depNum);
    double seek_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, vector<int> &tourAns);
public:
    double seek(treeNode* rootNode,double** weight,bool containRoot,vector<vector<int>>& tourAns,int depNum);
    double seek(int num, vector<vector<int> >& mst_info, double** weight,vector<int>& tourAns);
};


#endif //TTRPCMAKE_TOUR_H
