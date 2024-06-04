#ifndef TTRPCMAKE_TS_H
#define TTRPCMAKE_TS_H

#include <cstring>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
#include <algorithm>
#include "constant.h"
#include "tour.h"
#include "factory.h"

using namespace std;

class ts {
    int cnte;
    int tempDemand;
    int *h;
    double *dis;
    bool *vis;
    double res;
    int cnt;
    struct E {
        double w;
        int v, x;
    } *e;
    treeNode *tree;
    struct S {
        int u;
        double d;
        treeNode *tPtr;
    };
    priority_queue<S> q;

    friend bool operator<(const S &x, const S &y) { return x.d > y.d; }

    void adde(int u, int v, double w);

    void prim(int &n, int &depNum, vector<int> &demand);

public:
    ts() {
        cnte = 0, tempDemand = 0;
        h = new int[N];
        dis = new double[N];
        vis = new bool[N];
        for (int i = 0; i < N; i++) {
            h[i] = 0;
            dis[i] = MAX;
            vis[i] = false;
        }
        res = 0.0;
        cnt = 0;
        tree = nullptr;
        priority_queue<S>().swap(q);
        e = new E[M * 2];
    }
    ~ts(){
        free(h);
        free(dis);
        free(vis);
        free(tree);
        priority_queue<S>().swap(q);
        free(e);
    }

    static void freeTree(treeNode *root);

    static int calTotalDemand(treeNode *treeNode);

    void solve(int &n, int &cliNum, int &depNum, int &capacity,
               vector<int>& demand, vector<double>& openCost, double **dist);
};


#endif //TTRPCMAKE_TS_H
