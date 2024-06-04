#ifndef TTRPCMAKE_TSP_H
#define TTRPCMAKE_TSP_H

#include <vector>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <climits>

class tsp {
private:
    struct edge {
        int v, u;
        double w;

        edge() {}

        edge(const int &_v, const int &_u, const double &_w)
                : v(_v), u(_u), w(_w) {}
    };

    static const int MaxN = 400;
    static const int MaxM = 79800;
    static const int MaxNX = MaxN + MaxN;

    int N, M;
    edge **mat;

    int n_matches;
    double tot_weight;
    double lab[MaxNX + 1];

    int q_n, q[MaxN];
    int fa[MaxNX + 1], col[MaxNX + 1];
    int slackv[MaxNX + 1];

    int n_x;
    int bel[MaxNX + 1], blofrom[MaxNX + 1][MaxN + 1];
    std::vector<int>* bloch;
    inline double e_delta(const edge &e);
    inline void update_slackv(int v, int x);
    inline void calc_slackv(int x);
    inline void q_push(int x);
    inline void set_mate(int xv, int xu);
    inline void set_bel(int x, int b);
    inline void augment(int xv, int xu);
    inline int get_lca(int xv, int xu);
    inline void add_blossom(int xv, int xa, int xu);
    inline void expand_blossom1(int b);
    inline void expand_blossom_final(int b);
    inline bool on_found_edge(const edge &e);
    bool match();
    void calc_max_weight_match();
    void MST(std::vector<std::vector<double> > &Matrix, int n, std::vector<std::vector<int> > &info);
    void dfs(int x, int num, std::vector<std::vector<int> > &GG);
public:
    tsp() {
        mat = new edge *[MaxNX + 1];
        for (int i = 0; i < MaxNX + 1; i++) {
            mat[i] = new edge[MaxNX + 1];
        }
        bloch=new std::vector<int>[MaxNX + 1];
    };
    int mate[MaxNX + 1];
    double MWPM(int numberOfOddNode, std::vector<std::vector<double> > &oddMatrix);
    void seekTSP(int n, std::vector<std::vector<int>>& mst_info, double** dist);
};


#endif //TTRPCMAKE_TSP_H
