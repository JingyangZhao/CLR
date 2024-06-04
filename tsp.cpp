#include "tsp.h"

using namespace std;

const double INF = 2147483647;

template<class T>
inline void tension(T &a, const T &b) {
    if (b < a)
        a = b;
}

template<class T>
inline void relax(T &a, const T &b) {
    if (b > a)
        a = b;
}
//template <class T>
//inline int size(const T &a)
//{
//    return (int)a.size();
//}

//inline int getint() {
//    char c;
//    while (c = getchar(), '0' > c || c > '9');
//
//    int res = c - '0';
//    while (c = getchar(), '0' <= c && c <= '9')
//        res = res * 10 + c - '0';
//    return res;
//}


inline double tsp::e_delta(const edge &e) {
    return lab[e.v] + lab[e.u] - mat[e.v][e.u].w * 2;
}

inline void tsp::update_slackv(int v, int x) {
    if (!slackv[x] || e_delta(mat[v][x]) < e_delta(mat[slackv[x]][x]))
        slackv[x] = v;
}

inline void tsp::calc_slackv(int x) {
    slackv[x] = 0;
    for (int v = 1; v <= N; v++)
        if (mat[v][x].w > 0 && bel[v] != x && col[bel[v]] == 0)
            update_slackv(v, x);
}

inline void tsp::q_push(int x) {
    if (x <= N)
        q[q_n++] = x;
    else {
        for (int i = 0; i < size(bloch[x]); i++)
            q_push(bloch[x][i]);
    }
}

inline void tsp::set_mate(int xv, int xu) {
    mate[xv] = mat[xv][xu].u;
    if (xv > N) {
        edge e = mat[xv][xu];
        int xr = blofrom[xv][e.v];
        int pr = find(bloch[xv].begin(), bloch[xv].end(), xr) - bloch[xv].begin();
        if (pr % 2 == 1) {
            reverse(bloch[xv].begin() + 1, bloch[xv].end());
            pr = size(bloch[xv]) - pr;
        }

        for (int i = 0; i < pr; i++)
            set_mate(bloch[xv][i], bloch[xv][i ^ 1]);
        set_mate(xr, xu);

        rotate(bloch[xv].begin(), bloch[xv].begin() + pr, bloch[xv].end());
    }
}

inline void tsp::set_bel(int x, int b) {
    bel[x] = b;
    if (x > N) {
        for (int i = 0; i < size(bloch[x]); i++)
            set_bel(bloch[x][i], b);
    }
}

inline void tsp::augment(int xv, int xu) {
    while (true) {
        int xnu = bel[mate[xv]];
        set_mate(xv, xu);
        if (!xnu)
            return;
        set_mate(xnu, bel[fa[xnu]]);
        xv = bel[fa[xnu]], xu = xnu;
    }
}

inline int tsp::get_lca(int xv, int xu) {
    static bool book[MaxNX + 1];
    for (int x = 1; x <= n_x; x++)
        book[x] = false;
    while (xv || xu) {
        if (xv) {
            if (book[xv])
                return xv;
            book[xv] = true;
            xv = bel[mate[xv]];
            if (xv)
                xv = bel[fa[xv]];
        }
        swap(xv, xu);
    }
    return 0;
}

inline void tsp::add_blossom(int xv, int xa, int xu) {
    int b = N + 1;
    while (b <= n_x && bel[b])
        b++;
    if (b > n_x)
        n_x++;

    lab[b] = 0;
    col[b] = 0;

    mate[b] = mate[xa];

    bloch[b].clear();
    bloch[b].push_back(xa);
    for (int x = xv; x != xa; x = bel[fa[bel[mate[x]]]])
        bloch[b].push_back(x), bloch[b].push_back(bel[mate[x]]), q_push(bel[mate[x]]);
    reverse(bloch[b].begin() + 1, bloch[b].end());
    for (int x = xu; x != xa; x = bel[fa[bel[mate[x]]]])
        bloch[b].push_back(x), bloch[b].push_back(bel[mate[x]]), q_push(bel[mate[x]]);

    set_bel(b, b);

    for (int x = 1; x <= n_x; x++) {
        mat[b][x].w = mat[x][b].w = 0;
        blofrom[b][x] = 0;
    }
    for (int i = 0; i < size(bloch[b]); i++) {
        int xs = bloch[b][i];
        for (int x = 1; x <= n_x; x++)
            if (mat[b][x].w == 0 || e_delta(mat[xs][x]) < e_delta(mat[b][x]))
                mat[b][x] = mat[xs][x], mat[x][b] = mat[x][xs];
        for (int x = 1; x <= n_x; x++)
            if (blofrom[xs][x])
                blofrom[b][x] = xs;
    }
    calc_slackv(b);
}

inline void tsp::expand_blossom1(int b) // lab[b] == 1
{
    for (int i = 0; i < size(bloch[b]); i++)
        set_bel(bloch[b][i], bloch[b][i]);

    int xr = blofrom[b][mat[b][fa[b]].v];
    int pr = find(bloch[b].begin(), bloch[b].end(), xr) - bloch[b].begin();
    if (pr % 2 == 1) {
        reverse(bloch[b].begin() + 1, bloch[b].end());
        pr = size(bloch[b]) - pr;
    }

    for (int i = 0; i < pr; i += 2) {
        int xs = bloch[b][i], xns = bloch[b][i + 1];
        fa[xs] = mat[xns][xs].v;
        col[xs] = 1, col[xns] = 0;
        slackv[xs] = 0, calc_slackv(xns);
        q_push(xns);
    }
    col[xr] = 1;
    fa[xr] = fa[b];
    for (int i = pr + 1; i < size(bloch[b]); i++) {
        int xs = bloch[b][i];
        col[xs] = -1;
        calc_slackv(xs);
    }

    bel[b] = 0;
}

inline void tsp::expand_blossom_final(int b) // at the final stage
{
    for (int i = 0; i < size(bloch[b]); i++) {
        if (bloch[b][i] > N && lab[bloch[b][i]] == 0)
            expand_blossom_final(bloch[b][i]);
        else
            set_bel(bloch[b][i], bloch[b][i]);
    }
    bel[b] = 0;
}

inline bool tsp::on_found_edge(const edge &e) {
    int xv = bel[e.v], xu = bel[e.u];
    if (col[xu] == -1) {
        int nv = bel[mate[xu]];
        fa[xu] = e.v;
        col[xu] = 1, col[nv] = 0;
        slackv[xu] = slackv[nv] = 0;
        q_push(nv);
    } else if (col[xu] == 0) {
        int xa = get_lca(xv, xu);
        if (!xa) {
            augment(xv, xu), augment(xu, xv);
            for (int b = N + 1; b <= n_x; b++)
                if (bel[b] == b && lab[b] == 0)
                    expand_blossom_final(b);
            return true;
        } else
            add_blossom(xv, xa, xu);
    }
    return false;
}

bool tsp::match() {
    for (int x = 1; x <= n_x; x++)
        col[x] = -1, slackv[x] = 0;

    q_n = 0;
    for (int x = 1; x <= n_x; x++)
        if (bel[x] == x && !mate[x])
            fa[x] = 0, col[x] = 0, slackv[x] = 0, q_push(x);
    if (q_n == 0)
        return false;

    while (true) {
        for (int i = 0; i < q_n; i++) {
            int v = q[i];
            for (int u = 1; u <= N; u++)
                if (mat[v][u].w > 0 && bel[v] != bel[u]) {
                    double d = e_delta(mat[v][u]);
                    if (d < 1e-6 && d > -1e-6) {
                        if (on_found_edge(mat[v][u]))
                            return true;
                    } else if (col[bel[u]] == -1 || col[bel[u]] == 0)
                        update_slackv(v, bel[u]);
                }
        }

        double d = INF;
        for (int v = 1; v <= N; v++)
            if (col[bel[v]] == 0)
                tension(d, lab[v]);
        for (int b = N + 1; b <= n_x; b++)
            if (bel[b] == b && col[b] == 1)
                tension(d, lab[b] / 2);
        for (int x = 1; x <= n_x; x++)
            if (bel[x] == x && slackv[x]) {
                if (col[x] == -1)
                    tension(d, e_delta(mat[slackv[x]][x]));
                else if (col[x] == 0)
                    tension(d, e_delta(mat[slackv[x]][x]) / 2);
            }

        for (int v = 1; v <= N; v++) {
            if (col[bel[v]] == 0)
                lab[v] -= d;
            else if (col[bel[v]] == 1)
                lab[v] += d;
        }
        for (int b = N + 1; b <= n_x; b++)
            if (bel[b] == b) {
                if (col[bel[b]] == 0)
                    lab[b] += d * 2;
                else if (col[bel[b]] == 1)
                    lab[b] -= d * 2;
            }

        q_n = 0;
        for (int v = 1; v <= N; v++)
            if (lab[v] < 1e-6 && lab[v] > -1e-6) // all unmatched vertices' labels are zero! cheers!
                return false;
        for (int x = 1; x <= n_x; x++)
            if (bel[x] == x && slackv[x] && bel[slackv[x]] != x && e_delta(mat[slackv[x]][x]) < 1e-6 && e_delta(mat[slackv[x]][x]) > -1e-6 ) {
                if (on_found_edge(mat[slackv[x]][x]))
                    return true;
            }
        for (int b = N + 1; b <= n_x; b++)
            if (bel[b] == b && col[b] == 1 && lab[b] < 1e-6 && lab[b] > -1e-6)
                expand_blossom1(b);
    }
    return false;
}

void tsp::calc_max_weight_match() {
    for (int v = 1; v <= N; v++)
        mate[v] = 0;

    n_x = N;
    n_matches = 0;
    tot_weight = 0;

    bel[0] = 0;
    for (int v = 1; v <= N; v++)
        bel[v] = v, bloch[v].clear();
    for (int v = 1; v <= N; v++)
        for (int u = 1; u <= N; u++)
            blofrom[v][u] = v == u ? v : 0;

    double w_max = 0;
    for (int v = 1; v <= N; v++)
        for (int u = 1; u <= N; u++)
            relax(w_max, mat[v][u].w);
    for (int v = 1; v <= N; v++)
        lab[v] = w_max;

    while (match())
        n_matches++;

    for (int v = 1; v <= N; v++)
        if (mate[v] && mate[v] < v)
            tot_weight += mat[v][mate[v]].w;
}

double tsp::MWPM(int numberOfOddNode, vector<vector<double> > &oddMatrix) {
    N = numberOfOddNode, M = N * (N - 1) / 2;
    for (int v = 1; v <= N; v++)
        for (int u = 1; u <= N; u++)
            mat[v][u] = edge(v, u, 0);

    for (int i = 1; i <= N; i++)
        for (int j = i + 1; j <= N; j++) {
            double w = oddMatrix[i][j];
            //converting MWPM to maximum weight perfect matching problem
            //by adding a big number to each edge's weight
            w = -w + 1000000;
            mat[i][j].w = mat[j][i].w = w;
        }
    //maximum weight matching probelm
    calc_max_weight_match();

    return (double)N / 2 * 1000000 - tot_weight;
}

void tsp::MST(vector<vector<double> > &Matrix, int n, vector<vector<int> > &info) {
    //Using prime algorithm, start with vertex 1
    vector<bool> flag(n + 1);//recording whether vertex choosen or not
    for (int i = 1; i <= n; i++) flag[i] = false;
    vector<double > L(n + 1);
    for (int i = 1; i <= n; i++) L[i] = Matrix[1][i];
    vector<int> chosen;
    flag[1] = true;//start with vertex 1
    chosen.push_back(1);

    //searching minimum weight edge
    for (int i = 1; i <= n - 1; i++) {
        int now = 0;
        double min1 = INT_MAX;
        for (int j = 1; j <= n; j++) {
            if (flag[j] == false && L[j] < min1) {
                now = j;
                min1 = L[j];
            }
        }
        if (now == 0) {
            printf("mst generating error!\n");
            break;//in case that graph in not connected
        }

        flag[now] = true;
        chosen.push_back(now);
        for (int j = 1; j <= chosen.size(); j++)
            if (min1 == Matrix[now][chosen[j - 1]]) {
                info[now][chosen[j - 1]] = 1;
                info[chosen[j - 1]][now] = 1;
                break;
            }
        //updating vertices' minimum weight distance from source
        for (int j = 1; j <= n; j++) {
            if (flag[j] == false && L[j] > Matrix[now][j])
                L[j] = Matrix[now][j];
        }
    }
}

struct stack {
    int top, node[1000000];
} s;

void tsp::dfs(int x, int num, vector<vector<int> > &GG) {
    int i;
    s.node[++s.top] = x;

    for (i = 0; i < num; i++) {
        if (GG[i][x] > 0) {
            GG[i][x]--;
            GG[x][i]--; //deleting this edge
            dfs(i, num, GG);
            break;
        }
    }
}

void tsp::seekTSP(int n, vector<vector<int>> &mst_info, double **dist) {
//    int n = 10;
//    vector<vector<int> > mst_info(n + 1, vector<int>(n + 1));
    vector<vector<double> > Matrix(n + 1, vector<double>(n + 1));

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            Matrix[i][j] = dist[i - 1][j - 1];
        }
    }

    MST(Matrix, n, mst_info);

    //A simple test for MST(if there is a node with degree=0, then error)
    vector<int> odd_node;
    for (int i = 1; i <= n; i++) {
        int d = 0;
        for (int j = 1; j <= n; j++) {
            d += mst_info[i][j];
        }
        if (d == 0) printf("error!::::node degree = 0\n");
        else if (d % 2 == 1) odd_node.push_back(i);
    }

    int numberOfOddNode = odd_node.size();
    vector<vector<double> > oddMatrix(numberOfOddNode + 1, vector<double>(numberOfOddNode + 1));
    for (int i = 1; i <= numberOfOddNode; i++)
        for (int j = 1; j <= numberOfOddNode; j++)
            oddMatrix[i][j] = Matrix[odd_node[i - 1]][odd_node[j - 1]];

    //Seeking for a minimum perfect matching in graph HH
    MWPM(numberOfOddNode, oddMatrix);

    //Find a Eulerian Graph
    for (int i = 1; i <= numberOfOddNode; i++) mst_info[odd_node[i - 1]][odd_node[mate[i] - 1]]++;

}
