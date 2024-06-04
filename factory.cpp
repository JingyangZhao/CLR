#include "factory.h"
#include "constant.h"

using namespace std;

bool factory::myfunction(star_struct &a, star_struct &b) {
    return a.starcost < b.starcost;
}

void factory::openStatus(int m, int n, int capacity, double **connectionCosts, vector<double> &openCost,
                         vector<bool> &openFlag, vector<int> &demand) {

    //m: number of depots   n:number of customers
    //double connectionCosts[n+m+1][n+m+1]; // matrix
    double openCosts[m + 1];
    //input
    //   connectionCosts[n+m+1][n+m+1]   openCosts[i]: depot i

    vector<int> V;
    vector<bool> open_flag(m + 1);

    // receive data from parameters
    for (int i = 1; i < m + 1; i++) {
        openCosts[i] = openCost[i - 1];
        open_flag[i] = openFlag[i - 1];
    }

    for (int i = 1; i <= n; i++) {
        V.push_back(i + m);
    }

    while (V.size() > 0) {

        double ratio = 99999999.9;
        vector<int> Vplus;
        int openWhich = 0;

        for (int i = 1; i <= m; i++) {
            vector<star_struct> Star(V.size());

            for (int j = 0; j < V.size(); j++) {
                Star[j].starcost = connectionCosts[i][V[j]] * demand[V[j] - 1] * (2.0 / capacity);
                Star[j].index = V[j];
            }

            sort(Star.begin(), Star.end(), myfunction);

            double fi, tem = 0.0;
            int temd = 0;

            if (open_flag[i] == false)
                fi = openCosts[i];
            else
                fi = 0.0;
            for (int k = 0; k < Star.size(); k++) {
                tem += Star[k].starcost;
                temd += demand[Star[k].index-1];

                double r = ((constant::greedyAlpha/100.0) * fi + tem) / temd;
                if (r < ratio) {
//                    printf("%lf\n", r);
                    ratio = r;
                    Vplus.clear();
                    for (int kk = 0; kk <= k; kk++) {
                        Vplus.push_back(Star[kk].index);
                    }
                    openWhich = i;
                }
            }

            vector<star_struct>().swap(Star);
        }

        if (open_flag[openWhich] == false) open_flag[openWhich] = true;

        vector<int> ans;

        set_difference(V.begin(), V.end(), Vplus.begin(), Vplus.end(), back_inserter(ans));

        V.clear();

        for (int i = 0; i < ans.size(); i++) V.push_back(ans[i]);

        vector<int>().swap(ans);
        vector<int>().swap(Vplus);
    }

//    printf("open Status: \n");
    if (constant::openLog) cout << "greedy open status: ";
    for (int i = 1; i <= m; i++) {
        openFlag[i - 1] = open_flag[i];
        if (constant::openLog) cout << open_flag[i] << " ";
    }
    if (constant::openLog) cout << endl;

    vector<int>().swap(V);
    vector<bool>().swap(open_flag);
}
