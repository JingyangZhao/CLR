#include "graph.h"

using namespace std;

double **graph::buildGraphFromInstance(int &n, int &cliNum, int &depNum, int &capacity,
                                       vector<int> &demand, vector<double> &openCost) {
    ifstream instanceStream;
    instanceStream.open(constant::instanceFileName, ios::in);
    if (!instanceStream.is_open()) {
        cout << "open fail." << endl;
        exit(1);
    }
    double temp;
    vector<double> tempVector;
    vector<vector<double>> cli, dep;
    int tourCost;
    instanceStream >> cliNum >> depNum >> capacity >> tourCost;
    for (int i = 0; i < 4; i++)instanceStream >> temp;
    for (int i = 0; i < cliNum; i++) {
        for (int j = 0; j < 4; j++) {
            instanceStream >> temp;
            tempVector.push_back(temp);
        }
        cli.push_back(tempVector);
        tempVector.clear();
    }
    for (int i = 0; i < depNum; i++) {
        for (int j = 0; j < 6; j++) {
            instanceStream >> temp;
            tempVector.push_back(temp);
        }
        dep.push_back(tempVector);
        tempVector.clear();
    }
    n = cliNum + depNum + 1;
    double **completeG;
    completeG = new double *[n];
    for (int i = 0; i < n; i++) {
        completeG[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                completeG[i][j] = 0.0;
                continue;
            }
            if (j == 0) {
                if (i < depNum + 1) {
                    completeG[i][j] = 0.0;
                    completeG[j][i] = 0.0;
                    continue;
                } else {
                    completeG[i][j] = MAX;
                    completeG[j][i] = MAX;
                    continue;
                }
            }
            double uX, uY, vX, vY;
            if (i < depNum + 1) {
                uX = dep.at(i - 1).at(1);
                uY = dep.at(i - 1).at(2);
            } else {
                uX = cli.at(i - depNum - 1).at(1);
                uY = cli.at(i - depNum - 1).at(2);
            }
            if (j < depNum + 1) {
                vX = dep.at(j - 1).at(1);
                vY = dep.at(j - 1).at(2);
            } else {
                vX = cli.at(j - depNum - 1).at(1);
                vY = cli.at(j - depNum - 1).at(2);
            }
            double d = sqrt(pow(uX - vX, 2) + pow(uY - vY, 2));
            if ( j < depNum + 1 && i > depNum)d += 0.5 * tourCost;
            completeG[i][j] = d;
            completeG[j][i] = d;
        }
    }

    for (int i = 0; i < cliNum; i++) {
        demand.push_back(cli[i][3]);
    }
    for (int i = 0; i < depNum; i++) {
        openCost.push_back(dep[i][3]);
    }

//    for(int i=0;i<n;i++){
//        for(int j=0;j<n;j++){
//            cout<<completeG[i][j]<<" ";
//        }
//        cout<<endl;
//    }

    instanceStream.close();
    return completeG;
}