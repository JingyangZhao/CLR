#include "tour.h"
#include "tsp.h"

using namespace std;

void tour::dfs(int x, int num, vector<vector<int> > &GG) {
    int i;
    sEuler.node[++sEuler.top] = x;

    for (i = 0; i < num; i++) {
        if (GG[i][x] > 0) {
            GG[i][x]--;
            GG[x][i]--; //deleting this edge
            dfs(i, num, GG);
            break;
        }
    }
}

//Fleury algorithm for Euler Circuit
void tour::fleury(int num, vector<vector<int> > &info, vector<int> &euler) {
    vector<vector<int> > GG(num, vector<int>(num));
    for (int i = 0; i < num; i++)
        for (int j = 0; j < num; j++)
            GG[i][j] = info[i + 1][j + 1];

    int x = 0;
    int i, flag;
    sEuler.top = 0;
    sEuler.node[sEuler.top] = x;
    while (sEuler.top >= 0) {
        flag = 0;
        for (i = 0; i < num; i++) {
            if (GG[sEuler.node[sEuler.top]][i] > 0) {
                flag = 1;
                break;
            }
        }
        if (!flag) {
            euler.push_back(sEuler.node[sEuler.top--] + 1);
            //printf("%d ",s.node[s.top--]+1);
        } else dfs(sEuler.node[sEuler.top--], num, GG);
    }

    for (int i = 0; i < num; i++) vector<int>().swap(GG[i]);
    vector<vector<int> >().swap(GG);
}

//calculate the total distance of tsp
double
tour::cal_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, int delta, vector<int> &ans, bool containRoot,
        int rootToIndex, vector<int>& demand, map<int, int> &indexToTree, int depNum) {
    int len = euler.size();
    vector<bool> flagVisit(num + 1);
    // shortcut root node or not
    if (!containRoot)flagVisit[rootToIndex] = true;
    // shortcut demand==0 && not dep node
    for(int i=1;i<demand.size();i++){
        if(indexToTree[i]<=depNum)continue;
        if(demand[i]==0)flagVisit[i]=true;
    }
    for (int i = 0; i < len; i++) {
        int position = (i + delta) % (euler.size());
        if (flagVisit[euler[position]] == false) {
            flagVisit[euler[position]] = true;
            ans.push_back(euler[position]);
        }
    }

    double sumCost = 0;
    for (int i = 0; i < ans.size() - 1; i++) {
        double cost = matrix[ans[i]][ans[i + 1]];
        sumCost += cost;
    }
    int a = ans[ans.size() - 1], b = ans[0];
    sumCost += matrix[a][b];

    vector<bool>().swap(flagVisit);

    return sumCost;
}

//After the MWM in the 3/2-appro TSP, we can start with n different vertex
//thus seeking for a TSP with a smaller total weight
double
tour::seek_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, map<int, int> &indexToTree, bool containRoot,
         int rootToIndex, vector<vector<int>> &tourAns,vector<int>& demand, int depNum) {
    double minCostTour = MAX;
    vector<int> minAns;
    vector<int> curAns;

    for (int i = 0; i < euler.size(); i++) {
        curAns.clear();
        double res = cal_tsp(num, matrix, euler, i, curAns, containRoot, rootToIndex,demand,indexToTree,depNum);

        if (res < minCostTour) {
            minCostTour = res;
            minAns.clear();
            for (int curAn : curAns) {
                minAns.push_back(curAn);
            }
        }
    }

    // check and reserve
    vector<int> tour;
    bool containDep = false;
    for (int minAn : minAns) {
        if(constant::openLog) cout<<indexToTree[minAn]<<" ";
        tour.push_back(indexToTree[minAn]);
        if(indexToTree[minAn]<=depNum)containDep=true;
    }
    if(!containDep)cout<<"tour don't contain dep!"<<endl;
    if(constant::openLog) cout<<endl;
    tourAns.push_back(tour);
    if(constant::openLog) cout<<minCostTour<<endl;
    return minCostTour;
}

double tour::seek_tsp(int num, vector<vector<double> > &matrix, vector<int> &euler, vector<int> &tourAns) {
    double minCostTour = INT_MAX;
    vector<int> minAns;
    vector<int> curAns;
    vector<int> demand;
    map<int,int> indexToTree;

    for (int i = 0; i < euler.size(); i++) {
        curAns.clear();
        double res = cal_tsp(num, matrix, euler, i, curAns, true, 1,demand,indexToTree,0);

        if (res < minCostTour) {
            minCostTour = res;
            minAns.clear();
            for (int curAn : curAns) {
                minAns.push_back(curAn);
            }
        }
    }
    tourAns.clear();
    for (int minAn : minAns) {
        tourAns.push_back(minAn);
    }

    return minCostTour;
}

double tour::seek(treeNode *rootNode, double **weight, bool containRoot, vector<vector<int>> &tourAns, int depNum) {
    vector<int> euler;

    int num = 0;
    map<int, int> treeToIndex;
    map<int, int> indexToTree;
    vector<int> demand;demand.push_back(0);
    vector<treeNode *> queue;
    queue.push_back(rootNode);
    while (!queue.empty()) {
        treeNode *curNode = queue.at(queue.size() - 1);
        queue.pop_back();
        for (auto child:curNode->child) {
            queue.push_back(child);
        }
        demand.push_back(curNode->demand);
        num++;
        treeToIndex[curNode->index] = num;
        indexToTree[num] = curNode->index;
    }

    vector<vector<int> > mst_info(num + 1, vector<int>(num + 1)); //for edge numbers
    vector<vector<double> > matrix(num + 1, vector<double>(num + 1)); // for edge weights

    queue.push_back(rootNode);
    while (!queue.empty()) {
        treeNode *curNode = queue.at(queue.size() - 1);
        queue.pop_back();
        for (auto child:curNode->child) {
            queue.push_back(child);
        }
        if (curNode == rootNode)continue;
        int to = curNode->index, from = curNode->parent->index;
        mst_info[treeToIndex[to]][treeToIndex[from]] = 1;
        mst_info[treeToIndex[from]][treeToIndex[to]] = 1;
    }

    //A simple test for MST(if there is a node with degree=0, then error)
    vector<int> odd_node;
    for (int i = 1; i <= num; i++) {
        int d = 0;
        for (int j = 1; j <= num; j++) {
            d += mst_info[i][j];
        }
        if (d == 0) printf("error!::::node degree = 0\n");
        else if (d % 2 == 1) odd_node.push_back(i);
    }

    for (int i = 1; i <= num; i++) {
        for (int j = 1; j <= num; j++) {
            matrix[i][j] = weight[indexToTree[i]][indexToTree[j]];
        }
    }

    int numberOfOddNode = odd_node.size();
    vector<vector<double> > oddMatrix(numberOfOddNode + 1, vector<double>(numberOfOddNode + 1));
    for (int i = 1; i <= numberOfOddNode; i++)
        for (int j = 1; j <= numberOfOddNode; j++)
            oddMatrix[i][j] = matrix[odd_node[i - 1]][odd_node[j - 1]];
    tsp t;
    t.MWPM(numberOfOddNode,oddMatrix);
    //Find a Eulerian Graph
    for (int i = 1; i <= numberOfOddNode; i++) mst_info[odd_node[i - 1]][odd_node[t.mate[i] - 1]]++;

    fleury(num, mst_info, euler);
    euler.pop_back();
//    for(int i=0; i<euler.size(); i++) printf("%d ", indexToTree[euler[i]]); printf("\n");

    double res = seek_tsp(num, matrix, euler, indexToTree, containRoot, 1, tourAns,demand,depNum);

    for (int i = 0; i <= num; i++) vector<int>().swap(mst_info[i]), vector<double>().swap(matrix[i]);
    vector<vector<int> >().swap(mst_info);
    vector<vector<double> >().swap(matrix);
    vector<int>().swap(euler);

    return res;
}

double tour::seek(int num, vector<vector<int> > &mst_info, double **weight, vector<int> &tourAns) {
    vector<int> euler;

//    vector<vector<int> >mst_info(num+1, vector<int>(num+1)); //for edge numbers
    vector<vector<double> > matrix(num + 1, vector<double>(num + 1)); // for edge weights

    for (int i = 1; i <= num; i++) {
        for (int j = 1; j <= num; j++) {
            matrix[i][j] = weight[i - 1][j - 1];
        }
    }

    fleury(num, mst_info, euler);
    euler.pop_back();
//    for(int i=0; i<euler.size(); i++) printf("%d ", indexToTree[euler[i]]); printf("\n");

    double res = seek_tsp(num, matrix, euler, tourAns);

    for (int i = 0; i <= num; i++) vector<int>().swap(mst_info[i]), vector<double>().swap(matrix[i]);
    vector<vector<int> >().swap(mst_info);
    vector<vector<double> >().swap(matrix);
    vector<int>().swap(euler);

    return res;
}
