#include "ps.h"

void ps::findPath(int from, int to, int **edge, vector<int> &path, int depNum, int n) {
    if (to < depNum + 1) {
        path.push_back(to);
        return;
    }
    path.push_back(to);
    int nextPosition;
    for (nextPosition = 1; nextPosition < n; nextPosition++) {
        if (edge[to][nextPosition] && nextPosition != from)break;
    }
    if (nextPosition == n)nextPosition = from;
    findPath(to, nextPosition, edge, path, depNum, n);
}

void ps::solve(int &n, int &cliNum, int &depNum, int &capacity,
               vector<int> &demand, vector<double> &openCost, double **dist) {
    vector<bool> openFlag(depNum);
    for (int i = 0; i < depNum; i++)openFlag[i] = false;
    factory::openStatus(depNum, cliNum, capacity, dist, openCost, openFlag,demand);
    // construct G'
    double distPrime[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            distPrime[i][j] = dist[i][j];
            distPrime[j][i] = dist[j][i];
            if (!constant::plusDelta)continue;
            if (i == j)continue;
            if (j > 0 && j < depNum + 1 && i > depNum && !openFlag[j - 1]) {
                distPrime[i][j] += 0.25 * openCost[j - 1];
                distPrime[j][i] += 0.25 * openCost[j - 1];
            }
        }
    }
    // construct H
    // get c(r,v) means that the nearest distance
    // between costumer v and open deps
    int nearestDepForC[cliNum];
    int tempDep = -1;
    double tempDist;
    for (int i = 0; i < cliNum; i++) {
        tempDist = MAX;
        for (int j = 0; j < depNum; j++) {
            if (distPrime[i + depNum + 1][j + 1] < tempDist) {
                tempDist = distPrime[i + depNum + 1][j + 1];
                tempDep = j;
            }
        }
        nearestDepForC[i] = tempDep;
    }
    // record edge update or not
    bool updateEdgeInH[n][n];
    auto **distInH = new double *[cliNum + 1];
    for (int i = 0; i < cliNum + 1; i++) {
        distInH[i] = new double[cliNum + 1];
    }
    for (int i = 0; i < cliNum + 1; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                distInH[i][j] = 0;
                continue;
            }
            if (j == 0) {
                tempDist = distPrime[i + depNum][nearestDepForC[i - 1] + 1];
                distInH[i][j] = tempDist;
                distInH[j][i] = tempDist;
                continue;
            }
            tempDist = min(distPrime[i + depNum][j + depNum],
                           distPrime[i + depNum][nearestDepForC[i - 1] + 1] +
                           distPrime[j + depNum][nearestDepForC[j - 1] + 1]);
            if (tempDist != distPrime[i + depNum][j + depNum]) {
                updateEdgeInH[i + depNum][j + depNum] = true;
                updateEdgeInH[j + depNum][i + depNum] = true;
            }
            distInH[i][j] = tempDist;
            distInH[j][i] = tempDist;
        }
    }
    // seek tsp
    vector<vector<int> > mst_info(cliNum + 2, vector<int>(cliNum + 2));
    static tsp tspSeek;
    tspSeek.seekTSP(cliNum + 1, mst_info, distInH);
    vector<int> tspAns;
    static tour tourSeek;
    double tspCost = tourSeek.seek(cliNum + 1, mst_info, distInH, tspAns);
    if (constant::openLog) cout << "TSP cost: " << tspCost << endl;
//    sort(tspAns.begin(),tspAns.end());
//    for (int t:tspAns)cout << t << " ";
//    cout << endl;
    // get paths and cycles from tsp answer
    int **edge = new int *[n];
    for (int i = 0; i < n; i++) {
        edge[i] = new int[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            edge[i][j] = 0;
        }
    }
    int from, to;
    for (int i = 0; i < tspAns.size(); i++) {
        from = tspAns[i];
        if (i == tspAns.size() - 1)to = tspAns[0];
        else to = tspAns[i + 1];
        if (from == 1) {
            edge[nearestDepForC[to - 2] + 1][to + depNum - 1]++;
            edge[to + depNum - 1][nearestDepForC[to - 2] + 1]++;
        } else {
            if (to == 1) {
                edge[nearestDepForC[from - 2] + 1][from + depNum - 1]++;
                edge[from + depNum - 1][nearestDepForC[from - 2] + 1]++;
            } else {
                if (updateEdgeInH[from + depNum - 1][to + depNum - 1]) {
                    edge[nearestDepForC[to - 2] + 1][to + depNum - 1]++;
                    edge[to + depNum - 1][nearestDepForC[to - 2] + 1]++;
                    edge[nearestDepForC[from - 2] + 1][from + depNum - 1]++;
                    edge[from + depNum - 1][nearestDepForC[from - 2] + 1]++;
                } else {
                    edge[from + depNum - 1][to + depNum - 1]++;
                    edge[to + depNum - 1][from + depNum - 1]++;
                }
            }
        }
    }
//    for(int i=0;i<n;i++)cout<<edge[21][i]<<" ";cout<<endl;
//    int i = 0;
    int sum;
    int *degree = new int[n];
    if (constant::openLog) cout << "degree != 2 node: ";
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < n; j++) {
            sum += edge[i][j];
        }
        degree[i] = sum;
        if (constant::openLog) if (sum != 2)cout << i << ":" << sum << " ";
    }
    if (constant::openLog) cout << endl;
    vector<vector<int>> paths;
    vector<int> path;
    for (int i = 1; i <= depNum; i++) {
        for (int j = 1 + depNum; j < n;) {
            if (edge[i][j]) {
                path.clear();
                path.push_back(i);
                findPath(i, j, edge, path, depNum, n);
                paths.push_back(path);
                // update edge
                for (int k = 0; k < path.size() - 1; k++) {
                    edge[path[k]][path[k + 1]]--;
                    edge[path[k + 1]][path[k]]--;
                }
            } else {
                j++;
            }
        }
    }
    // check edge == 0
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (edge[i][j])cout << "edge " << i << " " << j << " = " << edge[i][j] << endl;
        }
    }
//    for (auto p:paths)for (int i:p)cout << i << " ";
//    cout << endl;
//    i=0;
//    for (auto m:paths){
//        for (int g:m){
//            cout<<g<<" ";
//        }
//        cout<<"size: "<<m.size()<<endl;
//    }
    if (constant::shortCutInPS) {
        // shortcut deps
        vector<vector<int>> depPaths(depNum);
        for (int i = 0; i < paths.size(); i++) {
            auto p = paths[i];
            if (p[0] == p[p.size() - 1]) {
                depPaths[p[0] - 1].push_back(i);
            } else {
                depPaths[p[0] - 1].push_back(i);
                depPaths[p[p.size() - 1] - 1].push_back(i);
            }
        }
        for (int i = 0; i < depNum;) {
            if (depPaths[i].size() > 1) {
                bool containCycle = false;
                int p_1 = depPaths[i][0];
                int p_2 = depPaths[i][1];
                // p_2 is cycle if exists
                for (int j = 1; j < depPaths[i].size(); j++) {
                    int p_2_cycle = depPaths[i][j];
                    if (paths[p_2_cycle][0] == paths[p_2_cycle][paths[p_2_cycle].size() - 1]) {
                        p_2 = p_2_cycle;
                        containCycle = true;
//                    cout<<"find a cycle in p_2"<<endl;
                        break;
                    }
                }
                if(paths[p_1][0] == paths[p_1][paths[p_1].size() - 1])containCycle= true;
                if(constant::shortCutCycleInPS){
                    if(!containCycle){
//                        cout<<"path to path is skip."<<endl;
                        i++;
                        continue;
                    }
                }
                if (paths[p_1][0] == i + 1)reverse(paths[p_1].begin(), paths[p_1].end());
                if (paths[p_2][0] != i + 1)reverse(paths[p_2].begin(), paths[p_2].end());
                paths[p_1].pop_back();
                paths[p_2].erase(paths[p_2].begin());
                for (int p:paths[p_2]) {
                    paths[p_1].push_back(p);
                }
                bool p_1_isCycle = paths[p_1][0] == i + 1;
                bool p_2_isCycle = paths[p_2][paths[p_2].size() - 1] == i + 1;
                if (p_1_isCycle || p_2_isCycle) {
                    for (auto rmP = depPaths[i].begin(); rmP != depPaths[i].end(); rmP++) {
                        if (*rmP == p_2) {
                            depPaths[i].erase(rmP);
                            break;
                        }
                    }
                } else {
                    depPaths[i].erase(depPaths[i].begin());
                    depPaths[i].erase(depPaths[i].begin());
                }
                int p_2_to = paths[p_2][paths[p_2].size() - 1] - 1;
                for (auto it = depPaths[p_2_to].begin(); it != depPaths[p_2_to].end();) {
                    if (*it == p_2) {
                        if (paths[p_1][0] != paths[p_1][paths[p_1].size() - 1]) {
                            *it = p_1;
                            break;
                        } else {
                            depPaths[p_2_to].erase(it);
                            break;
                        }
                    } else {
                        it++;
                    }
                }
                paths[p_2].clear();
            } else {
                i++;
            }
        }
//    for (auto m:paths) {
//        for (int g:m) {
//            cout << g << " ";
//        }
//        cout << "size: " << m.size() << endl;
//    }
    }

    // reserve the nearest dep
    for (int i = 0; i < paths.size(); i++) {
        auto &p = paths[i];
        if (p.empty())continue;
        double openFrom = openFlag[p[0] - 1] ? 0 : openCost[p[0] - 1];
        double openTo = openFlag[p[p.size() - 1] - 1] ?
                        0 : openCost[p[p.size() - 1] - 1];
        double distFrom = dist[p[0]][p[1]];
        double distTo = dist[p[p.size() - 1]][p[p.size() - 2]];
        if (openFrom == openTo) {
            if (distFrom <= distTo) {
                p.pop_back();
            } else {
                p.erase(p.begin());
                reverse(p.begin(), p.end());
            }
        } else {
            if (openFrom < openTo) {
                p.pop_back();
            } else {
                p.erase(p.begin());
                reverse(p.begin(), p.end());
            }
        }
    }

//    for (auto m:paths){
//        for (int g:m){
//            cout<<g<<" ";
//        }
//        cout<<"size: "<<m.size()<<endl;
//    }

    // construct Tree
    double totalCost = 0.0;
    auto *tree = new treeNode{nullptr, vector<treeNode *>(), 0, 0, 0};
    treeNode *parent;
    for (auto p:paths) {
        parent = tree;
        while (!p.empty()) {
            treeNode *child;
            if (p[0] >= depNum + 1)
                child = new treeNode{parent, vector<treeNode *>(), p[0], demand[p[0] - depNum - 1], 0};
            else
                child = new treeNode{parent, vector<treeNode *>(), p[0], 0, 0};
            parent->child.push_back(child);
            parent = child;
            p.erase(p.begin());
        }
    }
    // cal total demand
    ts::calTotalDemand(tree);
    if (constant::openLog) cout << "all open status: ";
    for (auto p:paths)if (!p.empty())openFlag[p[0] - 1] = true;
    for (int i = 0; i < depNum; i++) {
        if (constant::openLog) cout << openFlag[i] << " ";
    }
    if (constant::openLog) cout << endl;
    for (int i = 0; i < depNum; i++)if (openFlag[i])totalCost += openCost[i];
    // remove closed factory
    for (auto childIt = tree->child.begin(); childIt != tree->child.end();) {
        auto child = *childIt;
        if (!openFlag[child->index - 1]) {
            cout << "remove a closed dep!" << endl;
            tree->child.erase(childIt);
            ts::freeTree(child);
        } else {
            childIt++;
        }
    }
    // cal total cost
    static tour t;
    vector<vector<int>> tourAns;
    int tourNum = 0;
    if (constant::openLog) cout << "totalDemand: " << tree->totalDemand << endl;
    while (!tree->child.empty()) {
        treeNode *factoryNode = tree->child[0];
        auto factoryIt = tree->child.begin();
        if (capacity >= factoryNode->totalDemand) {
            if (!factoryNode->child.empty())
                totalCost += t.seek(factoryNode, dist, true, tourAns, depNum);
            if (constant::openLog) cout << "tour " << ++tourNum << " demand = " << factoryNode->totalDemand << endl;
            tree->child.erase(factoryIt);
            ts::freeTree(factoryNode);
        } else {
            treeNode *miniTree = nullptr;
            vector<treeNode *> scanQueue;
            scanQueue.push_back(factoryNode);
            while (!scanQueue.empty()) {
                miniTree = scanQueue.at(scanQueue.size() - 1);
                scanQueue.pop_back();
                if (miniTree->totalDemand <= capacity)continue;
                bool isMinimal = true;
                for (auto child:miniTree->child) {
                    if (child->totalDemand > capacity)isMinimal = false;
                    scanQueue.push_back(child);
                }
                if (isMinimal) {
                    break;
                }
            }
            int demandDelta = miniTree->demand - (miniTree->totalDemand - capacity);
            vector<treeNode *> allNodes;
            queue<treeNode *> optQ;
            optQ.push(miniTree);
            while (!optQ.empty()) {
                treeNode *cur = optQ.front();
                optQ.pop();
                for (auto child:cur->child) {
                    optQ.push(child);
                }
                if (cur == miniTree && demandDelta == 0)continue;
                allNodes.push_back(cur);
            }
            vector<int> openDeps;
            for (int i = 0; i < depNum; i++) {
                if (openFlag[i])openDeps.push_back(i);
            }
            // find minimum cost edge
            int u, v;
            int minU = 0, minV = 0;
            double minEdgeCost = MAX;
            for (auto node:allNodes) {
                for (int i:openDeps) {
                    u = node->index;
                    v = i + 1;
                    if (dist[u][v] < minEdgeCost) {
                        minEdgeCost = dist[u][v];
                        minU = u;
                        minV = v;
                    }
                }
            }
            bool rootIsFactory = false;
            if (miniTree->index < depNum + 1) {
                cout << "minimal ST's root is factoryNode" << endl;
                rootIsFactory = true;
            }
            auto minVNode = new treeNode{nullptr, vector<treeNode *>(), minV, 0, 0};
            // create a new root node
            vector<treeNode *> newChild;
            for (auto c:allNodes) {
                if (c->parent == miniTree)newChild.push_back(c);
                if (c->index == minU && !rootIsFactory) {
                    if (c->index != miniTree->index) {
                        c->child.push_back(minVNode);
                    } else {
                        newChild.push_back(minVNode);
                    }
                    minVNode->parent = c;
                }
            }
            auto *newRootNode = new treeNode{nullptr, newChild, miniTree->index, miniTree->demand,
                                             miniTree->totalDemand};
            totalCost += t.seek(newRootNode, dist, demandDelta != 0 && !rootIsFactory, tourAns, depNum);
            for (auto childIt = miniTree->child.begin(); childIt != miniTree->child.end();) {
                auto child = *childIt;
                for(auto newRootChildIt=newRootNode->child.begin();newRootChildIt!=newRootNode->child.end();){
                    auto newRootChild=*newRootChildIt;
                    if(child->index==newRootChild->index){
                        newRootNode->child.erase(newRootChildIt);
                    } else{
                        newRootChildIt++;
                    }
                }
                miniTree->child.erase(childIt);
                ts::freeTree(child);
            }
            ts::freeTree(newRootNode);
            miniTree->demand -= demandDelta;
            treeNode *parents = miniTree;
            while (parents != nullptr) {
                parents->totalDemand -= capacity;
                parents = parents->parent;
            }
            if (constant::openLog) cout << "tour " << ++tourNum << " demand = " << capacity << endl;
        }
    }

    // check correct or not
    double totalDistance = 0.0;
    for (int i = 0; i < depNum; i++) {
        if (openFlag[i]) {
            totalDistance += openCost[i];
        }
    }
    set<int> allCustomers;
    int num = 0;
    for (auto tour:tourAns) {
        double distance = 0.0;
        for (int i = 0; i < tour.size(); i++) {
            from = tour[i];
            if (i == tour.size() - 1) to = tour[0];
            else to = tour[i + 1];
            allCustomers.insert(from);
            distance += dist[from][to];
        }
        totalDistance += distance;
        if (constant::openLog) cout << "tour " << ++num << " distance = " << distance << endl;
    }
    for (int i = depNum + 1; i < n; i++) {
        if (allCustomers.find(i) == allCustomers.end())
            cout << "customer " << i << " is not satisfied!" << endl;
    }
    if (constant::openLog) cout << "Total Distance: " << totalDistance << endl;
    if (totalDistance != totalCost)cout << "Total Cost != Total Distance!" << endl;
    cout << totalCost << "\t";
    if (constant::openLog) cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
}

