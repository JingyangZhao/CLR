#include "ts.h"

using namespace std;

void ts::adde(int u, int v, double w) {
    e[++cnte] = E{w, v, h[u]}, h[u] = cnte;
}

void ts::prim(int &n, int &depNum, vector<int> &demand) {
    dis[0] = 0;
    auto *rootChild = new vector<treeNode *>;
    auto *rootNode = new treeNode{nullptr, *rootChild, 0, 0, 0};
    tree = rootNode;
    q.push({0, 0, rootNode});
    treeNode *parent;
    treeNode *curNode;
    vector<treeNode *> *child;
    treeNode *node;
    while (!q.empty()) {
        if (cnt >= n) break;
        int u = q.top().u;
        double d = q.top().d;
        parent = q.top().tPtr;
        q.pop();
        if (vis[u]) continue;
        vis[u] = true;
        ++cnt;
        res += d;
        child = new vector<treeNode *>;
        if (u >= 1 + depNum) {
            tempDemand = demand[u - depNum - 1];
        } else {
            tempDemand = 0;
        }
        node = new treeNode{parent, *child, u, tempDemand, 0};
        parent->child.push_back(node);
        curNode = node;
        for (int i = h[u]; i; i = e[i].x) {
            int v = e[i].v;
            double w = e[i].w;
            if (w < dis[v]) {
                dis[v] = w, q.push({v, w, curNode});
            }
        }
    }
    tree = tree->child.at(0);
}

int ts::calTotalDemand(treeNode *treeNode) {
    int totalDemand = 0;
    if (treeNode->child.empty()) {
        totalDemand = treeNode->demand;
    } else {
        totalDemand += treeNode->demand;
        for (auto tPtr:treeNode->child) {
            totalDemand += calTotalDemand(tPtr);
        }
    }
    treeNode->totalDemand = totalDemand;
    return totalDemand;
}

void ts::freeTree(treeNode *root) {
    if (root->child.empty()) {
        free(root);
        return;
    }
    for (auto c:root->child) {
        freeTree(c);
    }
}

void ts::solve(int &n, int &cliNum, int &depNum, int &capacity,
               vector<int> &demand, vector<double> &openCost, double **dist) {
    vector<bool> openFlag(depNum);
    for (int i = 0; i < depNum; i++)openFlag[i] = false;
    factory::openStatus(depNum, cliNum, capacity, dist, openCost, openFlag, demand);

    //int temopen = 0;
    //for (int i = 0; i < depNum; i++) if(openFlag[i]) temopen ++; cout<< temopen<<" ";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (constant::plusDelta) {
                if (j > 0 && j < depNum + 1 && i > depNum && openFlag[j - 1] == false) {
                    adde(i, j, dist[i][j] + 0.5 * openCost[j - 1]);
                    adde(j, i, dist[j][i] + 0.5 * openCost[j - 1]);
                } else {
                    adde(i, j, dist[i][j]);
                    adde(j, i, dist[j][i]);
                }
            } else {
                adde(i, j, dist[i][j]);
                adde(j, i, dist[j][i]);
            }
        }
    }
    prim(n, depNum, demand);
    calTotalDemand(tree);
    if (cnt == n) {
        if (constant::openLog) cout << "MST cost: " << res << endl;
    } else
        cout << "No MST." << endl;
    // CLR
    double totalCost = 0.0;
    if (constant::openLog) cout << "all open status: ";
    for (int i = 0; i < depNum; i++) {
        int j = 0;
        for (j = 0; j < tree->child.size(); j++) {
            if (i + 1 == tree->child[j]->index)break;
        }
        if (!openFlag[i])openFlag[i] = !tree->child[j]->child.empty();
        if (constant::openLog) cout << openFlag[i] << " ";
    }
    if (constant::openLog) cout << endl;
    for (int i = 0; i < depNum; i++)if (openFlag[i])totalCost += openCost[i];
    // remove closed factory
    for (auto childIt = tree->child.begin(); childIt != tree->child.end();) {
        auto child = *childIt;
        if (!openFlag[child->index - 1]) {
            tree->child.erase(childIt);
            freeTree(child);
        } else {
            childIt++;
        }
    }
    static tour t;
    vector<vector<int>> tourAns;
    while (!tree->child.empty()) {
        treeNode *factoryNode = tree->child[0];
        auto factoryIt = tree->child.begin();
        // containRoot = true means that cal root in shortcut
        bool containRoot;
        if (capacity >= factoryNode->totalDemand) {
            containRoot = true;
            if (!factoryNode->child.empty())
                totalCost += t.seek(factoryNode, dist, containRoot, tourAns, depNum);
            tree->totalDemand -= factoryNode->totalDemand;
            tree->child.erase(factoryIt);
            freeTree(factoryNode);
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
            // greedy partition the miniTree
            vector<vector<treeNode *>> splitTree;
            vector<treeNode *> leftNodes;
            for (auto node:miniTree->child)leftNodes.push_back(node);
            vector<treeNode *> optST;
            int optScore = 0;
//            if (miniTree->demand <= capacity) {
//                optST.push_back(miniTree);
//                optScore += miniTree->demand;
//            }
            bool containMiniTree = false;
            while (!leftNodes.empty() || !containMiniTree) {
                for (auto it = leftNodes.begin(); it != leftNodes.end();) {
                    auto node = *it;
                    if (optScore + node->totalDemand <= capacity) {
                        optST.push_back(node);
                        optScore += node->totalDemand;
                        leftNodes.erase(it);
                    } else {
                        it++;
                    }
                }
                if (!containMiniTree && optScore + miniTree->demand <= capacity) {
                    optST.push_back(miniTree);
                    containMiniTree = true;
                }
                splitTree.push_back(optST);
                optST.clear();
                optScore = 0;
            }
            for (vector<treeNode *> nodes:splitTree) {
                int score = 0;
                containRoot = false;
                vector<treeNode *> allNodes;
                queue<treeNode *> tempQueue;
                for (auto node:nodes) {
                    if (node == miniTree) {
                        score += node->demand;
                        containRoot = true;
                        allNodes.push_back(node);
                    } else {
                        score += node->totalDemand;
                        tempQueue.push(node);
                        while (!tempQueue.empty()) {
                            treeNode *curNode = tempQueue.front();
                            tempQueue.pop();
                            for (auto child:curNode->child) {
                                tempQueue.push(child);
                            }
                            allNodes.push_back(curNode);
                        }
                    }
                }
                if (score <= capacity / 2)continue;
                vector<int> openDeps;
                for (int i = 0; i < depNum; i++) {
                    if (openFlag[i])openDeps.push_back(i);
                }
                // find minimum cost edge
                int u, v;
                int minU = 0, minV = 0;
                double minEdgeCost = 0xffff;
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
//                    cout << "minimal ST's root is factoryNode" << endl;
                    rootIsFactory = true;
                    containRoot = true;
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
                totalCost += t.seek(newRootNode, dist, containRoot, tourAns, depNum);
                for (auto childIt = miniTree->child.begin(); childIt != miniTree->child.end();) {
                    auto child = *childIt;
                    if (find(nodes.begin(), nodes.end(), child) != nodes.end()) {
                        if (child->index == minV)cout << "child contains minV" << endl;
                        for (auto newRootChildIt = newRootNode->child.begin();
                             newRootChildIt != newRootNode->child.end();) {
                            auto newRootChild = *newRootChildIt;
                            if (child->index == newRootChild->index) {
                                newRootNode->child.erase(newRootChildIt);
                            } else {
                                newRootChildIt++;
                            }
                        }
                        miniTree->child.erase(childIt);
                        freeTree(child);
                    } else {
                        childIt++;
                    }
                }
                if (rootIsFactory)free(minVNode);
                freeTree(newRootNode);
                if (containRoot) {
                    miniTree->demand = 0;
                    treeNode *parents = miniTree;
                    while (parents != nullptr) {
                        parents->totalDemand -= score;
                        parents = parents->parent;
                    }
                } else {
                    treeNode *parents = miniTree;
                    while (parents != nullptr) {
                        parents->totalDemand -= score;
                        parents = parents->parent;
                    }
                }
                // remove rootNode if no child
                if (miniTree->totalDemand == 0) {
//                    cout << "miniTree is zero demand!" << endl;
                    for (auto it = miniTree->parent->child.begin(); it != miniTree->parent->child.end(); it++) {
                        if (*it == miniTree) {
                            miniTree->parent->child.erase(it);
                            freeTree(miniTree);
                            break;
                        }
                    }
                }
            }
        }
    }
    // check correct or not
    if (tree->totalDemand != 0)cout << "after totalDemand != 0" << endl;
    double totalDistance = 0.0;
    for (int i = 0; i < depNum; i++) {
        if (openFlag[i]) {
            totalDistance += openCost[i];
        }
    }
    set<int> allCustomers;
    int num = 0;
    int totalDemandSum = 0;
    for (auto tour:tourAns) {
        double distance = 0.0;
        int demandSum = 0;
        for (int i = 0; i < tour.size(); i++) {
            int from = tour[i];
            int to;
            if (i == tour.size() - 1) to = tour[0];
            else to = tour[i + 1];
            if (from > depNum) {
                demandSum += demand[from - depNum - 1];
                if (allCustomers.find(from) != allCustomers.end()) {
                    cout << "Duplicate To : " << from << "!!!!!!!!" << endl;
                }
                allCustomers.insert(from);
            }
            distance += dist[from][to];
        }
        totalDemandSum += demandSum;
        totalDistance += distance;
        if (constant::openLog)
            cout << "tour " << ++num << " demandSum = " << demandSum << " distance = " << distance << endl;
    }
    if (constant::openLog) cout << "total Demand: " << totalDemandSum << endl;
    for (int i = depNum + 1; i < n; i++) {
        if (allCustomers.find(i) == allCustomers.end())
            cout << "customer " << i - depNum << " is not satisfied!" << endl;
    }
    if (constant::openLog) cout << "Total Distance: " << totalDistance << endl;
    if (totalDistance != totalCost)cout << "Total Cost != Total Distance!" << endl;
    cout << totalCost << "\t";
    if (constant::openLog) cout << "---------------------------------" << endl;

}