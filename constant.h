#ifndef TTRPCMAKE_CONSTANT_H
#define TTRPCMAKE_CONSTANT_H

#include <string>

const int N = 5050;
const int M = 2e5 + 10;
const int MAX = 0xffff;

class constant {
public:
    static std::string instanceFileName;
    static bool plusDelta;
    static bool openLog;
    static bool shortCutInPS;
    static bool shortCutCycleInPS;
    static double greedyAlpha;
};

#endif //TTRPCMAKE_CONSTANT_H
