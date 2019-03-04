#ifndef Num_h
#define Num_h
#include <vector>
template <typename T>
class Num {
private:
    int num;
public:
    Num();
    Num(int n);
    int getNum();
    std::vector<T> test(int n);
};
#endif