#include "Num.h"

template <typename T> Num<T>::Num() : num(0) {};

template <typename T> Num<T>::Num(int n) : num(n) {};

template <typename T> int Num<T>::getNum(){ 
    return this->num;
}
template<typename T>
std::vector<T> test(int n) {
    return std::vector<T>(n, 0);
}
