#define LoadData_H
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stack>
#include <map>
#include <math.h>
#include <stdlib.h>

template<typename T> 
using Data = std::vector<std::vector<T>>;

template <typename T>
using RowData = std::vector<T>;

template <typename T>
using ColData = std::vector<T>;

template <typename T>
void loadData(Data<T> &data, const char *infile);

