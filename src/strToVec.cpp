//
// Created by WANGJade on 21/11/2018.
//

// change a string to vector<double>
#include "strToVec.h"
using namespace std;

void strToVec(string str, vector<double>& fea)
{
    stringstream ss(str);
    string buf;
    vector<double> vec;
    //vector<double> vec(1024);//可以初始化时指定vecrot容器大小
    while (ss >> buf)
        vec.push_back(atof(buf.c_str()));

    fea = vec;
}



