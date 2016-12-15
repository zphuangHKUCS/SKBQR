//
//  Tools.h
//  EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//

#ifndef EQFG_Tools_h
#define EQFG_Tools_h

#include <vector>
#include <cstring>
#include <iostream>
using namespace std;


vector<string> split(string str, string separator = "\t")
{
    vector<string> result;
    int cutAt;
    while( (cutAt = str.find_first_of(separator)) != str.npos )
    {
        if(cutAt > 0)
        {
            result.push_back(str.substr(0, cutAt));
        }else{
            result.push_back("");
        }
        str = str.substr(cutAt + 1);
    }
    if(str.length() > 0)
    {
        result.push_back(str);
    }else{
        result.push_back("");
    }
    return result;
}

vector<string> mysplit(string str, string separator = "\t")
{
    vector<string> result;
    int cutAt;
    while( (cutAt = str.find_first_of(separator)) != str.npos )
    {
        if(cutAt > 0)
        {
            result.push_back(str.substr(0, cutAt));
        }else{
            result.push_back("");
        }
        str = str.substr(cutAt + 1);
    }
    if(str.length() > 0)
    {
        result.push_back(str);
    }else{
        result.push_back("");
    }
    return result;
}



#endif
