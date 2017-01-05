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
#include <sstream>
using namespace std;



void tempsplit(string &s, string delim, std::vector<string> &elems) {
	stringstream ss;
	ss.str(s);
	std::string item;
	while (getline(ss, item, delim[0])) {
		elems.push_back(item);
	}
}

vector<string> split(string &s, string delim = "\t") {
	vector<std::string> elems;
	tempsplit(s, delim, elems);
	return elems;
}

/*
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
*/

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
