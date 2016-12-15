//
//  EQFG.h
//  EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//

#ifndef __EQFG__EQFG__
#define __EQFG__EQFG__

#include <stdio.h>
#include <cstring>
#include <vector>
#include <map>
#include <iostream>

using namespace std;


class EQFG_Edge
{
public:
    int sid_, eid_;
    double w_;
    EQFG_Edge(int id1, int id2, double w);
    EQFG_Edge(const EQFG_Edge & e);
};

class EQFG_QNode
{
public:
    int id_;
    vector<EQFG_Edge> toQueryEdges_; // outedges
    //vector<EQFG_Edge> inQueryEdges_; // inedges
    vector<EQFG_Edge> toEntityEdges_;
    EQFG_QNode(int id);
};


class EQFG_ENode
{
public:
    int id_;
    
    vector<EQFG_Edge> toQueryEdges_;
    vector<EQFG_Edge> toEntityEdges_;
    //vector<EQFG_Edge> inEntityEdges_;
    EQFG_ENode(int id);
};

class EQFG
{
public:
    
    map<string, int> query2id_;
    map<string, int> entity2id_;
    //map<string, int> doc2id_;
	
	//map<int, map<int, double>> entity2docPro_;
	//map<int, map<int, double>> doc2entityPro_;
	//map<int, map<int, double>> doc2queryPro_;
	
    vector<EQFG_QNode> QNodes_;
    vector<EQFG_ENode> ENodes_;
    
    vector<string> queries_;
    vector<string> entities_;
	//vector<string> documents_;
    
    
    EQFG(string queryFollowPath, string queryEntityPath, string queryCountPath);
	EQFG(string indexPath);
 
    void saveToFiles(string dirPath);
};


#endif /* defined(__EQFG__EQFG__) */
