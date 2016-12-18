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
#include "BoundHeap.h"

using namespace std;

#define NUMOFRELATEDENTITY 20
#define GAMMA 0.5

#define PPR_EPS 0.0001
#define PPR_IGNORE_INK 0.00001
#define LOAD_WEIGHT_IGNORE 0.00001


class EQFG_Edge
{
public:
    int sid_, eid_;
    double w_;
    EQFG_Edge(int id1, int id2, double w);
    EQFG_Edge(const EQFG_Edge & e);
};

class EQFG_Node
{
public:
    int id_;
    vector<EQFG_Edge> toQueryEdges_; // outedges
    //vector<EQFG_Edge> inQueryEdges_; // inedges
    vector<EQFG_Edge> toEntityEdges_;
    EQFG_Node(int id);
};




class EQFG
{
public:
	int k_;
	vector<pair<int, double> > rec_QFG(int qid);
	vector<pair<int, double> > rec_EQFG(int qid);

public:
    
    map<string, int> query2id_;
    map<string, int> entity2id_;
    //map<string, int> doc2id_;
	
	//map<int, map<int, double>> entity2docPro_;
	//map<int, map<int, double>> doc2entityPro_;
	//map<int, map<int, double>> doc2queryPro_;
	
    vector<EQFG_Node> QNodes_;
    vector<EQFG_Node> ENodes_;
    
    vector<string> queries_;
    vector<string> entities_;
	//vector<string> documents_;
    
    
	EQFG(string indexPath, int k = 5);
 
    void saveToFiles(string dirPath);

	void rec_QFG_fromfile(string inputPath, string outPath);
	void rec_EQFG_fromfile(string inputPath, string outPath);


};


#endif /* defined(__EQFG__EQFG__) */
