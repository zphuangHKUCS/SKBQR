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

#define NUMOFRELATEDENTITY 10
#define GAMMA 0.5

#define EQFG_PPR_ENTITY_ALPHA 0.5
#define EQFG_PPR_QUERY_ALPHA 0.5

#define PPR_EPS 0.001
#define PPR_IGNORE_INK 0.001
#define LOAD_WEIGHT_IGNORE 0.001

#define DIS_THRESHOLD 10000 // range equals 10km
//#define PRO_THRESHOLD 0.1   // the probability that falling within the range no less than 0.1

#define MAXCLICK 257402


class EQFG_Edge
{
public:
    int sid_, eid_;
    float w_;
    EQFG_Edge(int id1, int id2, float w);
    EQFG_Edge(const EQFG_Edge & e);
};

class EQFG_Node
{
public:
    int id_;
    vector<EQFG_Edge> toQueryEdges_; // outedges
    //vector<EQFG_Edge> inQueryEdges_; // inedges
    vector<EQFG_Edge> toEntityEdges_;
	vector<EQFG_Edge> toDocEdges_;
	map<pair<int, int>, map<int, float>> p2loc_;
	
	EQFG_Node(int id);

};




class EQFG
{
private:
	// The location of the user
	// HK  22.2833 114.15
	// California 40.0689	-79.8732
	// New York 40.7504	-73.9963
	// harverd 38.922558, -77.019416
	// Boston 42.3706	-71.027
	
	float Ulat, Ulon;
	//int UlocID = 5378;  // boston
	int UlocID; // new york

	

	vector<pair<int, double>> PPR_BCA(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType = 0);
	vector<pair<int, double>> PPR_BCA_lazy(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType = 0);
	vector<pair<int, double>> PPR_BCA_lazy_cache(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType = 0);
	double spatialAdjustWeight(int qid, double w, double beta, vector<double> & spCache);
	double getSpatialSim(int qid);
	double getSpatialSim_p(int qid);

public:
	int k_;
	vector<pair<int, double> > rec_QFG(int qid);
	vector<pair<int, double> > rec_EQFG(int qid);

public:
    
	//map<pair<int, int>, vector<int>> partition_;
	vector<pair<int, int>> loc2partition_;

    map<string, int> query2id_;
    map<string, int> entity2id_;
	map<string, int> term2id_;

	map<string, int> loc2id_;
	vector< pair<float, float> > loc2cor_;
	vector<string> locations_;

	map<int, map<int, float>> query2loc_;
    //map<string, int> doc2id_;
	
	//map<int, map<int, double>> entity2docPro_;
	//map<int, map<int, double>> doc2entityPro_;
	//map<int, map<int, double>> doc2queryPro_;
	
    vector<EQFG_Node> QNodes_;
    vector<EQFG_Node> ENodes_;
	vector<EQFG_Node> TNodes_;
    
    vector<string> queries_;
    //vector<string> entities_;
	//vector<string> documents_;
    
    
	EQFG(string indexPath, int k = 25);
 
    void saveToFiles(string dirPath);

	void rec_QFG_fromfile(string inputPath, string outPath);
	void rec_EQFG_fromfile(string inputPath, string outPath);

	void loadLocation(const string locDir);
	void loadTerm(const string queryDir);
	void loadQuery(string indexPath);
	void loadEntity(string indexPath);
};

class DQG
{
private:
	float Ulat, Ulon;
	int k_;
	double query2docWeight(int qid, int did, double w, double beta);
	double doc2queryWeight(int did, int qid, double w, double beta);
	void loadQuery(string indexPath);
	void loadDoc(string indexPath);
	vector<pair<int, double>> PPR_BCA(map<int, double> & initialInk, double alpha, double beta, int k);
	vector<pair<int, double> > rec_DQG(int qid);
public:
	map<string, int> query2id_;
	vector<EQFG_Node> QNodes_;
	vector<string> queries_;

	map<string, int> doc2id_;
	vector<EQFG_Node> DNodes_;
	vector<string> docs_;
	vector< pair<float, float>> doc2cor_;


	DQG(string indexPath, int k = 25);
	void loadLocation(const string locDir);

	void rec_DQG_fromfile(string inputPath, string outPath);
};
#endif /* defined(__EQFG__EQFG__) */
