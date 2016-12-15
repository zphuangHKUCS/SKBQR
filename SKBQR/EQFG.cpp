//
//  EQFG.cpp
//  EQFG
//
//  Created by 黄智鹏 on 16/4/10.
//  Copyright (c) 2016年 黄智鹏. All rights reserved.
//

#include "EQFG.h"
#include "BoundHeap.h"
#include <cstring>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/heap/fibonacci_heap.hpp>

#include "Tools.h"

using namespace std;



EQFG_Node::EQFG_Node(int id): id_(id)
{

}

EQFG_Edge::EQFG_Edge(int id1, int id2, double w): sid_(id1), eid_(id2), w_(w)
{

}
EQFG_Edge::EQFG_Edge(const EQFG_Edge & e)
{
    sid_ = e.sid_;
    eid_ = e.eid_;
    w_ = e.w_;
}




void EQFG::saveToFiles(string dirPath)
{
	cerr << "start saveing the query log to " << dirPath << endl;
    ofstream query2idOut(dirPath + "query2id.txt", ios::out);
    for(int i = 0; i < queries_.size(); ++i){
        query2idOut << queries_[i] << '\t' << i << endl;
    }
    query2idOut.close();
    
    ofstream entity2idOut(dirPath + "entity2id.txt", ios::out);
    for(int i = 0; i < entities_.size(); ++i){
        entity2idOut << entities_[i] << '\t' << i << endl;
    }
    entity2idOut.close();
    
    ofstream query2query_w_out(dirPath + "query2query_w.txt", ios::out);
    for(int i = 0; i < QNodes_.size(); ++i){
        query2query_w_out << i;
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            query2query_w_out << '\t' << QNodes_[i].toQueryEdges_[j].eid_ << '\t' << QNodes_[i].toQueryEdges_[j].w_;
        }
        query2query_w_out << endl;
    }
    query2query_w_out.close();
    
    ofstream entity2query_w_out(dirPath + "entity2query_w.txt", ios::out);
    for(int i = 0; i < ENodes_.size(); ++i){
        entity2query_w_out << i;
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            entity2query_w_out << '\t' << ENodes_[i].toQueryEdges_[j].eid_ << '\t' << ENodes_[i].toQueryEdges_[j].w_;
        }
        entity2query_w_out << endl;
    }
    entity2query_w_out.close();
    
    ofstream entity2entity_w_out(dirPath + "entity2entity_w.txt", ios::out);
    for(int i = 0; i < ENodes_.size(); ++i){
        entity2entity_w_out << i;
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            entity2entity_w_out << '\t' << ENodes_[i].toEntityEdges_[j].eid_ << '\t' << ENodes_[i].toEntityEdges_[j].w_;
        }
        entity2entity_w_out << endl;
    }
    entity2entity_w_out.close();
}






EQFG::EQFG(string indexPAth)
{
	string line;
	cerr << "start loading the query nodes." << endl;
	string temps = indexPAth + "query2id.txt";
	ifstream query2idIn(temps.c_str(), ios::in);
	while (getline(query2idIn, line)) {
		vector<string> strs = split(line, "\t");
		QNodes_.push_back(EQFG_Node(queries_.size()));
		query2id_[strs[0]] = queries_.size();
		queries_.push_back(strs[0]);
	}
	query2idIn.close();

	cerr << "start loading the entity nodes." << endl;
	temps = indexPAth + "entity2id.txt";
	ifstream entity2idIn(temps.c_str(), ios::in);
	while (getline(entity2idIn, line)) {
		vector<string> strs = split(line, "\t");
		ENodes_.push_back(EQFG_Node(entities_.size()));
		entity2id_[strs[0]] = entities_.size();
		entities_.push_back(strs[0]);
	}
	entity2idIn.close();
	cerr << "start loading the query2query edges." << endl;
	string tempPath = indexPAth + "query2query_w.txt";
	ifstream query2query_w_in(tempPath.c_str(), ios::in);
	while (getline(query2query_w_in, line)) {
		vector<string> strs = split(line, "\t");
		int sid = atoi(strs[0].c_str());
		for (int i = 1; i < strs.size(); i += 2) {
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			QNodes_[sid].toQueryEdges_.push_back(tempEdge);
		}
	}
	query2idIn.close();
	cerr << "start loading the entity2query edges." << endl;
	tempPath = indexPAth + "entity2query_w.txt";
	ifstream entity2queryIn(tempPath.c_str(), ios::in);
	while (getline(entity2queryIn, line)) {
		vector<string> strs = split(line, "\t");
		int sid = atoi(strs[0].c_str());
		for (int i = 1; i < strs.size(); i += 2) {
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			ENodes_[sid].toQueryEdges_.push_back(tempEdge);
			QNodes_[tempEdge.eid_].toEntityEdges_.push_back(tempEdge);
		}
	}
	entity2queryIn.close();
	cerr << "start loading the entity2entity edges" << endl;
	tempPath = indexPAth + "entity2entity_w.txt";
	ifstream entity2entityIn(tempPath.c_str(), ios::in);
	while (getline(entity2entityIn, line)) {
		vector<string> strs = split(line, "\t");
		int sid = atoi(strs[0].c_str());
		double sum = 0.0;
		for (int i = 1; i < strs.size(); i += 2) {
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			sum += atof(strs[i + 1].c_str());
			ENodes_[sid].toEntityEdges_.push_back(tempEdge);
		}
		for (int i = 0; i < ENodes_[sid].toEntityEdges_.size(); i++) {
			ENodes_[sid].toEntityEdges_[i].w_ /= sum;
		}
	}
	entity2entityIn.close();
	cerr << "end of building the graph." << endl;
}