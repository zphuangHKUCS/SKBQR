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



EQFG_QNode::EQFG_QNode(int id): id_(id)
{

}
EQFG_ENode::EQFG_ENode(int id): id_(id)
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


EQFG::EQFG(string graph_dataPath, string queryEntityPath, string queryCountPath)
{	
    if(graph_dataPath[graph_dataPath.size() -1] != '/')
	graph_dataPath += "/";
    string line;
    // read the query2id file
    string query2idPath = graph_dataPath + "query-id-map.tsv";
    cerr << "start reading " << query2idPath << endl;
    ifstream query2idIn(query2idPath.c_str(), ios::in);
    int maxID = -1;
    while(getline(query2idIn, line)){
	//cerr << line << endl;
        vector<string> strs = split(line, "\t");
        int id = atoi(strs[1].c_str());
        if(id > maxID)
            maxID = id;
        string query = strs[0];
        query2id_[query] = id;
    }
    query2idIn.close();
    queries_ = vector<string>(maxID + 1);
    for(map<string, int>::iterator i = query2id_.begin(); i != query2id_.end(); ++i){
        queries_[i->second] = i->first;
    }
    // add QNodes
    for(int i = 0; i < queries_.size(); ++i){
        QNodes_.push_back(EQFG_QNode(i));
    }
    // add Q->Q edges and weights
    
    string QQEdgePath = graph_dataPath + "query-rewrite-matrix.tsv";
    //cerr << "start reading " << QQEdgePath << endl;
    ifstream edgeIn(QQEdgePath.c_str(), ios::in);
    while(getline(edgeIn, line)){
        vector<string> strs = split(line, "\t");
        int sid = atoi(strs[0].c_str());
        for(int i = 1; i < strs.size(); i += 2){
            int eid = atoi(strs[i].c_str());
            double w = atof(strs[i+1].c_str());
            EQFG_Edge tempEdge(sid, eid, w);
            QNodes_[sid].toQueryEdges_.push_back(tempEdge);
            //QNodes_[eid].inQueryEdges_.push_back(tempEdge);
        }
    }
    edgeIn.close();
    cerr << "size of some_Query:\t" << QNodes_.size() << endl;
    // read the queryEntity file
    ifstream queryEntityIn(queryEntityPath.c_str(), ios::in);
    cerr << "start reading " << queryEntityPath << endl;
    while(getline(queryEntityIn, line)){
        vector<string> strs = split(line, "\t");
        string query = strs[0];
        if(query2id_.find(query) == query2id_.end()){
            int id = QNodes_.size();
            query2id_[query] = id;
            QNodes_.push_back(EQFG_QNode(id));
            queries_.push_back(query);
        }
        string entity = strs[1];
        if(entity2id_.find(entity) == entity2id_.end()){
            int id = ENodes_.size();
            entity2id_[entity] = id;
            ENodes_.push_back(EQFG_ENode(id));
            entities_.push_back(entity);
        }
        int qid = query2id_[query];
        int eid = entity2id_[entity];
        EQFG_Edge tempEdge(eid, qid, 1.0);
        ENodes_[eid].toQueryEdges_.push_back(tempEdge);
        QNodes_[qid].toEntityEdges_.push_back(tempEdge);
    }
    queryEntityIn.close();
    cerr << "size of all_Query:\t" << QNodes_.size() << endl;
    cerr << "start reading " << queryCountPath << endl;
    // set the weight of query-entity edges
    ifstream queryCountIn(queryCountPath.c_str(), ios::in);
    while(getline(queryCountIn, line)){
        vector<string> strs = split(line, "\t");
        int qid = query2id_[strs[0]];
        int count = atoi(strs[1].c_str());
        QNodes_[qid].count_ = count;
    }
    queryCountIn.close();
    
    for(int i = 0; i < ENodes_.size(); ++i){
        int totalCount = 0;
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            int qid = ENodes_[i].toQueryEdges_[j].eid_;
            totalCount += QNodes_[qid].count_;
        }
        for(int j = 0; j < ENodes_[i].toQueryEdges_.size(); ++j){
            int qid = ENodes_[i].toQueryEdges_[j].eid_;
            ENodes_[i].toQueryEdges_[j].w_ = (double)(QNodes_[qid].count_) / totalCount;
        }
    }
    cerr << "start setting the weight of e2e edges." << endl;
    // set the weight of entity-entity edges
    map<pair<int, int>, double> tempMap;
    for(int i = 0; i < QNodes_.size(); ++i){
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            int qid1 = i;
            int qid2 = QNodes_[i].toQueryEdges_[j].eid_;
            if(QNodes_[qid1].toEntityEdges_.size() == 0 || QNodes_[qid2].toEntityEdges_.size() == 0)
                continue;
            for(int k1 = 0; k1 < QNodes_[qid1].toEntityEdges_.size(); ++k1){
                for(int k2 = 0; k2 < QNodes_[qid2].toEntityEdges_.size(); ++k2){
                    int eid1 = QNodes_[qid1].toEntityEdges_[k1].sid_;
                    int eid2 = QNodes_[qid2].toEntityEdges_[k2].sid_;
                    EQFG_Edge tempEdge(eid1, eid2, 0.0);
                    bool flag = true;
                    for(int kk = 0; kk < ENodes_[eid1].toEntityEdges_.size(); ++ kk){
                        if(ENodes_[eid1].toEntityEdges_[kk].eid_ == tempEdge.eid_){
                            flag = false;
                            break;
                        }
                    }
                    if(flag){
                        ENodes_[eid1].toEntityEdges_.push_back(tempEdge);
                    }
                    
                    tempMap[make_pair(eid1, eid2)] = 1.0;
                }
            }
            
        }
    }
    for(int i = 0; i < QNodes_.size(); ++i){
        for(int j = 0; j < QNodes_[i].toQueryEdges_.size(); ++j){
            int qid1 = i;
            int qid2 = QNodes_[i].toQueryEdges_[j].eid_;
            if(QNodes_[qid1].toEntityEdges_.size() == 0 || QNodes_[qid2].toEntityEdges_.size() == 0)
                continue;
            double tobetimes = 1.0 - QNodes_[i].toQueryEdges_[j].w_ / (QNodes_[qid1].toEntityEdges_.size() * QNodes_[qid2].toEntityEdges_.size());
            for(int k1 = 0; k1 < QNodes_[qid1].toEntityEdges_.size(); ++k1){
                for(int k2 = 0; k2 < QNodes_[qid2].toEntityEdges_.size(); ++k2){
                    int eid1 = QNodes_[qid1].toEntityEdges_[k1].sid_;
                    int eid2 = QNodes_[qid2].toEntityEdges_[k2].sid_;
                    tempMap[make_pair(eid1, eid2)] *= tobetimes;
                }
            }
        }
    }
    for(map<pair<int, int>, double>::iterator i = tempMap.begin(); i != tempMap.end(); ++i){
        i->second = 1.0 - i->second;
    }
    for(int i = 0; i < ENodes_.size(); ++i){
        for(int j = 0; j < ENodes_[i].toEntityEdges_.size(); ++j){
            int id2 = ENodes_[i].toEntityEdges_[j].eid_;
            ENodes_[i].toEntityEdges_[j].w_ = tempMap[make_pair(i, id2)];
        }
    }
    cerr << "end building the graph." << endl;
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
		QNodes_.push_back(EQFG_QNode(queries_.size()));
		query2id_[strs[0]] = queries_.size();
		queries_.push_back(strs[0]);
	}
	query2idIn.close();

	cerr << "start loading the entity nodes." << endl;
	temps = indexPAth + "entity2id.txt";
	ifstream entity2idIn(temps.c_str(), ios::in);
	while (getline(entity2idIn, line)) {
		vector<string> strs = split(line, "\t");
		ENodes_.push_back(EQFG_ENode(entities_.size()));
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