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

double spatialAdjustWeight() {
	return 0.0;
}

vector<pair<int, double>> PPR_BCA(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType = 0)
{
	// edgeType = 0 for entity PPR
	// edgeType = 1 for query PPR
	BoundHeap heap(9999999);
	double activeInk = 1.0;
	map<int, double> result;

	// initialize the heap
	for (map<int, double>::iterator i = initialInk.begin(); i != initialInk.end(); ++i) {
		heap.push(*i);
	}
	while (heap.size() > 0 && activeInk > 0.001) {
		pair<int, double> topItem = heap.pop();
		cerr << activeInk << endl;
		if (topItem.second < 0.00001) {
			continue;
		}
		double increaseInk = topItem.second * alpha;
		if (result.find(topItem.first) == result.end()) {
			result[topItem.first] = 0.0;
		}
		result[topItem.first] += increaseInk;
		activeInk -= increaseInk;
		double distributedInk = (1.0 - alpha) * topItem.second;
		vector<EQFG_Edge> & edges = nodes[topItem.first].toEntityEdges_;
		if (edgeType == 1) {
			edges = nodes[topItem.first].toQueryEdges_;
		}
		for (int i = 0; i < edges.size(); ++i) {
			// No spatial adjustment now
			double tw = edges[i].w_;
			double addInk = distributedInk * tw;
			heap.push(make_pair(edges[i].eid_, addInk));
		}
	}
	// find the result
	BoundHeap topk(k);
	for (map<int, double>::iterator i = result.begin(); i != result.end(); ++i) {
		topk.push(make_pair(i->first, -i->second));
	}
	vector<pair<int, double> > reverseRet, ret;
	while (topk.size() > 0) {
		pair<int, double> item = topk.pop();
		reverseRet.push_back(make_pair(item.first, -item.second));
	}
	for (int i = 0; i < reverseRet.size(); ++i) {
		ret.push_back(reverseRet[reverseRet.size() - 1 - i]);
	}
	return ret;
}


EQFG_Node::EQFG_Node(int id): id_(id){}

EQFG_Edge::EQFG_Edge(int id1, int id2, double w): sid_(id1), eid_(id2), w_(w){}

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






EQFG::EQFG(string indexPAth, int k): k_(k)
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

vector<pair<int, double> > EQFG::rec_QFG(int qid)
{
	map<int, double> ink;
	ink[qid] = 1.0;
	return PPR_BCA(QNodes_, ink, 0.3, 1.0, k_, 1);
}

void EQFG::rec_QFG_fromfile(string inPath, string outPath)
{
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		cerr << line << endl;
		if (query2id_.find(line) != query2id_.end()) {
			int qid = query2id_[line];
			vector<pair<int, double> > ret = rec_QFG(qid);
			out << line;
			for (int i = 0; i < ret.size(); ++i) {
				out << '\t' << queries_[ret[i].first] << '\t' << ret[i].second;
			}
			out << endl;
		}
		else {
			out << line << endl;
		}
	}
	in.close();
	out.close();
}