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
#include <cmath>
#include <time.h>
#include <boost/heap/fibonacci_heap.hpp>

#include "Tools.h"


using namespace std;




double getDistance(double lat1, double lon1, double lat2, double lon2)
{
	/*
	C = sin(MLatA)*sin(MLatB)*cos(MLonA-MLonB) + cos(MLatA)*cos(MLatB)
	Distance = R*Arccos(C)*Pi/180
	*/
	double pi = 3.14159265358979323846;
	double Mlat1 = lat1 * pi / 180;
	double Mlat2 = lat2 * pi / 180;
	double Mlon1 = lon1 * pi / 180;
	double Mlon2 = lon2 * pi / 180;
	double C = sin(Mlat1) * sin(Mlat2) + cos(Mlon1 - Mlon2) * cos(Mlat1) * cos(Mlat2);
	double dis = acos(C) * 180 * 60 * 1.1515 * 1609.344 / pi;
	return dis;
}

double EQFG::getSpatialSim(int qid) // the user's location is stored in a global varible Ulat, Ulon
{
	double ret = 0.0;
	map<int, double> & locMap = this->query2loc_[qid];
	for (map<int, double>::iterator i = locMap.begin(); i != locMap.end(); ++i) {
		if (getDistance(Ulat, Ulon, loc2cor_[i->first].first, loc2cor_[i->first].second) <= DIS_THRESHOLD) {
			ret += i->second;
		}
	}
	return ret;
}

double EQFG::spatialAdjustWeight(int qid, double w, double beta) 
{
	return beta * w + (1 - beta) * getSpatialSim(qid);
}


vector<pair<int, double>> EQFG::PPR_BCA(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType)
{
	// edgeType = 0 for entity PPR
	// edgeType = 1 for query PPR
	BoundHeap heap(nodes.size());
	double activeInk = 0.0;
	map<int, double> result;

	// initialize the heap
	for (map<int, double>::iterator i = initialInk.begin(); i != initialInk.end(); ++i) {
		heap.push(*i);
		activeInk += i->second;
	}

	while (heap.size() > 0){ //&& activeInk > PPR_EPS) {
		//cerr << heap.size() << endl;
		pair<int, double> topItem = heap.pop();
		if (topItem.second < PPR_IGNORE_INK) {
			break;
		}
		double increaseInk = topItem.second * alpha;
		if (result.find(topItem.first) == result.end()) {
			result[topItem.first] = 0.0;
		}
		result[topItem.first] += increaseInk;
		activeInk -= increaseInk;
		double distributedInk = (1.0 - alpha) * topItem.second;
		
		// Bug checking
		//if (topItem.first > nodes.size()) continue;

		vector<EQFG_Edge> & edges = nodes[topItem.first].toEntityEdges_;
		if (edgeType == 1) {
			edges = nodes[topItem.first].toQueryEdges_;
		}
		for (int i = 0; i < edges.size(); ++i) {
			double tw = edges[i].w_;
			if (edgeType == 1) {
				// Adjust the weights for query2query edges
				tw = spatialAdjustWeight(edges[i].eid_, tw, beta);
			}
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

vector<pair<int, double>> EQFG::PPR_BCA_lazy(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType)
{
	// edgeType = 0 for entity PPR
	// edgeType = 1 for query PPR
	BoundHeap heap(nodes.size());
	double activeInk = 0.0;
	map<int, double> result;

	// initialize the heap
	for (map<int, double>::iterator i = initialInk.begin(); i != initialInk.end(); ++i) {
		heap.push(*i);
		activeInk += i->second;
	}

	vector<double> inkBuffer(nodes.size(), 0.0);

	while (heap.size() > 0) { //&& activeInk > PPR_EPS) {
		pair<int, double> topItem = heap.pop();
		if (topItem.second < PPR_IGNORE_INK) {
			break;
		}
		double increaseInk = topItem.second * alpha;
		if (result.find(topItem.first) == result.end()) {
			result[topItem.first] = 0.0;
		}
		result[topItem.first] += increaseInk;
		activeInk -= increaseInk;
		double distributedInk = (1.0 - alpha) * topItem.second;

		// Bug checking
		//if (topItem.first > nodes.size()) continue;

		vector<EQFG_Edge> & edges = nodes[topItem.first].toEntityEdges_;
		if (edgeType == 1) {
			edges = nodes[topItem.first].toQueryEdges_;
		}
		for (int i = 0; i < edges.size(); ++i) {
			double tw = edges[i].w_;
			if (edgeType == 1) {
				// Adjust the weights for query2query edges
				tw = spatialAdjustWeight(edges[i].eid_, tw, beta);
			}
			double addInk = distributedInk * tw;

			// lazy update
			inkBuffer[edges[i].eid_] += addInk;
			if (inkBuffer[edges[i].eid_] >= PPR_IGNORE_INK){
				heap.push(make_pair(edges[i].eid_, inkBuffer[edges[i].eid_]));
				inkBuffer[edges[i].eid_] = 0.0;
			}
			
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


void EQFG::loadLocation(const string locPath)
{
	string loc2corPath = locPath + "loc2cor";
	ifstream loc2corIn(loc2corPath.c_str(), ios::in);
	string line;
	cerr << "Starting reading " << loc2corPath << endl;
	while(getline(loc2corIn, line))	{
		vector<string> strs = split(line);
		if (strs.size() != 3)
			continue;
		loc2id_[strs[0]] = locations_.size();
		loc2cor_.push_back(make_pair(atof(strs[1].c_str()), atof(strs[2].c_str())));
		locations_.push_back(strs[0]);
	}
	loc2corIn.close();

	string query2locPath = locPath + "query2loc.txt";
	cerr << "Starting reading " << query2locPath << endl;
	ifstream query2locIn(query2locPath.c_str(), ios::in);
	while (getline(query2locIn, line)) {
		vector<string> strs = split(line);
		map<int, double> tempMap;
		int sum = 0;
		for (int i = 1; i < strs.size(); i += 2) {
			int count = atoi(strs[i + 1].c_str());
			tempMap[loc2id_[strs[i]]] = count;
			sum += count;
		}
		for (map<int, double>::iterator i = tempMap.begin(); i != tempMap.end(); ++i) {
			i->second /= sum;
		}
		query2loc_[query2id_[strs[0]]] = tempMap;
	}
	query2locIn.close();
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
		if (queries_[sid] == "-")
			continue;
		for (int i = 1; i < strs.size(); i += 2) {
			int eid = atoi(strs[i].c_str());
			if (queries_[eid] == "-")
				continue;
			EQFG_Edge tempEdge(sid, eid, atof(strs[i + 1].c_str()));
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
		if (sid > entities_.size()) continue;
		for (int i = 1; i + 1 < strs.size(); i += 2) {
			double w = atof(strs[i + 1].c_str());
			// Ignore too small weights
			if(w < LOAD_WEIGHT_IGNORE) 
				continue;
			if(atoi(strs[i].c_str()) > queries_.size()) 
				continue;
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			ENodes_[sid].toQueryEdges_.push_back(tempEdge);
			QNodes_[atoi(strs[i].c_str())].toEntityEdges_.push_back(tempEdge);
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
	return PPR_BCA(QNodes_, ink, EQFG_PPR_QUERY_ALPHA, 1.0, k_, 1);
}

vector<pair<int, double> > EQFG::rec_EQFG(int qid)
{
	// The first PPR
	map<int, double> eink;
	for (int i = 0; i < QNodes_[qid].toEntityEdges_.size(); ++i) {
		if(QNodes_[qid].toEntityEdges_[i].sid_ > entities_.size()) continue;
		eink[QNodes_[qid].toEntityEdges_[i].sid_] = 1.0 / QNodes_[qid].toEntityEdges_.size();
	}

	//vector<pair<int, double>> eidWeights = PPR_BCA(ENodes_, eink, EQFG_PPR_ENTITY_ALPHA, 1.0, NUMOFRELATEDENTITY, 0);
	vector<pair<int, double>> eidWeights = PPR_BCA_lazy(ENodes_, eink, EQFG_PPR_ENTITY_ALPHA, 1.0, NUMOFRELATEDENTITY, 0);
	
	// The second PPRs
	map<int, double> qink;
	for (int i = 0; i < eidWeights.size(); ++i) {
		int eid = eidWeights[i].first;
		double w_e = eidWeights[i].second;
		for (int j = 0; j < ENodes_[eid].toQueryEdges_.size(); ++j) {
			int id = ENodes_[eid].toQueryEdges_[j].eid_;
			if (qink.find(id) == qink.end()) {
				qink[id] = 0.0;
			}
			qink[id] += ENodes_[eid].toQueryEdges_[j].w_ * w_e * GAMMA;
		}
	}
	// Combine QFG with EQFG
	if (qink.find(qid) == qink.end()) {
		qink[qid] = 0.0;
	}
	qink[qid] += 1.0 - GAMMA;

	//return PPR_BCA(QNodes_, qink, EQFG_PPR_QUERY_ALPHA, 1.0, k_, 1);
	return PPR_BCA_lazy(QNodes_, qink, EQFG_PPR_QUERY_ALPHA, 0, k_, 1);
}

void EQFG::rec_QFG_fromfile(string inPath, string outPath)
{
	cerr << "Start running QFG reccommendation." << endl;
	clock_t t1 = clock();
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		//cerr << line << endl;
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
	clock_t t2 = clock();
	cerr << "QFG recommendation takes " << (t2 - t1 + 0.0) / CLOCKS_PER_SEC << "seconds" << endl;
}

void EQFG::rec_EQFG_fromfile(string inPath, string outPath)
{
	cerr << "Start running EQFG reccommendation." << endl;
	clock_t t1 = clock();
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		cerr << line << endl;
		if (query2id_.find(line) != query2id_.end()) {
			int qid = query2id_[line];
			vector<pair<int, double> > ret = rec_EQFG(qid);
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
	clock_t t2 = clock();
	cerr << "EQFG recommendation takes " << (t2 - t1 + 0.0) / CLOCKS_PER_SEC << "seconds" << endl;
}
