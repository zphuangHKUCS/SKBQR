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




double getDistance(float lat1, float lon1, float lat2, float lon2)
{
	/*
	C = sin(MLatA)*sin(MLatB)*cos(MLonA-MLonB) + cos(MLatA)*cos(MLatB)
	Distance = R*Arccos(C)*Pi/180
	*/
	float pi = 3.14159265358979323846;
	float Mlat1 = lat1 * pi / 180;
	float Mlat2 = lat2 * pi / 180;
	float Mlon1 = lon1 * pi / 180;
	float Mlon2 = lon2 * pi / 180;
	float C = sin(Mlat1) * sin(Mlat2) + cos(Mlon1 - Mlon2) * cos(Mlat1) * cos(Mlat2);
	float dis = acos(C) * 20014123.8528 / pi;
	return dis;
}
double getDistance_app(double lat1, double lon1, double lat2, double lon2)
{
	// return approximate distance
	double dis = 0.0;
	dis += (lat1 - lat2) * (lat1 - lat2);
	dis += (lon2 - lon1) * (lon2 - lon1);
	dis = sqrt(dis);
	dis *= 111.31955 * 1000;
	return dis;
}
double EQFG::getSpatialSim(int qid) // the user's location is stored in a global varible Ulat, Ulon
{
	double ret = 0.0;
	for (map<pair<int, int>, map<int, float>>::iterator iter = QNodes_[qid].p2loc_.begin(); iter != QNodes_[qid].p2loc_.end(); ++iter) {
		map<int, float> & locMap = iter->second;
		//cerr << locMap.size() << endl;
		for (map<int, float>::iterator i = locMap.begin(); i != locMap.end(); ++i) {
			if (getDistance(Ulat, Ulon, loc2cor_[i->first].first, loc2cor_[i->first].second) <= r_) {
				ret += i->second;
			}
		}
	}
	cerr << "spatial sim is: " << ret << endl;
	return ret;
}
double EQFG::getSpatialSim_p(int qid) // use the partition to compute
{
	if (QNodes_[qid].p2sims_.find(loc2partition_[UlocID]) == QNodes_[qid].p2sims_.end()) {
		//cerr << 0.0 << endl;
		return 0.0;
	}
	//cerr << QNodes_[qid].p2sims_[loc2partition_[UlocID]] << endl;
	return QNodes_[qid].p2sims_[loc2partition_[UlocID]];

	////////////
	double ret = 0.0;
	map<int, float> & locMap = QNodes_[qid].p2loc_[this->loc2partition_[UlocID]];
	
	for (map<int, float>::iterator i = locMap.begin(); i != locMap.end(); ++i) {
		if (getDistance(Ulat, Ulon, loc2cor_[i->first].first, loc2cor_[i->first].second) <= r_) {
			ret += i->second;
		}
	}
	//cerr << "spatial sim is: " << ret << endl;
	return ret;
}

double EQFG::spatialAdjustWeight(int qid, double w, double beta, vector<double> & spCache) 
{
	//return beta * w + (1 - beta) * getSpatialSim(qid);
	if (spCache.size() == 0) {
		return beta * w + (1 - beta) * getSpatialSim(qid);
		//return beta * w + (1 - beta) * getSpatialSim_p(qid);
	}
	else {
		if (spCache[qid] < 0) {
			double sptialSim = getSpatialSim(qid);
			//double sptialSim = getSpatialSim_p(qid);
			spCache[qid] = sptialSim;
		}
		//cerr << spCache[qid] << endl;
		return beta * w + (1 - beta) * spCache[qid];
	}
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
		vector<double> weights;
		double wSum = 0.0;
		if (edgeType == 1) {
			for (int i = 0; i < edges.size(); ++i) {
				vector<double> tempVec;
				double spatialWeight = spatialAdjustWeight(edges[i].eid_, edges[i].w_, beta, tempVec);
				//double spatialWeight = edges[i].w_;
				wSum += spatialWeight;
				weights.push_back(spatialWeight);
			}
			for (int i = 0; i < edges.size(); ++i) {
				weights[i] /= wSum;
			}
		}
		else {
			for (int i = 0; i < edges.size(); ++i) {
				weights.push_back(edges[i].w_);
			}
		}

		for (int i = 0; i < edges.size(); ++i) {
			double addInk = distributedInk * weights[i];
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
	clock_t t1, t2, t3, t4;
	t3 = clock();
	double timeforsptaial = 0.0;
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
		vector<double> weights;
		double wSum = 0.0;
		if (edgeType == 1) {
			for (int i = 0; i < edges.size(); ++i) {
				t1 = clock();
				vector<double> tempVec;
				double spatialWeight = spatialAdjustWeight(edges[i].eid_, edges[i].w_, beta, tempVec);
				t2 = clock();
				timeforsptaial += (t2 - t1 + 0.0) / CLOCKS_PER_SEC;
				//double spatialWeight = edges[i].w_;
				wSum += spatialWeight;
				weights.push_back(spatialWeight);
			}
			for (int i = 0; i < edges.size(); ++i) {
				weights[i] /= wSum;
			}
		}
		else {
			for (int i = 0; i < edges.size(); ++i) {
				weights.push_back(edges[i].w_);
			}
		}

		for (int i = 0; i < edges.size(); ++i) {
			double addInk = distributedInk * weights[i];

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
	t4 = clock();
	//if (edgeType == 1) {
		//cerr << "                      PPR takes:\t" << (t4 - t3 + 0.0) / CLOCKS_PER_SEC << " seconds" << endl;
		//cerr << "sptaial adjusting weights takes:\t" << timeforsptaial << " seconds" << endl;
	//}
	

	return ret;
}

vector<pair<int, double>> EQFG::PPR_BCA_lazy_cache(vector<EQFG_Node> & nodes, map<int, double> & initialInk, double alpha, double beta, int k, int edgeType)
{
	// edgeType = 0 for entity PPR
	// edgeType = 1 for query PPR
	clock_t t1, t2, t3, t4;
	t3 = clock();
	double timeforsptaial = 0.0;
	BoundHeap heap(nodes.size());
	double activeInk = 0.0;
	map<int, double> result;

	// initialize the heap
	for (map<int, double>::iterator i = initialInk.begin(); i != initialInk.end(); ++i) {
		heap.push(*i);
		activeInk += i->second;
	}

	vector<double> inkBuffer(nodes.size(), 0.0);
	vector<double> spatialCache(nodes.size(), -1);

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
		vector<double> weights;
		double wSum = 0.0;
		if (edgeType == 1) {
			for (int i = 0; i < edges.size(); ++i) {
				t1 = clock();
				double spatialWeight = spatialAdjustWeight(edges[i].eid_, edges[i].w_, beta, spatialCache);
				t2 = clock();
				timeforsptaial += (t2 - t1 + 0.0) / CLOCKS_PER_SEC;
				//double spatialWeight = edges[i].w_;
				wSum += spatialWeight;
				weights.push_back(spatialWeight);
			}
			for (int i = 0; i < edges.size(); ++i) {
				weights[i] /= wSum;
			}
		}
		else {
			for (int i = 0; i < edges.size(); ++i) {
				weights.push_back(edges[i].w_);
			}
		}

		for (int i = 0; i < edges.size(); ++i) {
			double addInk = distributedInk * weights[i];

			// lazy update
			inkBuffer[edges[i].eid_] += addInk;
			if (inkBuffer[edges[i].eid_] >= PPR_IGNORE_INK) {
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
	t4 = clock();
	//if (edgeType == 1) {
	//	cerr << "                     EQFG takes:\t" << (t4 - t3 + 0.0) / CLOCKS_PER_SEC << " seconds" << endl;
	//	cerr << "sptaial adjusting weights takes:\t" << timeforsptaial << " seconds" << endl;
	//}


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

	string indexPath = locPath + "partition.txt";
	cerr << "Starting reading " << indexPath << endl;
	loc2partition_ = vector<pair<int, int>>(locations_.size());
	ifstream indexIn(indexPath.c_str());
	while (getline(indexIn, line)) {
		vector<string> strs = split(line);
		int locID = loc2id_[strs[2]];
		int x = atoi(strs[0].c_str());
		int y = atoi(strs[1].c_str());
		pair<int, int> p = make_pair(x, y);
		//if (partition_.find(p) == partition_.end()) {
		//	partition_[p] = vector<int>();
		//}
		//partition_[p].push_back(locID);
		loc2partition_[locID] = p;
	}
	indexIn.close();

	string query2locPath = locPath + "query2loc.txt";
	int enum_q2l = 0;
	cerr << "Starting reading " << query2locPath << endl;
	ifstream query2locIn(query2locPath.c_str(), ios::in);
	while (getline(query2locIn, line)) {
		vector<string> strs = split(line);
		int qid = query2id_[strs[0]];
		map<int, float> tempMap;
		int sum = 0;
		for (int i = 1; i < strs.size(); i += 2) {
			int count = atoi(strs[i + 1].c_str());
			tempMap[loc2id_[strs[i]]] = count;
			sum += count;
		}
		for (map<int, float>::iterator i = tempMap.begin(); i != tempMap.end(); ++i) {
			i->second /= sum;
			enum_q2l += 1;
			pair<int, int> p = loc2partition_[i->first];
			if (QNodes_[qid].p2loc_.find(p) == QNodes_[qid].p2loc_.end()) {
				QNodes_[qid].p2loc_[p] = map<int, float>();
			}
			QNodes_[qid].p2loc_[p][i->first] = i->second;
			if (QNodes_[qid].p2sims_.find(p) == QNodes_[qid].p2sims_.end()) {
				QNodes_[qid].p2sims_[p] = 0.0;
			}
			QNodes_[qid].p2sims_[p] += i->second;
		}
		// remove the query2loc temparorily
		//query2loc_[qid] = tempMap;
		
	}
	query2locIn.close();
	cerr << "#location:" << '\t' << locations_.size() << endl;
	cerr << "#edge_q2l:" << '\t' << enum_q2l << endl;
}

EQFG_Node::EQFG_Node(int id): id_(id){}

EQFG_Edge::EQFG_Edge(int id1, int id2, float w): sid_(id1), eid_(id2), w_(w){}

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
    for(int i = 0; i < ENodes_.size(); ++i){
        entity2idOut << "```" << '\t' << i << endl;
//		entity2idOut << entities_[i] << '\t' << i << endl;
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

void EQFG::loadQuery(string indexPath)
{
	string line;
	cerr << "start loading the query nodes." << endl;
	string temps = indexPath + "query2id.txt";
	ifstream query2idIn(temps.c_str(), ios::in);
	queries_.reserve(9500000);
	QNodes_.reserve(9500000);
	while (getline(query2idIn, line)) {
		vector<string> strs = split(line, "\t");
		QNodes_.push_back(EQFG_Node(queries_.size()));
		query2id_[strs[0]] = queries_.size();
		queries_.push_back(strs[0]);
	}
	query2idIn.close();

	cerr << "start loading the query2query edges." << endl;
	string tempPath = indexPath + "query2query_w.txt";
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
}

void EQFG::loadEntity(string indexPath) {
	string line;
	cerr << "start loading the entity nodes." << endl;
	string temps = indexPath + "entity2id.txt";
	ENodes_.reserve(3500000);
	//entities_.reserve(3500000);
	ifstream entity2idIn(temps.c_str(), ios::in);
	while (getline(entity2idIn, line)) {
		vector<string> strs = split(line, "\t");
		int eid = ENodes_.size();
		ENodes_.push_back(EQFG_Node(eid));
		entity2id_[strs[0]] = eid;
		//entities_.push_back(strs[0]);
	}
	entity2idIn.close();

	cerr << "start loading the entity2query edges." << endl;
	string tempPath = indexPath + "entity2query_w.txt";
	ifstream entity2queryIn(tempPath.c_str(), ios::in);
	while (getline(entity2queryIn, line)) {
		vector<string> strs = split(line, "\t");
		int sid = atoi(strs[0].c_str());
		if (sid > ENodes_.size()) continue;
		for (int i = 1; i + 1 < strs.size(); i += 2) {
			double w = atof(strs[i + 1].c_str());
			// Ignore too small weights
			if (w < LOAD_WEIGHT_IGNORE)
				continue;
			if (atoi(strs[i].c_str()) > queries_.size())
				continue;
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			ENodes_[sid].toQueryEdges_.push_back(tempEdge);
			QNodes_[atoi(strs[i].c_str())].toEntityEdges_.push_back(tempEdge);
		}
	}
	entity2queryIn.close();
	cerr << "start loading the entity2entity edges" << endl;
	tempPath = indexPath + "entity2entity_w.txt";
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
}

void EQFG::loadTerm(string indexPath) {
	string line;
	cerr << "start loading the term nodes." << endl;
	string temps = indexPath + "term2id.txt";
	TNodes_.reserve(3000000);
	ifstream term2idIn(temps.c_str(), ios::in);
	while (getline(term2idIn, line)) {
		vector<string> strs = split(line, "\t");
		int tid = TNodes_.size();
		TNodes_.push_back(EQFG_Node(tid));
		term2id_[strs[0]] = tid;
	}
	term2idIn.close();

	cerr << "start loading the term2query edges." << endl;
	string tempPath = indexPath + "term2query_w.txt";
	ifstream term2queryIn(tempPath.c_str(), ios::in);
	int edgecount = 0;
	while (getline(term2queryIn, line)) {
		vector<string> strs = split(line, "\t");
		int sid = atoi(strs[0].c_str());
		if (sid > TNodes_.size()) continue;
		for (int i = 1; i + 1 < strs.size(); i += 2) {
			double w = atof(strs[i + 1].c_str());
			// Ignore too small weights
			if (w < LOAD_WEIGHT_IGNORE)
				continue;
			if (atoi(strs[i].c_str()) > TNodes_.size())
				continue;
			++edgecount;
			EQFG_Edge tempEdge(sid, atoi(strs[i].c_str()), atof(strs[i + 1].c_str()));
			TNodes_[sid].toQueryEdges_.push_back(tempEdge);
		}
	}
	term2queryIn.close();
	cerr << "#terms            " << '\t' << TNodes_.size() << endl;
	cerr << "#term2query edges " << '\t' << edgecount << endl;
}

EQFG::EQFG(string indexPAth, int k): k_(k)
{
	loadQuery(indexPAth);
	//loadEntity(indexPAth);
	cerr << "end of building the graph." << endl;
	cerr << "#query:" << '\t' << queries_.size() << endl;
}

vector<pair<int, double> > EQFG::rec_QFG(int qid)
{
	map<int, double> ink;
	ink[qid] = 1.0;
	return PPR_BCA(QNodes_, ink, EQFG_PPR_QUERY_ALPHA, 1.0, k_, 1);
}

vector<pair<int, double> > EQFG::rec_EQFG(int qid)
{
	Ulat = loc2cor_[UlocID].first;
	Ulon = loc2cor_[UlocID].second;

	// The first PPR
	map<int, double> eink;
	for (int i = 0; i < QNodes_[qid].toEntityEdges_.size(); ++i) {
		if(QNodes_[qid].toEntityEdges_[i].sid_ > ENodes_.size()) continue;
		eink[QNodes_[qid].toEntityEdges_[i].sid_] = 1.0 / QNodes_[qid].toEntityEdges_.size();
	}

	//vector<pair<int, double>> eidWeights = PPR_BCA(ENodes_, eink, EQFG_PPR_ENTITY_ALPHA, 1.0, NUMOFRELATEDENTITY, 0);
	//vector<pair<int, double>> eidWeights = PPR_BCA_lazy(ENodes_, eink, EQFG_PPR_ENTITY_ALPHA, 1.0, NUMOFRELATEDENTITY, 0);
	vector<pair<int, double>> eidWeights = PPR_BCA_lazy_cache(ENodes_, eink, EQFG_PPR_ENTITY_ALPHA, 1.0, NUMOFRELATEDENTITY, 0);

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

	//return PPR_BCA(QNodes_, qink, EQFG_PPR_QUERY_ALPHA, 0.5, k_, 1);
	//return PPR_BCA_lazy(QNodes_, qink, EQFG_PPR_QUERY_ALPHA, 0.5, k_, 1);
	return PPR_BCA_lazy_cache(QNodes_, qink, EQFG_PPR_QUERY_ALPHA, 0.5, k_, 1);
}

vector<pair<int, double> > EQFG::rec_TQG(int tid, double alpha, double beta)
{
	Ulat = loc2cor_[UlocID].first;
	Ulon = loc2cor_[UlocID].second;
	map<int, double> qink;
	for (int i = 0; i < TNodes_[tid].toQueryEdges_.size(); ++i) {
		qink[TNodes_[tid].toQueryEdges_[i].eid_] = TNodes_[tid].toQueryEdges_[i].w_;
	}
	//return PPR_BCA(QNodes_, qink, alpha, beta, 10 * k_, 1); // return more than k, so we can choose
	//return PPR_BCA_lazy(QNodes_, qink, alpha, beta, 10 * k_, 1); // return more than k, so we can choose
	return PPR_BCA_lazy_cache(QNodes_, qink, alpha, beta, 10 * k_, 1); // return more than k, so we can choose
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
	UlocID = loc2id_["New York"];
	cerr << "Start running EQFG reccommendation." << endl;
	clock_t t1 = clock();
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		//cerr << line << endl;
		vector<string> strs = split(line);
		string query = strs[1];
		if (query2id_.find(query) != query2id_.end()) {
			int qid = query2id_[query];
			vector<pair<int, double> > ret = rec_EQFG(qid);
			out << query;
			for (int i = 0; i < ret.size(); ++i) {
				if (queries_[ret[i].first] == query) {
					continue;
				}
				out << '\t' << queries_[ret[i].first] << '\t' << ret[i].second;
			}
			out << endl;
		}
		else {
			out << query << endl;
		}
	}
	in.close();
	out.close();
	clock_t t2 = clock();
	cerr << "EQFG recommendation takes " << (t2 - t1 + 0.0) / CLOCKS_PER_SEC << "seconds" << endl;
}

void EQFG::rec_TQG_fromfile(string inPath, string outPath, double alpha, double beta, double r)
{
	UlocID = loc2id_["New york"];
	r_ = r;
	cerr << "Start running TQG reccommendation." << endl;
	clock_t t1 = clock();
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		//cerr << line << endl;
		vector<string> strs = split(line);
		string query = strs[1];
		vector<string> terms = split(query, " ");
		map <int, double> result;
		bool firstTime = true;
		for (int i = 0; i < terms.size(); ++i) {
			//cerr << terms[i] << endl;
			if (term2id_.find(terms[i]) == term2id_.end()) {
				continue;
			}
			int tid = term2id_[terms[i]];
			vector<pair<int, double>> tempResult = rec_TQG(tid, alpha, beta);
			//cerr << "Finish rec_TQG for " << terms[i] << endl;
			if (firstTime) {
				firstTime = false;
				for (int j = 0; j < tempResult.size(); ++j) {
					result[tempResult[j].first] = tempResult[j].second;
				}
			}
			else {
				map<int, double> temptemp;
				for (int j = 0; j < tempResult.size(); ++j) {
					if (result.find(tempResult[j].first) != result.end()) {
						temptemp[tempResult[j].first] = result[tempResult[j].first] * tempResult[j].second;
					}
				}
				result.clear();
				result = temptemp;
			}
		}
		BoundHeap topk(k_);
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
		
		if (ret.size() > 0) {
			
			out << query;
			for (int i = 0; i < ret.size(); ++i) {
				if (queries_[ret[i].first] == query) {
					continue;
				}
				out << '\t' << queries_[ret[i].first] << '\t' << ret[i].second;
			}
			out << endl;
		}
		else {
			out << query << endl;
		}
	}
	in.close();
	out.close();
	clock_t t2 = clock();
	cerr << "EQFG recommendation takes " << (t2 - t1 + 0.0) / CLOCKS_PER_SEC << "seconds" << endl;
}


void DQG::loadDoc(string indexPath)
{
	string line;
	cerr << "start loading the query2Doc edges." << endl;
	string temps = indexPath + "doc2que_count.txt";
	ifstream doc2idIn(temps.c_str(), ios::in);

	while (getline(doc2idIn, line)) {
		vector<string> strs = split(line, "\t");
		if (doc2id_.find(strs[0]) == doc2id_.end()) {
			DNodes_.push_back(EQFG_Node(docs_.size()));
			doc2id_[strs[0]] = docs_.size();
			docs_.push_back(strs[0]);
		}
		int did = doc2id_[strs[0]];
		int qid = query2id_[strs[1]];
		double w = atof(strs[2].c_str()) / MAXCLICK;
		DNodes_[did].toQueryEdges_.push_back(EQFG_Edge(did, qid, w));
		QNodes_[qid].toDocEdges_.push_back(EQFG_Edge(qid, did, w));
	}
	doc2idIn.close();
}

void DQG::loadQuery(string indexPath)
{
	string line;
	cerr << "start loading the query nodes." << endl;
	string temps = indexPath + "query2id.txt";
	ifstream query2idIn(temps.c_str(), ios::in);
	queries_.reserve(9500000);
	QNodes_.reserve(9500000);
	while (getline(query2idIn, line)) {
		vector<string> strs = split(line, "\t");
		QNodes_.push_back(EQFG_Node(queries_.size()));
		query2id_[strs[0]] = queries_.size();
		queries_.push_back(strs[0]);
	}
	query2idIn.close();
}

DQG::DQG(string indexPAth, int k) : k_(k)
{
	loadQuery(indexPAth);
	loadDoc(indexPAth);

	cerr << "end of building the graph." << endl;
	cerr << "#query:" << '\t' << queries_.size() << endl;
	cerr << "#doc  :" << '\t' << docs_.size() << endl;
}

void DQG::loadLocation(const string locPath)
{
	map<string, int> loc2id_;
	vector< pair<double, double> > loc2cor_;
	vector<string> locations_;

	string loc2corPath = locPath + "loc2cor";
	ifstream loc2corIn(loc2corPath.c_str(), ios::in);
	string line;
	cerr << "Starting reading " << loc2corPath << endl;
	while (getline(loc2corIn, line)) {
		vector<string> strs = split(line);
		if (strs.size() != 3)
			continue;
		loc2id_[strs[0]] = locations_.size();
		loc2cor_.push_back(make_pair(atof(strs[1].c_str()), atof(strs[2].c_str())));
		locations_.push_back(strs[0]);
	}
	loc2corIn.close();

	doc2cor_ = vector< pair<float, float>>(docs_.size(), make_pair(float(0.0), float(0.0)));
	string url2locPath = locPath + "url2loc";
	cerr << "Starting reading " << url2locPath << endl;

	ifstream url2locIn(url2locPath, ios::in);
	while (getline(url2locIn, line)) {
		vector<string> strs = split(line);
		if (doc2id_.find(strs[0]) == doc2id_.end())
			continue;
		int did = doc2id_[strs[0]];

		int count = 0;
		double lon = 0.0;
		double lat = 0.0;

		for (int i = 1; i < strs.size(); i += 2) {
			int locid = loc2id_[strs[i]];
			int t = atoi(strs[i + 1].c_str());
			count += t;
			lon += loc2cor_[locid].first;
			lat += loc2cor_[locid].second;
		}
		if (count == 0) {
			continue;
		}
		lon /= count;
		lat /= count;
		doc2cor_[did] = make_pair(lon, lat);
	}
	url2locIn.close();

	cerr << "End loading location!" << endl;
}

void DQG::rec_DQG_fromfile(string inPath, string outPath)
{
	//40.7128° N, 74.0059° W new york
	Ulat = 40.7128;
	Ulon = -74.0059;
	cerr << "Start running DQG reccommendation." << endl;
	clock_t t1 = clock();
	ifstream in(inPath.c_str(), ios::in);
	ofstream out(outPath.c_str(), ios::out);
	string line;
	while (getline(in, line)) {
		//cerr << line << endl;
		vector<string> strs = split(line);
		string que = strs[1];
		if (query2id_.find(que) != query2id_.end()) {
			int qid = query2id_[que];
			vector<pair<int, double> > ret = rec_DQG(qid);
			out << que;
			for (int i = 0; i < ret.size(); ++i) {
				if (queries_[ret[i].first] == que) {
					continue;
				}
				out << '\t' << queries_[ret[i].first] << '\t' << ret[i].second;
			}
			out << endl;
		}
		else {
			out << que << endl;
		}
	}
	in.close();
	out.close();
	clock_t t2 = clock();
	cerr << "DQG recommendation takes " << (t2 - t1 + 0.0) / CLOCKS_PER_SEC << "seconds" << endl;
}

vector<pair<int, double> > DQG::rec_DQG(int qid)
{
	map<int, double> ink;
	ink.clear();
	ink[qid] = 1.0;
	return PPR_BCA(ink, EQFG_PPR_QUERY_ALPHA, 0.5, k_);
}

double DQG::query2docWeight(int qid, int did, double w, double beta)
{
	double lat = doc2cor_[did].first;
	double lon = doc2cor_[did].second;
	double dis = sqrt(pow(lat - Ulat, 2.0) + pow(lon - Ulon, 2.0));
	dis /= 180;
	if (dis > 1.0)
		dis = 1.0;
	return beta * w + (1 - beta) * (1.0 - dis);
}
double DQG::doc2queryWeight(int did, int qid, double w, double beta)
{
	double mindist = 999999;
	for (int i = 0; i < QNodes_[qid].toDocEdges_.size(); ++i) {
		int tempdid = QNodes_[qid].toDocEdges_[i].eid_;
		double lat = doc2cor_[tempdid].first;
		double lon = doc2cor_[tempdid].second;
		double dis = sqrt(pow(lat - Ulat, 2.0) + pow(lon - Ulon, 2.0));
		dis /= 180;
		if (dis > 1.0)
			dis = 1.0;
		if (dis < mindist)
			mindist = dis;
	}
	return beta * w + (1 - beta) * (1.0 - mindist);
}

vector<pair<int, double>> DQG::PPR_BCA(map<int, double> & initialInk, double alpha, double beta, int k)
{
	BoundHeap heap(QNodes_.size() + DNodes_.size());
	double activeInk = 0.0;
	map<int, double> result;

	vector<double> inkBuffer(QNodes_.size() + DNodes_.size(), 0.0);

	// initialize the heap
	for (map<int, double>::iterator i = initialInk.begin(); i != initialInk.end(); ++i) {
		heap.push(*i);
		activeInk += i->second;
	}

	while (heap.size() > 0) { //&& activeInk > PPR_EPS) {
		pair<int, double> topItem = heap.pop();
		if (topItem.second < PPR_IGNORE_INK) {
			break;
		}
		if (topItem.first >= QNodes_.size()) { // it is a doc
			int did = topItem.first - QNodes_.size();
			activeInk -= topItem.second * alpha;
			double distributedInk = (1.0 - alpha) * topItem.second;
			vector<EQFG_Edge> & edges = DNodes_[did].toQueryEdges_;
			vector<double> weights;
			double wSum = 0.0;
			for (int i = 0; i < edges.size(); ++i) {
				double spatialWeight = doc2queryWeight(did, edges[i].eid_, edges[i].w_, beta);
				wSum += spatialWeight;
				weights.push_back(spatialWeight);
			}
			for (int i = 0; i < edges.size(); ++i) {
				weights[i] /= wSum;
			}
			for (int i = 0; i < edges.size(); ++i) {
				double addInk = distributedInk * weights[i];

				// lazy update
				int tempQid = edges[i].eid_;
				inkBuffer[tempQid] += addInk;
				if (inkBuffer[tempQid] >= PPR_IGNORE_INK) {
					heap.push(make_pair(tempQid, inkBuffer[tempQid]));
					inkBuffer[tempQid] = 0.0;
				}
			}
		}
		else { // it is a query
			int qid = topItem.first;
			double increaseInk = topItem.second * alpha;
			activeInk -= topItem.second * alpha;
			if (result.find(qid) == result.end()) {
				result[qid] = 0.0;
			}
			result[qid] += increaseInk;
			activeInk -= increaseInk;
			double distributedInk = (1.0 - alpha) * topItem.second;
			vector<EQFG_Edge> & edges = QNodes_[qid].toDocEdges_;
			vector<double> weights;
			double wSum = 0.0;
			for (int i = 0; i < edges.size(); ++i) {
				double spatialWeight = query2docWeight(qid, edges[i].eid_, edges[i].w_, beta);
				wSum += spatialWeight;
				weights.push_back(spatialWeight);
			}
			for (int i = 0; i < edges.size(); ++i) {
				weights[i] /= wSum;
			}
			for (int i = 0; i < edges.size(); ++i) {
				double addInk = distributedInk * weights[i];
				
				// lazy update
				int tempdid = edges[i].eid_ + QNodes_.size();
				inkBuffer[tempdid] += addInk;
				if (inkBuffer[tempdid] >= PPR_IGNORE_INK) {
					heap.push(make_pair(tempdid, inkBuffer[tempdid]));
					inkBuffer[tempdid] = 0.0;
				}
				// the first QNodes.size() are query nodes, all the rest are doc nodes
				
				//heap.push(make_pair(edges[i].eid_ + QNodes_.size(), addInk)); // the first QNodes.size() are query nodes, all the rest are doc nodes
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