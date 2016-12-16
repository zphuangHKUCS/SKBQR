#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <cstdio>
#include "EQFG.h"


using namespace std;

void loadInputPaths(string pathPath, string & qlpath) 
{
	ifstream pathin = ifstream(pathPath.c_str(), ios::in);
	getline(pathin, qlpath);
}


int main()
{
	string qlPath;
	loadInputPaths("paths.txt", qlPath);
	EQFG eqfg(qlPath);
	vector<pair<int, double>> result = eqfg.rec_QFG(2);
	return 0;
}