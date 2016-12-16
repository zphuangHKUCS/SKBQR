#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <cstdio>
#include "EQFG.h"


using namespace std;

string qlPath, testPath, outPath;

void loadInputPaths(string pathPath, char** argv)
{
	ifstream pathin(pathPath.c_str(), ios::in);
	string line;
	getline(pathin, line);
	if (line == "WIN") {
		getline(pathin, qlPath);
		getline(pathin, testPath);
		getline(pathin, outPath);
	}
	else {
		getline(pathin, qlPath);
		getline(pathin, testPath);
		getline(pathin, outPath);
	}
	
}


int main(int argc, char** argv)
{
	loadInputPaths("paths.txt", argv);
	EQFG eqfg(qlPath);
	eqfg.rec_EQFG_fromfile(testPath, outPath);
	//vector<pair<int, double>> temp = eqfg.rec_EQFG(1079514);
	return 0;
}