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
	ifstream pathin = ifstream(pathPath.c_str(), ios::in);
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
	eqfg.rec_QFG_fromfile(testPath, outPath);
	return 0;
}