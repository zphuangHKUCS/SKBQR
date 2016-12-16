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
	cerr << line << endl;
	if (line == "WIN") {
		getline(pathin, qlPath);
		getline(pathin, testPath);
		getline(pathin, outPath);
		cerr << qlPath << endl;
		cerr << testPath << endl;
		cerr << outPath << endl;
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