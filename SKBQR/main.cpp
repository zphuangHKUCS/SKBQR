#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <cstdio>
#include "EQFG.h"

using namespace std;

string qlPath, testPath, locDir, outPath;

void loadInputPaths(string pathPath, char** argv)
{
	ifstream pathin(pathPath.c_str(), ios::in);
	string line;
	getline(pathin, line);
	if (line == "WIN") {
		getline(pathin, qlPath);
		getline(pathin, testPath);
		getline(pathin, locDir);
		getline(pathin, outPath);
	}
	else {
		getline(pathin, qlPath);
		getline(pathin, testPath);
		getline(pathin, locDir);
		getline(pathin, outPath);
	}
	
}


int main(int argc, char** argv)
{

	loadInputPaths("paths.txt", argv);

	/*    // EQFG
	EQFG eqfg(qlPath);
	eqfg.loadEntity(qlPath);
	eqfg.loadLocation(locDir);
	eqfg.rec_EQFG_fromfile(testPath, outPath);
	*/

	/*   //  DQG
	DQG dqg(qlPath);
	dqg.loadLocation(locDir);
	dqg.rec_DQG_fromfile(testPath, outPath);
	*/
	
	     // TQG
	EQFG eqfg(qlPath);
	eqfg.loadTerm(qlPath);
	eqfg.rec_TQG_fromfile(testPath, outPath);
	return 0;
}