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

void printUsage(char * argv[])
{
	cerr << "Usage:" << endl;
	cerr << argv[0] << " -r -T queryLogDir locationDir alpha beta dist input output " << "\t for TQG recommendation" <<endl;
	cerr << argv[0] << " -r -D queryLogDir locationDir alpha beta input output " << "\t for DQG recommendation" << endl;
}

int main(int argc, char** argv)
{
	if (argv[1][1] == 'r') {
		if (argv[2][1] == 'T') {
			EQFG eqfg(argv[3]);
			eqfg.loadTerm(argv[3]);
			eqfg.loadLocation(argv[4]);
			eqfg.rec_TQG_fromfile(testPath, outPath, atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
		}
		else if (argv[2][1] == 'D') {
		
		}
	}
	else if(argv[1][1] == 'h'){
		printUsage(argv);
		return 0;
	}
	return 0;

	//loadInputPaths("paths.txt", argv);

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
	
	return 0;
}