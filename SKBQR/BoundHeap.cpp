#include "BoundHeap.h"

#include <iostream>
using namespace std;

void BoundHeap::push(pair<int, double>  cand)
{
	int Id = cand.first;
	double scoreUpperBound = cand.second;
	if (handleMap.find(Id) != handleMap.end())
	{
		Heap::handle_type handle = handleMap[Id];
        double newValue = scoreUpperBound + (*handle).second;
		//boundHeap.update(handle, std::make_pair(Id, scoreUpperBound));
        //cerr << "old : " <<  (*handle).second << endl;
        boundHeap.update(handle, std::make_pair(Id, newValue));
        //cerr << "diff: " << scoreUpperBound << endl;
        //cerr << "new : " <<  (*handle).second << endl;
        //cerr << "top : " << boundHeap.top().second << endl;
	}
	else
	{
		if (boundHeap.size() < k)
		{
			Heap::handle_type handle = boundHeap.push(std::make_pair(Id, scoreUpperBound));
			handleMap[Id] = handle;
		}
		else if (scoreUpperBound < this->kthBound())
		{
			handleMap.erase(boundHeap.top().first);
			boundHeap.pop();
			handleMap[Id] = boundHeap.push(std::make_pair(Id, scoreUpperBound));
		}
	}
}

void BoundHeap::print() const
{
	for (Heap::const_iterator iter = this->boundHeap.begin(); iter != this->boundHeap.end(); ++iter)
	{
		int trajId = iter->first;
		double scoreUpperBound = iter->second;
		std::cout << trajId << ": " << scoreUpperBound << std::endl;
	}
}


