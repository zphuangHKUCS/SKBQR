#ifndef SRC_ALGORITHM_BOUNDHEAP_H_
#define SRC_ALGORITHM_BOUNDHEAP_H_
#include <cstdio>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <map>
#include <boost/heap/fibonacci_heap.hpp>

using namespace std;
class ComparatorScoreUpperBoundLess{
public:
	inline bool operator() (const pair<int, double> &a, const pair<int, double> &b) const
	{
		return a.second < b.second;
	}
};
typedef boost::heap::fibonacci_heap<pair<int, double>, boost::heap::compare<ComparatorScoreUpperBoundLess> > Heap;


class BoundHeap {
private:
	
	Heap boundHeap;
	unordered_map<int, Heap::handle_type> handleMap;
	int k;
public:
	BoundHeap(int kk): k(kk){}

	void push(pair<int, double> cand);
    pair<int, double> pop()
    {
        int tid = boundHeap.top().first;
        double value = boundHeap.top().second;
        boundHeap.pop();
        handleMap.erase(tid);
        return make_pair(tid, value);
    }
	double kthBound() const
	{
		if (this->boundHeap.size() < k)
			return std::numeric_limits<double>::max();
		else
			return this->boundHeap.top().second;
	}
    int size()
    {
        return handleMap.size();
    }
    bool find(int tid)
    {
        return handleMap.find(tid) != handleMap.end();
    }
	void print () const;
};



#endif /* SRC_ALGORITHM_BOUNDHEAP_H_ */
