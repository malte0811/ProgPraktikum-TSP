#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <vector>
#include <cstdint>

class UnionFind {
public:
	using index_t = uint32_t;
	using rank_t = uint32_t;
	explicit UnionFind(index_t size);
	index_t find(index_t start);
	//oder union, aber das ist kein gültiger Methodenname
	index_t merge(index_t startA, index_t startB);
	index_t mergeRoots(index_t rootA, index_t rootB);
private:
	typedef struct {
		index_t parent;
		rank_t rank;
	} entry;
	std::vector<entry> entries;
};
#endif
