#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <cstddef>
#include <vector>

using std::size_t;

class UnionFind {
public:
	using rank_t = size_t;

	explicit UnionFind(size_t size);

	size_t find(size_t start);

	//oder union, aber das ist kein gültiger Methodenname
	size_t merge(size_t startA, size_t startB);

	size_t mergeRoots(size_t rootA, size_t rootB);

private:
	struct entry {
		size_t parent;
		rank_t rank;
	};
	std::vector<entry> entries;
	size_t numSets;
};

#endif
