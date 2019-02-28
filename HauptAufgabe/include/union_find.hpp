#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

#include <vector>
#include <cstdint>
#include <cstddef>

class UnionFind {
public:
	using rank_t = uint32_t;

	explicit UnionFind(size_t size);

	size_t find(size_t start);

	//oder union, aber das ist kein gültiger Methodenname
	size_t merge(size_t startA, size_t startB);

	size_t mergeRoots(size_t rootA, size_t rootB);

private:
	typedef struct {
		size_t parent;
		rank_t rank;
	} entry;
	std::vector<entry> entries;
};

#endif
