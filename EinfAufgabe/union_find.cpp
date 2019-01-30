#include <vector>
#include "union_find.hpp"

/*
Union-Find-Struktur auf der Menge 0,...,size-1.
Entspricht praktisch direkt den Algorithmen aus der Vorlesung
*/

UnionFind::UnionFind(index_t size): entries(size) {
	for (index_t i = 0;i<size;++i) {
		entry& e = entries[i];
		e.rank = 0;
		e.parent = i;
	}
}

UnionFind::index_t UnionFind::find(index_t start) {
	entry& currEntry = entries[start];
	if (currEntry.parent==start) {
		return start;
	}
	//In den meisten Aufrufen wird parent die Wurzel sein, durch diesen
	//"Doppelschritt" ist in diesem Fall keine Rekursion nötig.
	entry& parent = entries[currEntry.parent];
	if (parent.parent==currEntry.parent) {
		return currEntry.parent;
	}
	const index_t root = find(parent.parent);
	currEntry.parent = root;
	parent.parent = root;
	return root;
}

UnionFind::index_t UnionFind::merge(index_t startA, index_t startB) {
	return mergeRoots(find(startA), find(startB));
}

UnionFind::index_t UnionFind::mergeRoots(index_t rootA, index_t rootB) {
	if (rootA==rootB) {
		return rootA;
	}
	entry& entryA = entries[rootA];
	entry& entryB = entries[rootB];
	if (entryA.rank>entryB.rank) {
		entryB.parent = rootA;
		return rootA;
	} else {
		entryA.parent = rootB;
		if (entryA.rank==entryB.rank) {
			++entryA.rank;
		}
		return rootB;
	}
}
