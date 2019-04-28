#include <union_find.hpp>
#include <vector>

/*
Union-Find-Struktur auf der Menge 0,...,size-1.
Entspricht praktisch direkt den Algorithmen aus der Vorlesung
*/
UnionFind::UnionFind(size_t size) : entries(size), numSets(size) {
	for (size_t i = 0; i < size; ++i) {
		entry& e = entries[i];
		e.rank = 0;
		e.parent = i;
	}
}

size_t UnionFind::find(size_t start) {
	if (numSets <= 1) {
		return 0;
	}
	entry& currEntry = entries[start];
	if (currEntry.parent == start) {
		return start;
	}
	//In den meisten Aufrufen wird parent die Wurzel sein, durch diesen
	//"Doppelschritt" ist in diesem Fall keine Rekursion nötig.
	entry& parent = entries[currEntry.parent];
	if (parent.parent == currEntry.parent) {
		return currEntry.parent;
	}
	const size_t root = find(parent.parent);
	currEntry.parent = root;
	parent.parent = root;
	return root;
}

size_t UnionFind::merge(size_t startA, size_t startB) {
	return mergeRoots(find(startA), find(startB));
}

size_t UnionFind::mergeRoots(size_t rootA, size_t rootB) {
	if (rootA == rootB) {
		return rootA;
	}
	--numSets;
	entry& entryA = entries[rootA];
	entry& entryB = entries[rootB];
	if (entryA.rank > entryB.rank) {
		entryB.parent = rootA;
		return rootA;
	} else {
		entryA.parent = rootB;
		if (entryA.rank == entryB.rank) {
			++entryA.rank;
		}
		return rootB;
	}
}
