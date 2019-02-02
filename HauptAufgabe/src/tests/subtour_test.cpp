#include <tsp_instance.hpp>
#include <subtour_cut_gen.hpp>

int main() {
	Graph g;
	Graph::Node n0 = g.addNode();
	Graph::Node n1 = g.addNode();
	Graph::Node n2 = g.addNode();
	Graph::Node n3 = g.addNode();
	Graph::Node n4 = g.addNode();
	Graph::Node n5 = g.addNode();
	g.addEdge(n0, n1);
	g.addEdge(n0, n2);
	g.addEdge(n2, n1);
	g.addEdge(n4, n1);
	g.addEdge(n3, n2);
	g.addEdge(n3, n4);
	g.addEdge(n3, n5);
	g.addEdge(n4, n5);
	SubtourCutGen sct(g);
	std::vector<double> sol{1, 1, 1, 1, 0, 1, 1, 1};
	LinearProgram lp("test", LinearProgram::maximize);
	for (unsigned i = 0; i<8; ++i) {
		lp.addVariable(0, 0, 1);
	}
	std::cout << sct.validate(lp, sol) << std::endl;
	sol[4] = 1;
	std::cout << sct.validate(lp, sol) << std::endl;
}