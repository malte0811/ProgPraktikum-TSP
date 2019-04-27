#ifndef CONTRACTION_RULE_HPP
#define CONTRACTION_RULE_HPP

#include <utility>
#include <cstddef>
#include <vector>
#include <tsp_instance.hpp>
#include <functional>
#include <tsp_utils.hpp>


class ContractionRule {
public:
	using Contraction = std::vector<std::vector<Graph::Node>>;
	using Validator = std::function<Contraction(const Graph&, const std::vector<Graph::Node>&,
												const std::vector<Graph::Edge>&, const Graph::EdgeMap<double>&,
												const Graph::NodeMap<bool>&)>;

	explicit ContractionRule(ContractionRule::Validator rule);

	bool contractAll(Graph& g, lemon::GraphExtender<lemon::ListGraphBase>::NodeMap<bool>& used,
					 Graph::EdgeMap<double>& costs,
					 tsp_util::ContractionMapTSP& contrMap) const;

private:
	const Validator validate;

	bool contract(const Contraction& contr, Graph& g, Graph::EdgeMap<double>& costs,
				  tsp_util::ContractionMapTSP& map) const;
};

#endif