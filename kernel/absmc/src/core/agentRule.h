#ifndef ABSMC_AGENT_RULE_H
#define ABSMC_AGENT_RULE_H

#include <vector>
#include <cassert>
#include <memory>

namespace absmc {

/// Base class for classes defining agent rules.
template<class Agent>
class AgentRule {
public:
    typedef Agent agent_t;
    typedef typename Agent::base_t agent_base_t;

    virtual ~AgentRule() { }

    /// Apply the rule to the given agent.
    virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const = 0;
};

/// Composite agent rule, chaining several agent rules.
/// This provides a mechanism to compose an agent rule set as a sequence of several simple rules.
/// Rules are applied in the same order as they are added to the composite rule object.
template<class Agent>
class CompositeRule : public AgentRule<Agent> {
public:
    typedef Agent agent_t;
    typedef typename Agent::base_t agent_base_t;
    typedef AgentRule<Agent> rule_t;

    void add(std::shared_ptr<rule_t> rule) {
        assert(rule);
        ruleVector.push_back(rule);
    }
    void clear() {
        ruleVector.clear();
    }
    virtual void apply(size_t iAgent, agent_t & agent, std::vector<agent_base_t*> & agentsCreated, std::vector<agent_base_t*> & agentsDeleted) const {
        for (size_t i=0; i<ruleVector.size(); i++) {
            (ruleVector[i])->apply(iAgent, agent, agentsCreated, agentsDeleted);
        }
    }

private:
    std::vector< std::shared_ptr<rule_t> > ruleVector;
};

} // end namespace absmc

#endif
