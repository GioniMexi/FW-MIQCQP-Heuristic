/**
 * @file transformers.cpp
 * @brief Solution Transformers (a.k.a. rounders) Source
 */

#include <map>
#include <list>
#include <iterator>
#include <set>
#include <signal.h>
#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <utils/floats.h>
#include <utils/fileconfig.h>
#include <utils/consolelog.h>

#include "fixandprop/transformers.h"

using namespace dominiqs;


// macro type savers
#define READ_FROM_CONFIG( what, defValue ) what = gConfig().get("fp."#what, defValue)
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)
#define LOG_CONFIG( what ) LOG_ITEM("fp."#what, what)

// Rounders

static bool DEF_RANDOMIZED_ROUNDING = true;
static bool DEF_LOG_DETAILS = true;
static uint64_t DEF_SEED = 0;

PropagatorRounding::PropagatorRounding() : randomizedRounding(DEF_RANDOMIZED_ROUNDING), logDetails(DEF_LOG_DETAILS)
{
}

void PropagatorRounding::readConfig()
{
	READ_FROM_CONFIG( randomizedRounding, DEF_RANDOMIZED_ROUNDING );
	READ_FROM_CONFIG( logDetails, DEF_LOG_DETAILS );
	consoleInfo("[config rounder]");
	LOG_CONFIG( randomizedRounding );
	LOG_CONFIG( logDetails );
	uint64_t seed = gConfig().get<uint64_t>("seed", DEF_SEED);
	roundGen.setSeed(seed);
	roundGen.warmUp();
	filterConstraints = gConfig().get("fp.filterConstraints", true);
	LOG_ITEM("fp.filterConstraints", filterConstraints);
}

void PropagatorRounding::init(
    int ncols,
    const std::vector<double>& xLb,
    const std::vector<double>& xUb,
    const std::vector<char>& xType,
    int nrows,
    const std::vector<dominiqs::SparseVector>& sparseRows,
    const std::vector<char>& senses,
    const std::vector<double>& rhs)
{
	domain = std::make_shared<Domain>();

    // Add variables to the domain using the provided variable names.
    for (int j = 0; j < ncols; j++) {
        domain->pushVar("var" + std::to_string(j), xType[j], xLb[j], xUb[j]);
    }

    // Connect domain to engine
    prop.setDomain(domain);

    // Generate propagators with analyzers.
    std::list<std::string> fNames;
    PropagatorFactories::getInstance().getIDs(std::back_insert_iterator<std::list<std::string>>(fNames));
    for (const std::string& name : fNames) {
        PropagatorFactoryPtr fact(PropagatorFactories::getInstance().create(name));
        factories[fact->getPriority()] = fact;
    }

    int filteredOut = 0;
    for (int i = 0; i < nrows; i++) {
        auto itr = factories.begin();
        auto end = factories.end();
        ConstraintPtr c = std::make_shared<Constraint>();

        // Fill the constraint using the provided vectors.
        c->row   = sparseRows[i];
        c->sense = senses[i];
        c->rhs   = rhs[i];

        DOMINIQS_ASSERT(c->sense == 'E' || c->sense == 'L' || c->sense == 'G');

        // Optionally filter constraints.
        if (filterConstraints) {
            const int* idx = c->row.idx();
            const double* coef = c->row.coef();
            unsigned int size = c->row.size();
            bool allCont = true;
            double largest = std::numeric_limits<double>::min();
            double smallest = std::numeric_limits<double>::max();
            for (unsigned int k = 0; k < size; k++) {
                if (!domain->isVarFixed(idx[k]) && (domain->varType(idx[k]) != 'C')) {
                    allCont = false;
                    break;
                }
                double tmp = fabs(coef[k]);
                largest = std::max(largest, tmp);
                smallest = std::min(smallest, tmp);
            }
            double dynamism = (largest / smallest);
            if ((allCont && greaterThan(dynamism, 10.0)) || greaterThan(dynamism, 1000.0)) {
                filteredOut++;
                continue;
            }
        }

        // Try to analyze the constraint with each available propagator.
        while (itr != end) {
            PropagatorPtr p = itr->second->analyze(*(domain.get()), c.get());
            if (p) {
                prop.pushPropagator(p);
                break;
            }
            ++itr;
        }
    }

    // Log propagator stats.
    consoleInfo("[propagator stats]");
    for (const auto& kv : factories) {
        consoleLog("{}: {}", kv.second->getName(), kv.second->created());
    }
    consoleLog("#filtered out: {}\n", filteredOut);

    // Dump state.
    state = prop.getStateMgr();
    state->dump();
}

void PropagatorRounding::ignoreGeneralIntegers(bool flag)
{
}

bool PropagatorRounding::apply_fix_and_prop(const std::vector<double>& in, std::vector<double>& out,
							   std::vector<double>& newLbs, std::vector<double>& newUbs,
							   const std::vector<int>& toFix)
{

	state->restore();

	double t = getRoundingThreshold(0.5, roundGen);

	// main loop
	int next;
	int nToFix = toFix.size();
	int i = 0;

	if (nToFix == 0)
	{
		bool stillFeasible = prop.propagate();
		if (!stillFeasible) return false;
		// update bounds for all variables
		for (int j = 0; j < in.size(); j++)
		{
			newLbs[j] = domain->varLb(j);
			newUbs[j] = domain->varUb(j);
		}
		return true;
	}

	while(i < nToFix) 
	{
		next = toFix[i];
 		// standard rounding
		double newval = in[next];
		if (domain->varType(next) == 'B') doRound(in[next], newval, t);
		else if(domain->varType(next) == 'I')
		{
			// general integer variable: here we take into account also the tightened domain
			// in order to have a more clever rounding!
			if (lessEqualThan(in[next], domain->varLb(next))) newval = domain->varLb(next);
			else if (greaterEqualThan(in[next], domain->varUb(next))) newval = domain->varUb(next);
			else doRound(in[next], newval, t);
		}
		else
		{
			// continuous variables: here we also take into account the tightened domain
			if (lessEqualThan(in[next], domain->varLb(next))) newval = domain->varLb(next);
			else if (greaterEqualThan(in[next], domain->varUb(next))) newval = domain->varUb(next);
			else newval = in[next];
		}
		// propagate
		bool stillFeasible = prop.propagate(next, in[next]);
		// consoleLog("Still feasible after fixing: {}", stillFeasible);
		if (!stillFeasible) return false;
		DOMINIQS_ASSERT( domain->isVarFixed(next) );
		i++;
	}
	// update bounds for all variables
	for (int j = 0; j < in.size(); j++)
	{
		newLbs[j] = domain->varLb(j);
		newUbs[j] = domain->varUb(j);
	}
	return true;
}

void PropagatorRounding::clear()
{
	// clear
	// delete state;
	prop.clear();
	factories.clear();
}

// auto registration
// please register your class here with an appropriate name

class TR_FACTORY_RECORDER
{
public:
	TR_FACTORY_RECORDER()
	{
		std::cout << "Registering SolutionTransformers...";
		TransformersFactory::getInstance().registerClass<PropagatorRounding>("propround");
		std::cout << "done" << std::endl;
	}
};

TR_FACTORY_RECORDER my_tr_factory_recorder;
