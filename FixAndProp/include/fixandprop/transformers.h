/**
 * @file transformers.h
 * @brief Solution Transformers (a.k.a. rounders) Header
 */

#ifndef TRANSFORMERS_H
#define TRANSFORMERS_H

#include <utils/randgen.h>

#include <propagator/domain.h>
#include <propagator/prop_engine.h>

#include "fp_interface.h"

/**
 * Rounding helpers
 */

inline void doRound(const double& in, double& out, const double& thr) { out = floor(in + thr); }

inline double getRoundingThreshold(bool isRandom, dominiqs::RandGen& gen)
{
	if (isRandom)
	{
		//double t = gen.getFloat();
		//if (t <= 0.5) return 2 * t * (1 - t);
		//else return 2 * t * (t - 1) + 1;
		return gen.getFloat();
	}
	else return 0.5;
}

/**
 * rounding + constraint propagation
 */

class PropagatorRounding : public dominiqs::SolutionTransformer
{
public:
	PropagatorRounding();
	~PropagatorRounding() { clear(); }
	void readConfig();
	void init(
		int ncols,
		const std::vector<double>& xLb,
		const std::vector<double>& xUb,
		const std::vector<char>& xType,
		int nrows,
		const std::vector<dominiqs::SparseVector>& sparseRows,
		const std::vector<char>& senses,
		const std::vector<double>& rhs);
	void ignoreGeneralIntegers(bool flag);
	bool apply_fix_and_prop(const std::vector<double>& in, std::vector<double>& out,
		       				std::vector<double>& newLbs, std::vector<double>& newUbs, const std::vector<int>& toFix);
	void clear();
protected:
	dominiqs::RandGen roundGen;
	bool randomizedRounding;
	bool logDetails;
	// data
	DomainPtr domain;
	StatePtr state;
	PropagationEngine prop;
	std::map<int, PropagatorFactoryPtr> factories;
	bool filterConstraints;
};

#endif /* TRANSFORMERS_H */
