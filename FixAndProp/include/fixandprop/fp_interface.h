/**
 * @file fp_interface.h
 * @brief Solution Transformers (a.k.a. rounders) Interface
 */

#ifndef FP_INTERFACE_H
#define FP_INTERFACE_H

#include <vector>
#include <memory>

#include <utils/singleton.h>
#include <utils/factory.h>
#include <utils/fileconfig.h>

namespace dominiqs {

/**
 * Solution Transformer interface
 * Base class for frac->int (i.e. rounding) transformations
 */

class SolutionTransformer
{
public:
	virtual ~SolutionTransformer() {}
	virtual void readConfig() {}
	/**
	 * Read needed information (if any) about the problem (@param pinfo)
	 */
	virtual void init(
		int ncols,
		const std::vector<double>& xLb,
		const std::vector<double>& xUb,
		const std::vector<char>& xType,
		const std::vector<std::string>& xNames,
		int nrows,
		const std::vector<SparseVector>& sparseRows,
		const std::vector<char>& senses,
		const std::vector<double>& rhs) {}
	virtual void ignoreGeneralIntegers(bool flag) {}
	/**
	 * Trasform the vector given as input @param in and store the result in @param out
	 */
	virtual bool apply_fix_and_prop(const std::vector<double>& in, std::vector<double>& out, 
		                            std::vector<double>& newLbs, std::vector<double>& newUbs, 
									const std::vector<int>& toFix) = 0;
	/**
	 *
	 */
	virtual void clear() {}
};

typedef std::shared_ptr<SolutionTransformer> SolutionTransformerPtr;

typedef SingletonHolder<Factory<SolutionTransformer, std::string>> TransformersFactory;

#ifdef LIBFP_STATIC
// manual registration
void registerTransformers();
#endif // LIBFP_STATIC

} // namespace dominiqs

#endif /* FP_INTERFACE_H */
