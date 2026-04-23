/**
 * @file main.cpp
 * @brief Main App using the C API wrapper for the transformer.
 */

#include <iostream>
#include "fixandprop/c_interface.h"  // The C API header
#include <utils/consolelog.h>

using namespace std;

int main (int argc, char const *argv[])
{
	// Create a transformer handle via the C API.
	TransformerHandle handle = create_transformer("propround");
	
	// Define problem data:
	// 3 binary variables: x1, x2, x3.
	int ncols = 3;
	double xLb[] = {0.0, 0.0, 0.0};
	double xUb[] = {1.0, 1.0, 1.0};
	// xType is a string representing each variable type. For 3 binary vars, "BBB".
	char xType[] = "BBB";
	// Array of variable names.
	const char* xNames[] = {"x1", "x2", "x3"};

	int ntoFix = 2;  // Number of variables to fix.
	int toFix[] = {0, 1};  // Fix the first two variables.
	// Define 3 constraints:
	// Constraint 0: x1 + x2 = 1   (indices 0 and 1)
	// Constraint 1: x1 + x3 = 1   (indices 0 and 2)
	// Constraint 2: x2 + x3 = 0   (indices 1 and 2)
	int nrows = 3;
	// For each constraint, specify the number of nonzeros.
	int rowNNZ[] = {2, 2, 2};
	// Flattened indices: for constraint 0: 0,1; for 1: 0,2; for 2: 1,2
	int sparseIndices[] = {0, 1, 0, 2, 1, 2};
	// Flattened coefficients (all 1.0).
	double sparseCoeffs[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	// Constraint senses: "EEE" means all equality constraints.
	char senses[] = "EEE";
	// Right-hand side values.
	double rhs[] = {1.0, 1.0, 0.0};

	// Initialize the transformer with the problem data.
	transformer_init(handle, ncols, xLb, xUb, xType, nrows, rowNNZ, sparseIndices, sparseCoeffs,
					senses, rhs);

	// Prepare input vector and arrays for the result.
	double in[]  = {0.6, 0.2, 0.55};  // fractional solution input
	double out[3] = {0.0, 0.0, 0.0};    // solution output
	double newLbs[3] = {0.0, 0.0, 0.0};
	double newUbs[3] = {0.0, 0.0, 0.0};

	// Call the transformer’s apply_fix_and_prop method via the C API.
	transformer_apply_fix_and_prop(handle, in, out, newLbs, newUbs, ncols, ntoFix, toFix);

	for (int j = 0; j < ncols; j++)
	{
		consoleLog("Variable {} has bounds [{}, {}]", j, newLbs[j], newUbs[j]);
	}

	int nNoFix = 0;  // Number of variables to fix.
	int NoFix[] = {};  // Fix the first two variables.

	// Call the transformer’s apply_fix_and_prop method via the C API.
		transformer_apply_fix_and_prop(handle, in, out, newLbs, newUbs, ncols, nNoFix, NoFix);

	for (int j = 0; j < ncols; j++)
	{
		consoleLog("Variable {} has bounds [{}, {}]", j, newLbs[j], newUbs[j]);
	}

	// When finished, destroy the transformer.
	destroy_transformer(handle);

	return 0;
}
