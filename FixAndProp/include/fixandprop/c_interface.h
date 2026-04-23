#ifndef C_INTERFACE_H
#define C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Opaque handle to the transformer.
typedef void* TransformerHandle;

// Create and destroy a transformer.
TransformerHandle create_transformer(const char* name);
void destroy_transformer(TransformerHandle handle);

// Initialize the transformer with problem data.
// Here, for sparse rows we pass flattened data along with an array
// rowNNZ that tells how many nonzeros are in each row.
void transformer_init(
    TransformerHandle handle,
    int ncols,
    const double* xLb,
    const double* xUb,
    const char* xType,         // e.g. "BBB" for 3 binary variables
    int nrows,
    const int* rowNNZ,         // array of length nrows: number of nonzeros per row
    const int* sparseIndices,  // flattened indices for all constraints
    const double* sparseCoeffs,// flattened coefficients
    const char* senses,        // e.g. "EEE" for 3 equality constraints
    const double* rhs);

// Apply the propagation method.
// in, out, newLbs, newUbs are all arrays of length ncols.
bool transformer_apply_fix_and_prop(
    TransformerHandle handle,
    const double* in,
    double* out,
    double* newLbs,
    double* newUbs,
    int ncols,
    int ntofix,
    const int* toFix);

#ifdef __cplusplus
}
#endif

#endif // C_INTERFACE_H
