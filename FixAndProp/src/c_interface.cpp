#include "fixandprop/c_interface.h"
#include "fixandprop/transformers.h"  // Provides PropagatorRounding etc.
#include "fixandprop/fp_interface.h"
#include <vector>
#include <string>

using namespace dominiqs;

extern "C" {

// For simplicity we wrap our transformer in a simple struct.
struct TransformerWrapper {
    // We assume that "propround" creates an object of type PropagatorRounding.
    PropagatorRounding* transformer;
};

void print_vector_double(const std::vector<double>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void print_vector_int(const std::vector<int>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void print_vector_char(const std::vector<char>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void print_vector_string(const std::vector<std::string>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

TransformerHandle create_transformer(const char* name) {
    TransformerWrapper* wrapper = new TransformerWrapper;
    // Use the factory to create the transformer.
    std::cout << "Creating transformer with name: " << name << std::endl;
    wrapper->transformer = static_cast<PropagatorRounding*>(
        TransformersFactory::getInstance().create(std::string(name))
    );
    return static_cast<TransformerHandle>(wrapper);
}

void destroy_transformer(TransformerHandle handle) {
    TransformerWrapper* wrapper = static_cast<TransformerWrapper*>(handle);
    delete wrapper->transformer;
    delete wrapper;
}

void transformer_init(
    TransformerHandle handle,
    int ncols,
    const double* xLb,
    const double* xUb,
    const char* xType,         // e.g. "BBB"
    int nrows,
    const int* rowNNZ,         // length nrows
    const int* sparseIndices,  // flattened indices
    const double* sparseCoeffs,// flattened coefficients
    const char* senses,        // e.g. "EEE"
    const double* rhs)
{
    TransformerWrapper* wrapper = static_cast<TransformerWrapper*>(handle);
    // Convert C arrays to std::vector types.
    std::vector<double> vec_xLb(xLb, xLb + ncols);
    std::vector<double> vec_xUb(xUb, xUb + ncols);
    std::vector<char> vec_xType(xType, xType + ncols);
    
    // Build the sparseRows vector.
    std::vector<SparseVector> sparseRows;
    int indexOffset = 0;
    for (int i = 0; i < nrows; i++) {
        SparseVector sv;
        int nnz = rowNNZ[i];
        std::vector<int> idx(nnz);
        std::vector<double> coef(nnz);
        for (int j = 0; j < nnz; j++) {
            idx[j] = sparseIndices[indexOffset];
            coef[j] = sparseCoeffs[indexOffset];
            indexOffset++;
        }
       sv.copy(idx.data(), coef.data(), nnz);
        sparseRows.push_back(sv);
    }
    
    std::vector<char> vec_senses(senses, senses + nrows);
    std::vector<double> vec_rhs(rhs, rhs + nrows);
    
    wrapper->transformer->readConfig();
    
    // Call the C++ init method.
    wrapper->transformer->init(ncols, vec_xLb, vec_xUb, vec_xType, nrows, sparseRows, vec_senses, vec_rhs);
}

bool transformer_apply_fix_and_prop(
    TransformerHandle handle,
    const double* in,
    double* out,
    double* newLbs,
    double* newUbs,
    int ncols,
    int ntofix,
    const int* toFix)
{
    TransformerWrapper* wrapper = static_cast<TransformerWrapper*>(handle);
    
    // Convert input array to vector.
    std::vector<double> vec_in(in, in + ncols);
    std::vector<double> vec_out(ncols, 0.0);
    std::vector<double> vec_newLbs(ncols, 0.0);
    std::vector<double> vec_newUbs(ncols, 0.0);
    std::vector<int> vec_toFix(toFix, toFix + ntofix);

    // print_vector_double(vec_in);
    // print_vector_double(vec_out);
    // print_vector_double(vec_newLbs);
    // print_vector_double(vec_newUbs);
    // print_vector_int(vec_toFix);
    // Call the C++ method.
    bool feasible = wrapper->transformer->apply_fix_and_prop(vec_in, vec_out, vec_newLbs, vec_newUbs, vec_toFix);
    
    // Copy results back to the provided arrays.
    for (int i = 0; i < ncols; i++) {
        out[i] = vec_out[i];
        newLbs[i] = vec_newLbs[i];
        newUbs[i] = vec_newUbs[i];
    }
    return feasible;
}

} // extern "C"
