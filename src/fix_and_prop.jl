using SparseArrays
using MathOptInterface
const MOI = MathOptInterface

# Path to the shared library - platform-aware
const _lib_path = if Sys.isapple()
    "FixAndProp/build/libfp.dylib"
elseif Sys.iswindows()
    "FixAndProp/build/libfp.dll"
else
    "FixAndProp/build/libfp.so"
end

# Check if the library exists
const FIXANDPROP_AVAILABLE = isfile(_lib_path)
const lib = _lib_path
const TransformerHandle = Ptr{Nothing}

# Debug flag (shared with propagation.jl if both are loaded)
const DEBUG_FIX_AND_PROP = false

if !FIXANDPROP_AVAILABLE
    @warn "FixAndProp library not found at $lib. External propagation will be disabled. " *
          "To enable, compile FixAndProp: cd FixAndProp && mkdir -p build && cd build && cmake .. && make -j"
end

# ----------------------------------------------------------------
# Transformer Functions
# ----------------------------------------------------------------
mutable struct TransformerHandleWrapper
    handle::Ptr{Cvoid}
    function TransformerHandleWrapper(handle::Ptr{Cvoid})
        obj = new(handle)  # Use `new` to create the object without recursion.
        finalizer(obj) do o
            ccall((:destroy_transformer, lib), Cvoid, (Ptr{Cvoid},), o.handle)
        end
        return obj
    end
end

"""
Initialize a transformer with the problem data.
Requires FixAndProp library to be compiled.
"""
function init_transformer(xLb::Vector{Float64}, xUb::Vector{Float64}, xType::String,
    rowNNZ::Vector{Int}, sparseIndices::Vector{Int},
    sparseCoeffs::Vector{Float64}, senses::String, rhs::Vector{Float64})
    if !FIXANDPROP_AVAILABLE
        error("FixAndProp library not available. Set use_external_propagation=false or compile FixAndProp.")
    end
    ncols = Cint(length(xLb))  # number of variables
    nrows = Cint(length(rhs))  # number of constraints

    # Convert integer arrays to Cint.
    rowNNZ_c = Cint.(rowNNZ)
    sparseIndices_c = Cint.(sparseIndices)

    # Create the transformer handle using the transformer name "propround".
    handle = ccall((:create_transformer, lib), TransformerHandle, (Cstring,), "propround")

    # Initialize the transformer with the problem data.
    ccall((:transformer_init, lib), Cvoid,
        (TransformerHandle, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar},
            Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}),
        handle, ncols, pointer(xLb), pointer(xUb), pointer(xType), nrows,
        pointer(rowNNZ_c), pointer(sparseIndices_c), pointer(sparseCoeffs),
        pointer(senses), pointer(rhs)
    )
    return handle
end

"""
Apply the fix-and-propagate operation and return the new bounds and feasibility flag.
"""
function apply_fix_and_prop(handle::TransformerHandle, in_arr::Vector{Float64}, toFix::Vector{Cint})
    ncols = Cint(length(in_arr))
    ntoFix = Cint(length(toFix))

    # Allocate output arrays.
    out_arr = Vector{Cdouble}(undef, ncols)
    newLbs = Vector{Cdouble}(undef, ncols)
    newUbs = Vector{Cdouble}(undef, ncols)

    # Call the fix-and-propagate function from the shared library.
    feasible = ccall((:transformer_apply_fix_and_prop, lib), Bool,
        (TransformerHandle, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Cint, Cint, Ptr{Cint}),
        handle, pointer(in_arr), pointer(out_arr), pointer(newLbs), pointer(newUbs),
        ncols, ntoFix, pointer(toFix)
    )

    return newLbs, newUbs, feasible
end

"""
Destroy the transformer and free its resources.
"""
function destroy_transformer(handle::TransformerHandle)
    ccall((:destroy_transformer, lib), Cvoid, (TransformerHandle,), handle)
end

# ----------------------------------------------------------------
# Utility Functions for Constraint Processing
# ----------------------------------------------------------------


function build_sparse_linear_function_data(f::MOI.ScalarAffineFunction{Float64}, nvars::Int)
    c0 = f.constant
    b = spzeros(nvars)
    for term in f.terms
        b[term.variable.value] = term.coefficient
    end
    return b, c0
end


# ----------------------------------------------------------------
# Data Extraction for Propagation
# ----------------------------------------------------------------

"""
Extract and format the problem data needed for the transformer.
"""
function problem_data_for_propagator(global_variable_bounds, model)
    num_vars = num_variables(model)
    xLb = [global_variable_bounds.lower_bounds[i] for i in 1:num_vars]
    xUb = [global_variable_bounds.upper_bounds[i] for i in 1:num_vars]
    integer_variables = get_integer_variables_indices(model)
    binary_variables = get_binary_variables_indices(model)
    xType = ""
    for i in 1:num_vars
        if i in binary_variables
            xType *= "B"
        elseif i in integer_variables
            xType *= "I"
        else
            xType *= "C"
        end
    end

    rowNNZ = Vector{Int}()
    sparseIndices = Vector{Int}()
    sparseCoeffs = Vector{Float64}()
    all_senses = ""
    all_rhs = Vector{Float64}()

    # Process equality constraints.
    eq_ids = MOI.get(model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}())
    for id in eq_ids
        f = MOI.get(model, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(model, MOI.ConstraintSet(), id)
        b, c0 = build_sparse_linear_function_data(f, num_vars)
        rhs, _ = get_rhs_and_type(cons_set)
        push!(rowNNZ, length(findall(!iszero, b)))
        append!(sparseIndices, findall(!iszero, b))
        append!(sparseCoeffs, b.nzval)
        all_senses *= "E"
        push!(all_rhs, rhs - c0)
    end

    # Process less-than constraints.
    lt_ids = MOI.get(model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}())
    for id in lt_ids
        f = MOI.get(model, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(model, MOI.ConstraintSet(), id)
        rhs, _ = get_rhs_and_type(cons_set)
        b, c0 = build_sparse_linear_function_data(f, num_vars)
        push!(rowNNZ, length(findall(!iszero, b)))
        append!(sparseIndices, findall(!iszero, b))
        append!(sparseCoeffs, b.nzval)
        all_senses *= "L"
        push!(all_rhs, rhs - c0)
    end

    # Process greater-than constraints.
    gt_ids = MOI.get(model, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}}())
    for id in gt_ids
        f = MOI.get(model, MOI.ConstraintFunction(), id)
        cons_set = MOI.get(model, MOI.ConstraintSet(), id)
        rhs, _ = get_rhs_and_type(cons_set)
        b, c0 = build_sparse_linear_function_data(f, num_vars)
        push!(rowNNZ, length(findall(!iszero, b)))
        append!(sparseIndices, findall(!iszero, b))
        append!(sparseCoeffs, b.nzval)
        all_senses *= "G"
        push!(all_rhs, rhs - c0)
    end

    # Adjust indices to be zero-based.
    sparseIndices = sparseIndices .- 1
    return xLb, xUb, xType, rowNNZ, sparseIndices, sparseCoeffs, all_senses, all_rhs
end

function initialize_propagation(global_variable_bounds, model)
    # Extract formatted data using individual components.
    xLb, xUb, xType, rowNNZ, sparseIndices, sparseCoeffs, senses, rhs =
        problem_data_for_propagator(global_variable_bounds, model)
    # Initialize transformer.
    handle = init_transformer(xLb, xUb, xType, rowNNZ, sparseIndices, sparseCoeffs, senses, rhs)
    return handle
end

# ----------------------------------------------------------------
# Propagation Wrapper Function
# ----------------------------------------------------------------

"""
Perform propagation using the transformer.
Extracts problem data, calls the transformer, and updates bounds.
"""
function apply_prop_external(model, global_variable_bounds, fix_list, x_ref)
    propagation_handle = initialize_propagation(global_variable_bounds, model)

    wrapped_handle = TransformerHandleWrapper(propagation_handle)

    # Adjust fix_list indices to zero-based.
    fix_list_c = Cint.(fix_list) .- Cint(1)

    # Apply fix-and-propagate.
    newLbs, newUbs, feasible = apply_fix_and_prop(wrapped_handle.handle, x_ref, fix_list_c)

    if DEBUG_FIX_AND_PROP
        println("Propagation result feasible: $feasible")
    end

    # Update the global bounds.
    new_bounds = copy(global_variable_bounds)
    for i in 1:length(x_ref)
        new_bounds.lower_bounds[i] = newLbs[i]
        new_bounds.upper_bounds[i] = newUbs[i]
    end

    num_reduced_domains = 0
    for i in 1:length(x_ref)
        if i in fix_list
            continue
        else
            if new_bounds.lower_bounds[i] > global_variable_bounds.lower_bounds[i] + 1e-5 ||
               new_bounds.upper_bounds[i] < global_variable_bounds.upper_bounds[i] - 1e-5
                num_reduced_domains += 1
            end
        end
    end
    if DEBUG_FIX_AND_PROP
        println("Number of reduced domains: $num_reduced_domains")
    end

    # Verification: Ensure new bounds are within original limits.
    if DEBUG_FIX_AND_PROP && feasible
        for i in 1:length(x_ref)
            if !(global_variable_bounds.lower_bounds[i] - 1e-5 <= new_bounds.lower_bounds[i] <=
                 global_variable_bounds.upper_bounds[i] + 1e-5)
                @assert false
            end
            if !(global_variable_bounds.lower_bounds[i] - 1e-5 <= new_bounds.upper_bounds[i] <=
                 global_variable_bounds.upper_bounds[i] + 1e-5)
                @assert false
            end
            # Check that variables in the min cover are fixed.
            if i in fix_list && new_bounds.lower_bounds[i] != new_bounds.upper_bounds[i]
                @assert false
            end
        end
    end

    return new_bounds.lower_bounds, new_bounds.upper_bounds, !feasible, num_reduced_domains
end
