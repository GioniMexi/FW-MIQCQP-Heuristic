#!/bin/bash
#
# Usage: run_instance.sh <binary_path> <instance_result_dir> <instance> <parameter_file>
#
# Environment variables:
#   EXECUTION_MODE - "cluster" (uses srun) or "local" (direct execution)
#   NUM_THREADS    - Number of Julia threads (default: 8)

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <binary_path> <instance_result_dir> <instance> <parameter_file>"
    exit 1
fi

binary_path="$1"
instance_result_dir="$2"
instance="$3"
parameter_file="$4"

# Configuration with defaults
EXECUTION_MODE="${EXECUTION_MODE:-cluster}"
NUM_THREADS="${NUM_THREADS:-8}"

# Unset GUROBI_HOME to use the bundled Gurobi from the binary
# This ensures consistent Gurobi version regardless of system installation
unset GUROBI_HOME

# Validate inputs
if [ ! -x "$binary_path" ]; then
    echo "Binary not found or not executable: $binary_path"
    exit 1
fi

if [ ! -f "$instance" ]; then
    echo "Instance file not found: $instance"
    exit 1
fi

if [ ! -f "$parameter_file" ]; then
    echo "Parameter file not found: $parameter_file"
    exit 1
fi

# Create the results directory if it doesn't exist.
mkdir -p "$instance_result_dir"

# Derive the base name (without .* extension) and set the log file path.
base_name=$(basename "$instance")
log_file="${instance_result_dir}/${base_name%.*}.log"

echo "Processing instance: $instance"
echo "Result directory: $instance_result_dir"
echo "Binary log file: $log_file"
echo "Execution mode: $EXECUTION_MODE"

# Record the start time (in nanoseconds)
start_time=$(date +%s%N)

# Build command arguments
CMD_ARGS="--file_path=$instance --result_path=$instance_result_dir --parameter_file=$parameter_file"
JULIA_ARGS="--julia-args --threads=$NUM_THREADS --check-bounds=no"

"$binary_path" $CMD_ARGS $JULIA_ARGS > "$log_file" 2>&1

exit_code=$?

# Record the end time and compute elapsed time in seconds
end_time=$(date +%s%N)
elapsed_time=$(echo "scale=2; ($end_time - $start_time) / 1000000000" | bc)
echo "TotalTime: $elapsed_time" | tee -a "$log_file"

exit $exit_code

