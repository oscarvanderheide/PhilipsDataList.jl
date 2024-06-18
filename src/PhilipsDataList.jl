module PhilipsDataList

using CSV
using DataFrames
using ProgressBars

# Raw data from Philip's MR system is stored as single precision (complex) floats
const COMPLEX_ELTYPE = ComplexF32

# Types of "complex data vectors" in the .list file
const COMPLEX_DATA_VECTOR_TYPES = ["STD", "REJ", "PHX", "FRX", "NOI", "NAV", "DNA"]

include("reader.jl")
include("utils.jl")

export read_data_list

end
