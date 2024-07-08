module PhilipsDataList

using CSV
using DataFrames
using FFTW
using NamedDims
using OffsetArrays
using ProgressBars

# Raw data from Philip's MR system is stored as single precision (complex) floats
const COMPLEX_ELTYPE = ComplexF32

# Types of "complex data vectors" in the .list file
const COMPLEX_DATA_VECTOR_TYPES = ["STD", "REJ", "PHX", "FRX", "NOI", "NAV", "DNA"]

# Dimensions of the data vectors in the .list file
const DIMENSIONS_STD = (:kx, :ky, :kz, :loca, :chan, :aver, :dyn, :mix, :card, :echo, :extr1, :extr2, :rf, :grad)

include("reader.jl")
include("kspace.jl")
include("utils.jl")

export read_data_list
export to_kspace

end
