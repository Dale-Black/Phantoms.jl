module Phantoms

using ImageFiltering
using Statistics
using ImageComponentAnalysis
using DataFrames

include("segmentation_QRM.jl")
include("utils.jl")

export 
    # segmentation_QRM 
    mask_heart,
    get_indices,
    find_edges,
    poppable_keys,
    compute_CCI,
    CCI_calcium_image,

    # utils.jl 
    find_circle,
    get_pixel_size

end
