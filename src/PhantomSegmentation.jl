module PhantomSegmentation

using ImageFiltering
using Statistics
using LinearAlgebra
using ImageComponentAnalysis
using ImageSegmentation
using DataFrames

include("QRM.jl")
include("QRM_simulation.jl")
include("QRM_motion.jl")
include("utils.jl")

export
    #QRM_simulation.jl
    calc_output_motion,
    mask_inserts_motion,

    #QRM_simulation.jl
    center_points_simulation,
    calc_centers_simulation,
    mask_inserts_simulation,

    # QRM.jl
    mask_heart,
    get_calcium_slices,
    get_calcium_center_slices,
    poppable_keys,
    compute_CCI,
    mask_rod,
    calc_output,
    center_points,
    calc_centers,
    mask_inserts,

    # utils.jl 
    find_circle,
    get_pixel_size,
    angle_calc,
    create_circular_mask,
    find_triangle_points,
    mass_calibration

end
