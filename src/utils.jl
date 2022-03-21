"""
    findCircle(point_1, point_2, point_3)

Given three points, return center points of the circle with [y, x] order instead of 
the more common [x, y] ordering.
"""
function find_circle(point_1, point_2, point_3)
    x1, y1 = point_1
    x2, y2 = point_2
    x3, y3 = point_3

    x12 = x1 - x2
    x13 = x1 - x3
    y12 = y1 - y2
    y13 = y1 - y3
    y31 = y3 - y1
    y21 = y2 - y1
    x31 = x3 - x1
    x21 = x2 - x1

    sx13 = x1^2 - x3^2
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2
    sy21 = y2^2 - y1^2

    f = (
        ((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) ÷
        (2 * ((y31) * (x12) - (y21) * (x13)))
    )

    g = (
        ((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) ÷
        (2 * ((x31) * (y12) - (x21) * (y13)))
    )

    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where center is (h = -g, k = -f)  
    center_insert = [-g, -f]

    return center_insert
end

"""
    get_pixel_size(header)

Get the pixel size of a DICOM image given the DICOM header
"""
function get_pixel_size(header)
    pixel_size = try
        header[(0x0028, 0x0030)]
    catch
        FOV = header[(0x0018, 0x1100)]
        matrix_size = header[(0x0028, 0x0010)]

        pixel_size = FOV / matrix_size
    end
    return pixel_size
end

"""
    angle_calc(side1, side2)

Calculate angle between two sides of rectangular triangle
"""
function angle_calc(side1, side2)
    if side1 == 0
        angle = 0
    elseif side2 == 0
        angle = π / 2
    else
        angle = atan(side1 / side2)
    end

    return angle
end

"""
    create_circular_mask(h, w, center_circle, radius_circle)

Create a circle mask
"""
function create_circular_mask(h, w, center_circle, radius_circle)
    Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]) .^ 2 .+ (Y .- center_circle[2]) .^ 2)
    mask = dist_from_center .<= radius_circle
    return mask
end

"""

"""
function mass_calibration(
    dcm_array, center_large_LD, center, cal_rod_slice, rows, cols, spacing
)
    center_LD = center_large_LD
    dist_x = abs(center_LD[1] - center[1])
    dist_y = abs(center_LD[2] - center[2])

    if dist_x == 0
        mass_center_x = center[1]
        if center_LD[2] > center[2]
            mass_center_y = round(center[1] - round(23 / spacing[1], RoundUp), RoundUp)
        else
            mass_center_y = round(center[1] + round(23 / spacing[1], RoundUp), RoundUp)
        end
    elseif dist_y == 0
        mass_center_y = center[2]
        if center_LD[1] > center[1]
            mass_center_x = round(center[1] - round(23 / spacing[1], RoundUp), RoundUp)
        else
            mass_center_x = round(center[0] + round(23 / spacing[1], RoundUp), RoundUp)
        end

    else
        mass_angle = atan(dist_y / dist_x)
        dist_x = (23 / spacing[1]) * cos(mass_angle)
        dist_y = (23 / spacing[1]) * sin(mass_angle)

        if (center_LD[1] < center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] < center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        end
    end

    mass_cal_center = [mass_center_y, mass_center_x]
    x_distance = abs(center[1] - mass_cal_center[2])
    angled_distance = sqrt(
        (center[1] - mass_cal_center[2])^2 + (center[2] - mass_cal_center[2])^2
    )
    angle_0_200HA = acos(x_distance / angled_distance) * 180 / π
    mask_0HU = PhantomSegmentation.create_circular_mask(
        cols, rows, mass_cal_center, Int(round(6.9 / spacing[1]))
    )
    masked_0HU = mask_0HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count = length(findall(x -> x != 0, masked_0HU))
    mean_0HU = sum(masked_0HU) / nonzero_count
    std_0HU = 0

    for voxel in vec(masked_0HU)
        if voxel != 0
            std_0HU += (voxel - mean_0HU)^2
        end
    end

    std_0HU = sqrt(std_0HU / nonzero_count)
    mask_200HU = PhantomSegmentation.create_circular_mask(
        cols, rows, (center[2], center[1]), Int(round(6.9 / spacing[1]))
    )
    masked_200HU = mask_200HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count_200HU = length(findall(x -> x != 0, masked_200HU))
    mean_200HU = sum(masked_200HU) / nonzero_count_200HU
    mass_cal_factor = 0.2 / (mean_200HU - mean_0HU)
    water_rod_metrics = mean_0HU, std_0HU

    return mass_cal_factor, angle_0_200HA, water_rod_metrics
end