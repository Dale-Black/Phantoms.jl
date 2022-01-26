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
