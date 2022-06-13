"""
    mask_heart(
        header; 
        array_used=nothing, 
        radius_val=95, 
        slice_used_center=nothing
        )
Given a QRM Phantom with heart insert, this function will create a mask of the whole heart for
image processing purposes.
"""
function mask_heart(header, array_used, slice_used_center; radius_val=95)
    pixel_size = get_pixel_size(header)

    radius = (radius_val / 2) / pixel_size[1]
    central_image = copy(array_used[:, :, slice_used_center])
    central_image = Int.(central_image .< -200)
    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
    end
    central_image = mapwindow(median, central_image, (kern, kern))
    center = [size(central_image, 1) ÷ 2, size(central_image, 2) ÷ 2]
    a = copy(central_image)
    local point_1
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] + index] == 1 &&
            central_image[center[1] + index, center[2] + index + 5] == 1
        )
            point_1 = [center[1] + index, center[2] + index]
            break
        else
            a[center[1] + index, center[2] + index] = 2
        end
    end

    local point_2
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] - index] == 1 &&
            central_image[center[1] + index, center[2] - index - 5] == 1
        )
            point_2 = [center[1] + index, center[2] - index]
            break
        else
            a[center[1] + index, center[2] - index] = 2
        end
    end

    local point_3
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] - index, center[2] - index] == 1 &&
            central_image[center[1] - index, center[2] - index - 5] == 1
        )
            point_3 = [center[1] - index, center[2] - index]
            break
        else
            a[center[1] - index, center[2] - index] = 2
        end
    end

    center_insert = find_circle(point_1, point_2, point_3)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y - center_insert[1])^2)

    mask = dist_from_center .<= radius[1]
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
    end

    return masked_array, center_insert, mask
end

"""
	get_calcium_slices(dcm_array, header; calcium_threshold=130)

	Returns the slices that contain the calcium vessel inserts `slice_dict` 
	and the slices that contain the calcium calibration rod insert `large_index`.

	The calcium rod slice `large_index` usually omits the first and last slice 
	of the rod. So to account for all of the slices containing the calcium rod, 
	one would want a range like so: `(large_index - 1):(large_index + 1)`
"""
function get_calcium_slices(dcm_array, header; calcium_threshold=130)
    array = copy(dcm_array)
    array = Int.(array .> (1.1 * calcium_threshold))

    pixel_size = get_pixel_size(header)
    CCI_5mm_num_pixels = Int(round(π * (5 / 2)^2 / pixel_size[1]^2))
    cal_rod_num_pixels = Int(round(π * (20 / 2)^2 / pixel_size[1]^2))

    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
    end

    slice_dict = Dict()
    large_index = []
    cal_rod_dict = Dict()
    for idx in 1:size(array, 3)
        array_filtered = mapwindow(median, array[:, :, idx], (kern, kern))
        components = ImageComponentAnalysis.label_components(array_filtered)
        a1 = analyze_components(components, BasicMeasurement(; area=true, perimeter=true))
        a2 = analyze_components(components, BoundingBox(; box_area=true))
        df = leftjoin(a1, a2; on=:l)

        count_5mm = 0
        count = 0
        for row in eachrow(df)
            count += 1
            df_area = Int(round(row[:area]))

            r1_1 = Int(round(CCI_5mm_num_pixels * 0.6))
            r1_2 = Int(round(CCI_5mm_num_pixels * 1.5))
            r2_1 = Int(round(cal_rod_num_pixels * 0.7))
            r2_2 = Int(round(cal_rod_num_pixels * 1.3))

            if df_area in r1_1:r1_2
                count_5mm += 1
            elseif df_area in r2_1:r2_2
                indices = row[:box_indices]
                x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
                y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
                cal_rod_dict[count] = [x_point, y_point]
            end
        end

        if count_5mm > 0 && count_5mm < 4
            slice_dict[idx] = count_5mm
        end

        poppable_keys = []
        for key in cal_rod_dict
            start_coordinate = [key[2][1], key[2][2]]

            x_right = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] + x_right] == 1
                x_right += 1
            end

            x_left = 0
            while array_filtered[start_coordinate[1], start_coordinate[2] - x_left] == 1
                x_left += 1
            end

            y_top = 0
            while array_filtered[start_coordinate[1] + y_top, start_coordinate[2]] == 1
                y_top += 1
            end

            y_bottom = 0
            while array_filtered[start_coordinate[1] - y_bottom, start_coordinate[2]] == 1
                y_bottom += 1
            end

            x_dist = x_right + x_left
            y_dist = y_top + y_bottom

            range1 = round(0.7 * y_dist):round(1.2 * y_dist)
            if ((x_dist in range1) == false) ||
                ((round(0.7 * y_dist) == 0) && (round(1.2 * y_dist) == 0))
                push!(poppable_keys, key)
            else
                nothing
            end
        end

        for key in poppable_keys
            pop!(cal_rod_dict)
        end

        if length(cal_rod_dict) == 0
            nothing
        else
            append!(large_index, idx)
        end
    end
    return slice_dict, large_index
end

"""
	get_calcium_center_slices(dcm_array, slice_dict, large_index)

	Returns the slices that contain the calcium vessel inserts `slice_dict`
	and the center of the calcium calibration rod slice `flipped_index`
"""
function get_calcium_center_slices(dcm_array, slice_dict, large_index)
    flipped_index = Int(round(median(large_index)))
    edge_index = []
    if flipped_index < (size(dcm_array, 3) / 2)
        flipped = -1
        for element in large_index
            if element > (size(dcm_array, 3) / 2)
                append!(edge_index, element)
            end
        end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in minimum(edge_index):size(dcm_array, 3)
                try
                    delete!(slice_dict, index_edge)
                catch
                    nothing
                end
            end
            for element2 in edge_index
                deleteat!(large_index, findall(x -> x == element2, large_index))
            end
        end

        for element in 1:maximum(large_index)
            try
                delete!(slice_dict, element)
            catch
                nothing
            end
        end
    else
        flipped = 1
        for element in large_index
            if element < (size(dcm_array, 3) / 2)
                append!(edge_index, element)
            end
        end
        if length(edge_index) == 0
            nothing
        else
            for index_edge in 1:maximum(edge_index)
                try
                    delete!(slice_dict, index_edge)
                catch
                    nothing
                end
            end
            for element2 in edge_index
                deleteat!(large_index, findall(x -> x == element2, large_index))
            end
        end
        for element in minimum(large_index):size(dcm_array, 3)
            try
                delete!(slice_dict, element)
            catch
                nothing
            end
        end
    end
    return slice_dict, flipped, flipped_index
end

"""
    poppable_keys(flipped, flipped_index, header, slice_dict)

No idea what the function is
"""
function poppable_keys(flipped, flipped_index, header, slice_dict)
    SliceThickness = header[(0x0018, 0x0050)]
    poppable_keys = []
    if flipped == -1
        for key in slice_dict
            if key[1] > (flipped_index + (55 / SliceThickness))
                append!(poppable_keys, key)
            elseif flipped == 1
                for key in slice_dict
                    if key[1] < (flipped_index - (55 / SliceThickness))
                        append!(poppable_keys, key)
                    end
                end
            end
        end
    end
    for key in poppable_keys
        pop!(slice_dict)
    end
    return slice_dict
end

"""
    compute_CCI(dcm_array, header, slice_dict; calcium_threshold=130)

Main `CCI_calcium_image` function
"""
function compute_CCI(dcm_array, header, slice_dict, flipped; calcium_threshold=130)
    SliceThickness = header[(0x0018, 0x0050)]
    max_key, _ = maximum(zip(values(slice_dict), keys(slice_dict)))
    max_keys = []
    for key in slice_dict
        if key[2] == max_key
            append!(max_keys, key[1])
        end
    end
    slice_CCI = Int(floor(median(max_keys)))

    array = copy(dcm_array)
    array = Int.(array .> calcium_threshold)

    calcium_image = array .* dcm_array
    quality_slice = Int.(round(slice_CCI - flipped * (20 / SliceThickness)))

    cal_rod_slice = slice_CCI + (flipped * Int(30 / SliceThickness))

    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

"""
	mask_rod(dcm_array, header; calcium_threshold=130)
"""
function mask_rod(dcm_array, header; calcium_threshold=130)
    slice_dict, large_index = get_calcium_slices(
        dcm_array, header; calcium_threshold=calcium_threshold
    )
    slice_dict, flipped, flipped_index = get_calcium_center_slices(
        dcm_array, slice_dict, large_index
    )
    slice_dict = poppable_keys(flipped, flipped_index, header, slice_dict)
    calcium_image, slice_CCI, quality_slice, cal_rod_slice = compute_CCI(
        dcm_array, header, slice_dict, flipped; calcium_threshold=calcium_threshold
    )
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

"""
	calc_output(dcm_array, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3))

Calculate the output of a dcm_array
"""
function calc_output(
    dcm_array, header, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3)
)
    # Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
    PixelSpacing = get_pixel_size(header)
    SliceThickness = header[(0x0018, 0x0050)]
    CCI_min = Int((CCI_slice - round(5 / SliceThickness, RoundUp)) - 1)
    CCI_max = Int((CCI_slice + round(5 / SliceThickness, RoundUp)) + 1)
    central_CCI = Int(round(CCI_max - CCI_min) / 2)

    if CCI_min ≤ 0
        CCI_min = 1
    end
    if CCI_max > size(dcm_array, 3)
        CCI_max = size(dcm_array, 3)
    end

    CCI_array = copy(dcm_array[:, :, CCI_min:CCI_max])

    image_kernel = Int(round(3 / PixelSpacing[1]))
    if image_kernel % 2 == 0
        image_kernel += 1
    end

    CCI_array_binary = copy(CCI_array)
    CCI_array_binary = Int.(CCI_array_binary .> 1.0 * calcium_threshold)
    inp =
        CCI_array_binary[:, :, central_CCI - 1] +
        CCI_array_binary[:, :, central_CCI] +
        CCI_array_binary[:, :, central_CCI + 1]
    components = ImageComponentAnalysis.label_components(inp, comp_connect)
    a1 = analyze_components(components, BasicMeasurement(; area=true, perimeter=true))
    a2 = analyze_components(components, BoundingBox(; box_area=true))
    df = leftjoin(a1, a2; on=:l)
    centroids = []
    for row in eachrow(df)
        indices = row[:box_indices]
        x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
        y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
        push!(centroids, (x_point, y_point))
    end

    centroids = deleteat!(centroids, 1)

    i1 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI - 1], (image_kernel, image_kernel)
    )
    i2 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI], (image_kernel, image_kernel)
    )
    i3 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI + 1], (image_kernel, image_kernel)
    )

    image_for_center = i1 + i2 + i3

    components2 = ImageComponentAnalysis.label_components(image_for_center, comp_connect)
    components2 = Int.(components2 .> 0)
    components2 = ImageComponentAnalysis.label_components(components2, comp_connect)

    b1 = analyze_components(components2, BasicMeasurement(; area=true, perimeter=true))
    b2 = analyze_components(components2, BoundingBox(; box_area=true))
    df2 = leftjoin(b1, b2; on=:l)
    centroids2 = []
    for row in eachrow(df2)
        indices = row[:box_indices]
        x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
        y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
        push!(centroids2, (y_point, x_point))
    end

    output = length(unique(components2)), components2, df2, centroids2
    return output
end

"""
	center_points(output, header, tmp_center, CCI_slice)

Function ...
"""
function center_points(dcm_array, output, header, tmp_center, CCI_slice)
    PixelSpacing = PhantomSegmentation.get_pixel_size(header)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    sizes = []
    for row in eachrow(output[3])
        area = row[:area]
        append!(sizes, area)
    end

    centroids = output[4]
    largest = Dict()
    for index in 1:length(centroids)
        x = centroids[index][1]
        y = centroids[index][2]
        dist_loc = sqrt((tmp_center[2] - x)^2 + (tmp_center[1] - y)^2)
        dist_loc *= PixelSpacing[1]
        if dist_loc > 31
            largest[index] = [round(y), round(x)]
        else
            nothing
        end
    end

    max_dict = Dict()
    radius = round(2.5 / PixelSpacing[1], RoundUp)
    for key in largest
        tmp_arr = create_circular_mask(rows, cols, (key[2][2], key[2][1]), radius)
        tmp_arr = @. abs(tmp_arr * dcm_array[:, :, CCI_slice]) +
            abs(tmp_arr * dcm_array[:, :, CCI_slice - 1]) +
            abs(tmp_arr * dcm_array[:, :, CCI_slice + 1])
        tmp_arr = @. ifelse(tmp_arr == 0, missing, tmp_arr)
        max_dict[key[1]] = median(skipmissing(tmp_arr))
    end
    large1_index, large1_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large1_key)
    large2_index, large2_key = maximum(zip(values(max_dict), keys(max_dict)))
    pop!(max_dict, large2_key)
    large3_index, large3_key = maximum(zip(values(max_dict), keys(max_dict)))

    center1 = largest[large1_key]
    center2 = largest[large2_key]
    center3 = largest[large3_key]

    center = find_circle(center1, center2, center3)
    return center, center1, center2, center3
end

"""
	calc_centers(dcm_array, output, header, tmp_center, CCI_slice)

Function ...
"""
function calc_centers(dcm_array, output, header, tmp_center, CCI_slice; angle_factor=0)
    PixelSpacing = get_pixel_size(header)
    center, center1, center2, center3 = center_points(
        dcm_array, output, header, tmp_center, CCI_slice
    )
    centers = Dict()
    for size_index4 in (center1, center2, center3)
        center_index = size_index4
        side_x = abs(center[1] - center_index[1]) + angle_factor
        side_y = abs(center[2] - center_index[2]) + angle_factor

        angle = angle_calc(side_x, side_y)
        if (center_index[1] < center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] < center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] + (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] > center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (center_index[1] > center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] - (12.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (12.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (25 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (25 / PixelSpacing[2]) * cos(angle)),
            ]
        elseif (side_x == 0 && center_index[2] < center[2])
            medium_calc = [center_index[1], center_index[2] + (12.5 / PixelSpacing[2])]
            low_calc = [center_index[1], center_index[2] + (25 / PixelSpacing[2])]
        elseif (side_x == 0 && center_index[2] > center[2])
            medium_calc = [center_index[1], center_index[2] - (12.5 / PixelSpacing[2])]
            low_calc = [center_index[1], center_index[2] - (25 / PixelSpacing[2])]
        elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] - (12.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [center_index[1] - (25 / PixelSpacing[1]), center_index[2]]
        elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] + (12.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [(center_index[1] + (25 / PixelSpacing[1])), center_index[1]]
        else
            error("unknown angle")
        end

        if size_index4 == center1
            centers[:Large_HD] = Int.(round.(center_index))
            centers[:Medium_HD] = Int.(round.(medium_calc))
            centers[:Small_HD] = Int.(round.(low_calc))

        elseif size_index4 == center2
            centers[:Large_MD] = Int.(round.(center_index))
            centers[:Medium_MD] = Int.(round.(medium_calc))
            centers[:Small_MD] = Int.(round.(low_calc))

        elseif size_index4 == center3
            centers[:Large_LD] = Int.(round.(center_index))
            centers[:Medium_LD] = Int.(round.(medium_calc))
            centers[:Small_LD] = Int.(round.(low_calc))

        else
            nothing
        end
    end
    return centers
end

"""
	mask_inserts(
		dcm_array, masked_array, header, CCI_slice, center_insert; 
		calcium_threshold=130, comp_connect=trues(3, 3)
		)

Function ...
"""
function mask_inserts(
    dcm_array,
    masked_array,
    header,
    CCI_slice,
    center_insert;
    calcium_threshold=130,
    comp_connect=trues(3, 3),
    angle_factor=0
)
    output = calc_output(masked_array, header, CCI_slice, calcium_threshold, comp_connect)
    insert_centers = calc_centers(dcm_array, output, header, center_insert, CCI_slice; angle_factor=angle_factor)

    PixelSpacing = PhantomSegmentation.get_pixel_size(header)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])

    mask_L_HD = create_circular_mask(
        cols, rows, insert_centers[:Large_HD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_L_MD = create_circular_mask(
        cols, rows, insert_centers[:Large_MD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_L_LD = create_circular_mask(
        cols, rows, insert_centers[:Large_LD], (round(5 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_M_HD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_HD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_MD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_MD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_M_LD = create_circular_mask(
        cols,
        rows,
        insert_centers[:Medium_LD],
        (round(3 / PixelSpacing[1], RoundUp) / 2) + 1,
    )
    mask_S_HD = create_circular_mask(
        cols, rows, insert_centers[:Small_HD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_MD = create_circular_mask(
        cols, rows, insert_centers[:Small_MD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )
    mask_S_LD = create_circular_mask(
        cols, rows, insert_centers[:Small_LD], (round(1 / PixelSpacing[1], RoundUp) / 2) + 1
    )

    return transpose(mask_L_HD),
    transpose(mask_M_HD),
    transpose(mask_S_HD),
    transpose(mask_L_MD),
    transpose(mask_M_MD),
    transpose(mask_S_MD),
    transpose(mask_L_LD),
    transpose(mask_M_LD),
    transpose(mask_S_LD)
end
