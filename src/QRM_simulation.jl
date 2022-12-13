"""
    center_points_simulation(dcm_array, output, header, tmp_center, CCI_slice)

Function ...
"""
function center_points_simulation(dcm_array, output, header, tmp_center, CCI_slice)
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
        tmp_arr = @. abs(tmp_arr * dcm_array[:,:,CCI_slice]) + abs(tmp_arr * dcm_array[:,:,CCI_slice - 1]) + abs(tmp_arr * dcm_array[:,:,CCI_slice + 1])
        tmp_arr = @. ifelse(tmp_arr == 0, missing, tmp_arr)
        max_dict[key[1]] = median(skipmissing(tmp_arr))
	end
    _, large1_key = maximum(zip(values(max_dict), keys(max_dict)))
    center1 = largest[large1_key]
	
	center = vec(tmp_center')
	p1 = vec(center1')
	offset = -2π/180
	offset2 = -5π/180
	p2, p3 = find_triangle_points(p1, center; offset=offset, offset2=offset2)
	cent1, cent2, cent3 = vec(p1'), vec(p2'), vec(p3')
	
    center = find_circle(cent1, cent2, cent3)
	return center, cent1, cent2, cent3
end

"""
    calc_centers_simulation(dcm_array, output, header, tmp_center, CCI_slice)

Function ...
"""
function calc_centers_simulation(dcm_array, output, header, tmp_center, CCI_slice)
    PixelSpacing = PhantomSegmentation.get_pixel_size(header)
    center, center1, center2, center3 = center_points_simulation(
        dcm_array, output, header, tmp_center, CCI_slice
    )
    centers = Dict()
    for center_index in (center1, center2, center3)
        side_x = abs(center[1] - center_index[1])
        side_y = abs(center[2] - center_index[2])
        angle = angle_calc(side_x, side_y)
        if (center_index[1] < center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] + (10.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (10.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (17 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (17 / PixelSpacing[2]) * cos(angle)),
            ]

        elseif (center_index[1] < center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] + (10.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (10.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] + (17 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (17 / PixelSpacing[2]) * cos(angle)),
            ]

        elseif (center_index[1] > center[1] && center_index[2] < center[2])
            medium_calc = [
                center_index[1] - (10.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (10.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (17 / PixelSpacing[1]) * sin(angle),
                (center_index[2] + (17 / PixelSpacing[2]) * cos(angle)),
            ]

        elseif (center_index[1] > center[1] && center_index[2] > center[2])
            medium_calc = [
                center_index[1] - (10.5 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (10.5 / PixelSpacing[2]) * cos(angle)),
            ]
            low_calc = [
                center_index[1] - (17 / PixelSpacing[1]) * sin(angle),
                (center_index[2] - (17 / PixelSpacing[2]) * cos(angle)),
            ]

        elseif (side_x == 0 && center_index[2] < center[2])
            medium_calc = [center_index[1], center_index[2] + (10.5 / PixelSpacing[2])]
            low_calc = [center_index[1], center_index[2] + (17 / PixelSpacing[2])]

        elseif (side_x == 0 && center_index[2] > center[2])
            medium_calc = [center_index[1], center_index[2] - (10.5 / PixelSpacing[2])]
            low_calc = [center_index[1], center_index[2] - (17 / PixelSpacing[2])]

        elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] - (10.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [center_index[1] - (17 / PixelSpacing[1]), center_index[2]]

        elseif (center_index[1] > center[1] && side_y == 0)
            medium_calc = [center_index[1] + (10.5 / PixelSpacing[1]), center_index[2]]
            low_calc = [(center_index[1] + (17 / PixelSpacing[1])), center_index[1]]

        else
            error("unknown angle")
        end

        if center_index == center1
            centers[:Large_HD] = Int.(round.(center_index))
            centers[:Medium_HD] = Int.(round.(medium_calc))
            centers[:Small_HD] = Int.(round.(low_calc))

        elseif center_index == center2
            centers[:Large_MD] = Int.(round.(center_index))
            centers[:Medium_MD] = Int.(round.(medium_calc))
            centers[:Small_MD] = Int.(round.(low_calc))

        elseif center_index == center3
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
    mask_inserts_simulation(
        dcm_array,
        masked_array,
        header,
        CCI_slice,
        center_insert;
        calcium_threshold=130,
        comp_connect=trues(3, 3),
    )

Function ...
"""
function mask_inserts_simulation(
    dcm_array,
    masked_array,
    header,
    CCI_slice,
    center_insert;
    calcium_threshold=130,
    comp_connect=trues(3, 3),
)
    output = calc_output(masked_array, header, CCI_slice, calcium_threshold, comp_connect)
    insert_centers = calc_centers_simulation(
        dcm_array, output, header, center_insert, CCI_slice
    )

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