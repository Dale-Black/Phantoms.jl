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
function mask_heart(
	header, array_used, slice_used_center; 
	radius_val=95, 
	)
	
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
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] + index] == 1 && central_image[center[1] + index, center[2] + index + 5] == 1) 
			point_1 = [center[1] + index, center[2] + index]
			break
		else
			a[center[1] + index, center[2] + index] = 2
		end
	end
    
	local point_2
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] + index, center[2] - index] == 1 && central_image[center[1] + index, center[2] - index - 5] == 1) 
			point_2 = [center[1] + index, center[2] - index]
			break
		else
			a[center[1] + index, center[2] - index] = 2
		end
	end
	
	local point_3
	for index in 1:size(central_image, 2) ÷ 2
		if (central_image[center[1] - index, center[2] - index] == 1 && central_image[center[1] - index, center[2] - index - 5] == 1)
			point_3 = [center[1] - index, center[2] - index]
			break
		else
			a[center[1] - index, center[2] - index] = 2
		end
	end
	
    center_insert = find_circle(point_1, point_2, point_3)
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y-center_insert[1])^2)

    mask = dist_from_center .<= radius[1]  
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
	end

    return masked_array, center_insert, mask
end

"""
	get_indices(dcm_array, header; calcium_threshold=130)

Get the indices of a dcm_array
"""
function get_indices(dcm_array, header; calcium_threshold=130)
    array = copy(dcm_array)
    array = Int.(array .> (1.1 * calcium_threshold))

	pixel_size = get_pixel_size(header)
    CCI_5mm_num_pixels = Int(round(π * (5/2)^2 / pixel_size[1]^2))
    cal_rod_num_pixels = Int(round(π * (20/2)^2 / pixel_size[1]^2))
    
	kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
	end
    
    slice_dict = Dict()
    large_index = []
    cal_rod_dict = Dict()
    for idx in 1:size(array, 3)
		array_filtered = mapwindow(median, array[:,:,idx], (kern, kern))
        components = ImageComponentAnalysis.label_components(array_filtered)
        count_5mm = 0
		a1 = analyze_components(
			components, BasicMeasurement(area=true, perimeter=true)
		)
		a2 = analyze_components(components, BoundingBox(box_area = true))
		df = leftjoin(a1, a2, on = :l)
		count = 0
		for row in eachrow(df)
			count += 1
			df_area = Int(round(row[:area]))
			if df_area in Int(round(CCI_5mm_num_pixels * 0.6)):Int(round(CCI_5mm_num_pixels * 1.5))
				count_5mm += 1
			elseif df_area in Int(round(cal_rod_num_pixels * 0.7)):Int(round(cal_rod_num_pixels * 1.3))
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
            
            if ((x_dist in round(0.7 * y_dist):round(1.2 * y_dist)) == false) || ((round(0.7 * y_dist) == 0) && (round(1.2 * y_dist) == 0))
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
	find_edges(dcm_array, slice_dict, large_index)

Find the edges of a dcm_array
"""
function find_edges(dcm_array, slice_dict, large_index)
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
				deleteat!(large_index, findall(x->x==element2, large_index))
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
				deleteat!(large_index, findall(x->x==element2, large_index))
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
	SliceThickness = header[(0x0018,0x0050)]
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
	SliceThickness = header[(0x0018,0x0050)]
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
    quality_slice = round(slice_CCI - flipped * (20 / SliceThickness))

    cal_rod_slice = slice_CCI + (flipped * Int(30 / SliceThickness))
    
    return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end

"""
    CCI_calcium_image(dcm_array, header; calcium_threshold=130)
"""
function CCI_calcium_image(dcm_array, header; calcium_threshold=130)
    slice_dict, large_index = get_indices(
		dcm_array, header; 
		calcium_threshold=calcium_threshold
	)
    slice_dict, flipped, flipped_index = find_edges(
		dcm_array, slice_dict, large_index
	)
    slice_dict = poppable_keys(flipped, flipped_index, header, slice_dict)
    calcium_image, slice_CCI, quality_slice, cal_rod_slice = compute_CCI(
		dcm_array, header, slice_dict, flipped; calcium_threshold=calcium_threshold
	)
	return calcium_image, slice_CCI, quality_slice, cal_rod_slice
end