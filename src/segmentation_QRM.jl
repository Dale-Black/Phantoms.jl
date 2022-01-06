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

