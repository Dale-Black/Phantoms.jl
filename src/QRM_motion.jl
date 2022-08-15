function calc_output_motion(
    dcm_array, header, slices, calcium_threshold=130, comp_connect=trues(3, 3)
)
    # Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
    PixelSpacing = PhantomSegmentation.get_pixel_size(header)
    SliceThickness = header[(0x0018, 0x0050)]
    CCI_min = slices[1] - 1
    CCI_max = slices[2] + 2
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

function mask_inserts_motion(dcm_array, slices1, slices2; threshold=115, radius=5)
	output1 = calc_output_motion(masked_array, header, slices1, threshold)
	output2 = calc_output_motion(masked_array, header, slices2, threshold)
	center1 = output1[4][1]
	center2 = output2[4][1]

	mask1 = create_circular_mask(size(dcm_array)[1:2]..., center1, radius)
	mask2 = create_circular_mask(size(dcm_array)[1:2]..., center2, radius)
	return mask1, mask2
end
