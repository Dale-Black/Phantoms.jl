function mask_inserts_motion(dcm_array, slices1, slices2; threshold=115, radius=5)
	output1 = calc_output(dcm_array, header, slices1, threshold)
	output2 = calc_output(dcm_array, header, slices2, threshold)
	center1 = output1[4][1]
	center2 = output2[4][1]

	mask1 = create_circular_mask(size(dcm_array)[1:2]..., center1, radius)
	mask2 = create_circular_mask(size(dcm_array)[1:2]..., center2, radius)
	return mask1, mask2
end
