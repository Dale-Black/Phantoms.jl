@testset ExtendedTestSet "mask_heart" begin
    path = string(cd(pwd, "..") , "/images/Large_rep1")
    dcms = dcmdir_parse(path)
    dcm_array = load_dcm_array(dcms)
    header = dcms[1].meta

    @testset ExtendedTestSet "mask_heart" begin
        masked_array, center_insert, mask = mask_heart(
            header, dcm_array, size(dcm_array, 3) ÷ 2
            )
        rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
        Y, X = collect(1:512), collect(1:512)'
        dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y-center_insert[1])^2)
        answer = dist_from_center .<= 97.3360655737705
        @test mask == answer
    end
end

@testset ExtendedTestSet "get_indices" begin
    path = string(cd(pwd, "..") , "/images/Large_rep1")
    dcms = dcmdir_parse(path)
    dcm_array = load_dcm_array(dcms)
    header = dcms[1].meta
    masked_array, center_insert, mask = mask_heart(
            header, dcm_array, size(dcm_array, 3) ÷ 2
            )

    @testset ExtendedTestSet "get_indices" begin
        slice_dict, large_index = get_indices(masked_array, header);
        answer_dict = Dict(25=>2, 26=>2, 24=>1)
        @test slice_dict == answer_dict
    end

    @testset ExtendedTestSet "get_indices" begin
        slice_dict, large_index = get_indices(masked_array, header);
        answer_index = [10, 11, 12, 13, 14, 15, 16]
        @test large_index == answer_index
    end
end