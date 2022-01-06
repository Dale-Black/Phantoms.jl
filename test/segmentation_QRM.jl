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
