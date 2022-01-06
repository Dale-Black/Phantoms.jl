@testset ExtendedTestSet "find_circle" begin
    @testset ExtendedTestSet "find_circle" begin
        test = find_circle([309, 309], [312, 200], [155, 155])
        answer = [212, 251]
        @test test == answer
    end

    @testset ExtendedTestSet "find_circle" begin
        test = find_circle([20, 20], [30, 10], [10, 20])
        answer = [15, 5]
        @test test == answer
    end
end

@testset ExtendedTestSet "angle_calc" begin
    @testset ExtendedTestSet "angle_calc" begin
        test = angle_calc(4, 3)
        answer = 0.9272952180016122
        @test test ≈ answer
    end

    @testset ExtendedTestSet "angle_calc" begin
        test = angle_calc(1, 2)
        answer = 0.4636476090008061
        @test test ≈ answer
    end
end

@testset ExtendedTestSet "create_circular_mask" begin
    @testset ExtendedTestSet "create_circular_mask" begin
        test = create_circular_mask(4, 4, [2, 2], 1);
        answer = Bool.([
            0  1  0  0
            1  1  1  0
            0  1  0  0
            0  0  0  0
        ])
        @test test == answer
    end
end

