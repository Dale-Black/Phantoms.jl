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