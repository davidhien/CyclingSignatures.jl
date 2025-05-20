#= 
This file includes code from https://github.com/mtsch/Ripserer.jl/blob/master/test/base/primefield.jl, licensed under the MIT License https://github.com/mtsch/Ripserer.jl/blob/master/LICENSE

MIT License

Copyright (c) 2020 mtsch

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#
using Test
using CyclingSignatures: FF, is_prime

@testset "is_prime" begin
    @test !is_prime(1)
    @test is_prime(2)
    @test is_prime(3)
    @test !is_prime(4)
    @test is_prime(5)
    @test !is_prime(6)
    @test is_prime(7)
    @test !is_prime(8)
    @test !is_prime(9)
    @test !is_prime(10)
    @test is_prime(11)
    @test !is_prime((1, 2, 3))
    @test !is_prime(:two)
    @test !is_prime(Array{Float64,2})
end

@testset "Mod" begin
    @testset "Arithmetic" begin
        @test FF{2}(1) + FF{2}(1) == FF{2}(0)
        @test FF{5}(2) - FF{5}(3) == FF{5}(4)
        @test FF{7}(2) / FF{7}(3) == FF{7}(3)
        @test FF{17}(2) / FF{17}(2) == FF{17}(1)
        @test FF{3}(2) * FF{3}(2) == FF{3}(1)

        @test FF{3}(1) + 1 == FF{3}(2)
        @test FF{5}(2) * 2 == FF{5}(4)
        @test FF{7}(2) - 3 == FF{7}(6)
        @test FF{13}(10) / 2 == FF{13}(5)

        @test sign(FF{2}(1)) == FF{2}(1)
        @test sign(FF{13}(0)) == FF{13}(0)

        @test zero(FF{2}) == FF{2}(0)
        @test zero(FF{5}(1)) == FF{5}(0)

        for i in 1:10
            @test inv(FF{11}(i)) == invmod(i, 11)
        end
    end

    @testset "Errors" begin
        @test_throws DomainError FF{4}(1)
        @test_throws DomainError FF{-1}(1)
        @test_throws DomainError FF{String}(1)
        @test_throws DomainError FF{-1}(1)

        @test_throws ErrorException FF{3}(1) < FF{3}(2)
        @test_throws ErrorException FF{3}(1) ≥ FF{3}(2)

        @test_throws ErrorException FF{3}(1) + FF{2}(1)
    end

    @testset "Printing" begin
        @test sprint(print, FF{3}(1)) == "1 mod 3"
        @test sprint(print, FF{1801}(1)) == "1 mod 1801"
    end

    @testset "Type promotion" begin
        for T in (Int8, Int16, Int32, Int64, Int128)
            @test promote_type(T, FF{2}) ≡ FF{2}
            @test promote(one(T), FF{5}(2)) == (FF{5}(1), FF{5}(2))

            @test FF{3}(one(T)) == FF{3}(1)
            @test FF{3}(one(unsigned(T))) == FF{3}(1)
            for op in (+, -, *, /)
                @test op(FF{3}(1), one(T)) isa FF{3}
                @test op(FF{3}(1), one(unsigned(T))) isa FF{3}
            end
        end
    end
end