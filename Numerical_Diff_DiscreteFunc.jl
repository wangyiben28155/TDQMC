module Discrete_Func_diff

export Three_Point, Five_Point, Derivative_2

#using ..TDQMC
using SparseArrays


function Three_Point(y::Vector{<:Union{Complex{T},T}}; dL::T) where {T<:AbstractFloat}   #在程序中我们的空间步长都是一致的
    local coefficient = 1.0 / (2.0 * dL)
    local num = length(y)
    local derivative_numrical_func = zeros(eltype(y), num)
    local Transform_Matrix::SparseMatrixCSC = spzeros(T, num, num)
    local Section_1 = [-3.0, 4.0, -1.0]
    local Section_2 = [-1.0, 1.0]
    local I, J, V = (zeros(2 * (num + 1)),
        zeros(2 * (num + 1)),
        zeros(T, 2 * (num + 1)))

    V = [Section_1; repeat(Section_2, outer = num - 2); -reverse(Section_1)]
    I = [[1, 1, 1]; repeat(collect(2:num-1), inner = 2); [num, num, num]]
    J = [[1, 2, 3]; vcat([[i, i + 2] for i = 1:num-2]...); [num - 2, num - 1, num]]
    Transform_Matrix = sparse(I, J, V)

    derivative_numrical_func = coefficient .* (Transform_Matrix * y)

    return derivative_numrical_func
end

function Five_Point(y::Vector{<:Union{Complex{T},T}}; dL::T) where {T<:AbstractFloat}        #五步法计算一阶导数
    local coefficient = 1.0 / (12.0 * dL)
    local num = length(y)
    local derivative_numrical_func = zeros(eltype(y), num)
    local Transform_Matrix::SparseMatrixCSC = spzeros(T, num, num)
    local Section_1 = [-25.0, 48.0, -36.0, 16.0, -3.0]
    local Section_2 = [1.0, -8.0, 8.0, -1.0]
    local I, J, V = (zeros(4 * (num + 1)),
        zeros(4 * (num + 1)),
        zeros(T, 4 * (num + 1)))

    V = [repeat(Section_1, outer = 2); repeat(Section_2, outer = num - 4); repeat(-reverse(Section_1), outer = 2)]
    I = [repeat([1, 2], inner = 5); repeat(collect(3:num-2), inner = 4); repeat([num - 1, num], inner = 5)]
    J = [collect(1:5); collect(2:6); repeat([1, 2, 4, 5], outer = num - 4) .+ repeat(collect(0:num-5), inner = 4); collect(num-5:num-1); collect(num-4:num)]

    Transform_Matrix = sparse(I, J, V)

    derivative_numrical_func = coefficient .* (Transform_Matrix * y)

    return derivative_numrical_func
end


function Derivative_2(y::Vector{<:Union{Complex{T},T}}; dL::T) where {T<:AbstractFloat}
    local coefficient = 1 / (12 * dL^2)
    local num = length(y)
    local derivative_two_func = zeros(eltype(y), num)
    local Transform_Matrix::SparseMatrixCSC = spzeros(T, num, num)
    local Section_1 = [35.0, -104.0, 114.0, -56.0, 11.0]
    local Section_2 = [11.0, -20.0, 6.0, 4.0, -1.0]
    local Section_3 = [-1.0, 16.0, -30.0, 16.0, -1.0]
    local I, J, V = (zeros(5 * num),
        zeros(5 * num),
        zeros(T, 5 * num))

    V = [Section_1; Section_2; repeat(Section_3, outer = num - 4); reverse(Section_2); reverse(Section_1)]
    I = repeat(collect(1:num), inner = 5)
    J = [repeat(collect(1:5), outer = 2); repeat(collect(1:5), outer = num - 4) .+ repeat(collect(0:num-5), inner = 5); repeat(collect(num-4:num), outer = 2)]

    Transform_Matrix = sparse(I, J, V)

    derivative_two_func = coefficient .* (Transform_Matrix * y)

    return derivative_two_func
end

end