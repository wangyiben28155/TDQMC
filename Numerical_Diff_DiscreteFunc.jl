module Discrete_Func_diff

export Three_Point, Five_Point

using SparseArrays

function Dimension_Judge(x::AbstractArray, y::AbstractArray)
    if length(x) != length(y)
        @error "The dimension of x and y is not equal." x_length = length(x) y_length = length(y)
    elseif length(x) <= 2
        @error "The length is too short."
    end
end

function Three_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}   #在程序中我们的空间步长都是一致的, 我们只需要计算
    local coefficient = 1 / (2 * dL)
    local num = length(x)
    local derivative_numrical_func = zeros(T, num)
    local Transform_Matrix::SparseMatrixCSC
    local Section_1 = [-3.0, 4.0, -1.0]
    local Section_2 = [-1.0, 1.0]
    local I, J, V = (zeros(2 * (num + 1)),
        zeros(2 * (num + 1)),
        zeros(2 * (num + 1)))

    V = [Section_1; repeat(Section_2, outer = num - 2); -reverse(Section_1)]
    I = [[1, 1, 1]; repeat(collect(2:num-1), inner = 2); [num, num, num]]
    J = [[1, 2, 3]; vcat([[i, i + 2] for i = 1:num-2]...); [num - 2, num - 1, num]]
    Transform_Matrix = sparse(I, J, V)

    derivative_numrical_func = coefficient .* (Transform_Matrix * y)

    return x, derivative_numrical_func
end

function Five_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}
    local coefficient = 1 / (12 * dL)
    local num = length(x)
    local derivative_numrical_func = zeros(Float64, num)
    local Transform_Matrix::SparseMatrixCSC
    local Section_1 = [-25.0, 48.0, -36.0, 16.0, -3.0]
    local Section_2 = [1.0, -8.0, 8.0, -1.0]
    local I, J, V = (zeros(4 * (num + 1)),
        zeros(4 * (num + 1)),
        zeros(4 * (num + 1)))

    V = [repeat(Section_1, outer = 2); repeat(Section_2, outer = num - 4); repeat(-reverse(Section_1), outer = 2)]
    I = [repeat([1, 2], inner = 5); repeat(collect(3:num-2), inner = 4); repeat([num - 1, num], inner = 5)]
    J = [collect(1:5); collect(2:6); vcat([[1, 2, 4, 6] .+ i for i = 0:num-5]...); collect(num-5:num-1); collect(num-4:num)]

    Transform_Matrix = sparse(I, J, V)
    println(size(Transform_Matrix))

    derivative_numrical_func = coefficient .* (Transform_Matrix * y)

    return x, derivative_numrical_func
end

end