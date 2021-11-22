module Discrete_Func_diff

export Three_Point, Five_Point

using SparseArrays

function Dimension_Judge(x::Vector{<:Any}, y::Vector{<:Any})
    if length(x) != length(y)
        @error "The dimension of x and y is not equal." x_length = length(x) y_length = length(y)
    elseif length(x) <= 2
        @error "The length is too short."
    end
end

function Three_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}           #在程序中我们的空间步长都是一致的, 我们只需要计算
    local coefficient = 1 / (2 * dL)
    local num = length(x)
    local derivative_numrical_func = zeros(num)
    local Transform_Matrix::SparseMatrixCSC
    local Section_1 = [-3.0, 4.0, -1.0]
    local Section_2 = [-1.0, 1.0]
    local I, J, V = (zeros(2 * (num + 1)), 
                     zeros(2 * (num + 1)), 
                     zeros(2 * (num + 1)))

    V = [Section_1;repeat(Section_2, outer = 2 * (n-2));-reverse(Section_1)]
    I = [[1,1,1];repeat(collect(2:num-1), inner = 2);[num, num, num]]
    J = [[1,2,3];]






    derivative_numrical_func = Transform_Matrix * y

    return x, coefficient .* derivative_numrical_func
end

function Five_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}
    local coefficient = 1 / (2 * dL)
    local derivative_numrical_func = zeros(length(x))
    local Transform_Matrix = sparse()



end