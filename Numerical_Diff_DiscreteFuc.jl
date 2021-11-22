module Discrete_Func_diff

export Three_Point, Five_Point

using SparseArrays

macro Dimension_Judge(x::Vector{Any}, y::Vector{Any})
    quote
        if length($x) != length($y)
            @error "The dimension of x and y is not equal." x_length = length(x) y_length = length(y)
        end
    end
end

function Three_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}           #在程序中我们的空间步长都是一致的, 我们只需要计算
    local coefficient = 1 / (2 * dL)
    local derivative_numrical_func = zeros(length(x))
    local Transform_Matrix = sparse(1)


    return x, derivative_numrical_func
end

function Five_Point(x::Vector{T}, y::Vector{T}; dL = Parameter.dx) where {T<:AbstractFloat}
    local coefficient = 1 / (2 * dL)
    local derivative_numrical_func = zeros(length(x))
    local Transform_Matrix = sparse()

    if length(x) != length(y)
        @error "The dimension of x and y is not equal."
    end



end