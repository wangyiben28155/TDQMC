module Stitching

export extend_num, stitch!, stitch_Matrix


function extend_num(s_range::T, l_range::T, s_num::Int) where {T<:AbstractFloat}
    local l_num::Int = floor(Int, (l_range / s_range) * (s_num - 1) + 1)      #代表的是数组的长度,与接下来插入的索引有关

    if isodd(l_num) && isodd(s_num)
        return l_num
    else
        return @error "wrong l_range to extend"
    end
end


function stitch(A::Vector, l_num::Int)                             # A的长度发生改变变为l_num
    local L = length(A)
    local stitch_Vec = zeros(eltype(A), (l_num - L) >> 1)
    local B = deepcopy(A)

    splice!(B, L+1:L, stitch_Vec)
    splice!(B, 1:0, reverse!(stitch_Vec))

    return B

end


function stitch_Matrix(B::Array{T,3}, l_num::Int) where {T<:Complex{<:AbstractFloat}}      #目的是将小的数组转换为Matrix{Vector}的结构
    local a, b, c = size(B)
    local Raw_Matrix::Matrix{<:Vector} = [zeros(eltype(B), l_num) for i in 1:b, j in 1:c]

    if l_num < a
        return @error "l_num must be larger than a"
    end

    for j in 1:c, i in 1:b
        Raw_Matrix[i, j][:] = stitch(B[:, i, j], l_num)
    end

    return Raw_Matrix
end


end