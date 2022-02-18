module Stitching

export extend_num, stitch!


function extend_num(s_range::T, l_range::T, s_num::Int) where {T<:AbstractFloat}
    local l_num::Int = floor(Int, (l_range / s_range) * (s_num - 1) + 1)      #代表的是数组的长度,与接下来插入的索引有关

    if isodd(l_num) && isodd(s_num)
        return l_num
    else
        return @error "wrong l_range to extend"
    end
end


function stitch!(A::Vector; s_num::Int, l_num::Int)                             # A的长度发生改变变为l_num
    local stitch_Vec = zeros(eltype(A), (l_num - s_num) >> 1)
    local L = length(A)

    if L != s_num
        return @error "Length is not equal."
    else
        splice!(A, L+1:L, stitch_Vec)
        splice!(A, 1:0, reverse!(stitch_Vec))

        return A
    end
end


end