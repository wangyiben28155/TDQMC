module Find_nearest_k

export find_k_index


function find_eg_index(target::T; x::AbstractVector{T}) where {T<:Number}                #这里x为从小到大排序好的数组
    local Start::Integer = 1
    local End::Integer = length(x)
    local Middle::Integer = (Start + End) >> 1

    while true
        if Start > End
            return End + 1
        end

        Middle = (Start + End) >> 1

        if x[Middle] == target
            return Middle
        elseif x[Middle] < target
            Start = Middle + 1
        else
            End = Middle - 1
        end
    end
end


function find_closest_index(target::T; x::AbstractVector{T}) where {T<:Number}
    local eg_index::Integer = find_eg_index(target, x = x)     #这里直接计算得到大于target值得坐标

    if eg_index == length(x) + 1
        return eg_index - 1
    elseif eg_index == 1
        return 1
    else
        if target - x[eg_index-1] <= x[eg_index] - target
            return eg_index - 1
        else
            return eg_index
        end
    end

end

function get_ele(x::AbstractVector{T}, index::Integer) where {T<:Number}
    if index < 1 || index > length(x)
        return Inf
    else
        return x[index]
    end
end


function find_k_index(target::T1; x::AbstractVector{T1}, k::T2) where {T1<:Number,T2<:Integer}       #找到最接近得k个元素,并返回一个索引的数组
    local k_index::Vector{T2} = zeros(T2, k)
    local closest_index::T2 = find_closest_index(target, x = x)
    local left_index = closest_index - 1
    local right_index = closest_index + 1

    k_index[1] = closest_index

    for i = 2:k
        if abs(target - get_ele(x, left_index)) <= abs(get_ele(x, right_index) - target)
            k_index[i] = left_index
            left_index -= 1
        else
            k_index[i] = right_index
            right_index += 1
        end
    end

    sort!(k_index)

    return k_index
end


end