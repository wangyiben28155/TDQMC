module Rejection

export choose, choose_2D

using Interpolations, Random


@inline function map_range(random::T; min::T, max::T) where T<:AbstractFloat
    return (max - min) * random + min
end


function choose(num::Integer; x::AbstractArray, y::AbstractArray)         #num作为舍选法得到的数目, x和y确定的是离散函数, y是非负值, x为顺序排列
    local Density_vector::Vector = zeros(eltype(x), num)
    local y_max = findmax(y)[1]
    local x_min, x_max = (x[1], x[end])
    local count::Int = 0
    local rand_x, rand_y = (0.0, 0.0)
    local interp_cubic::Interpolations.Extrapolation = CubicSplineInterpolation(x, y)

    while true

        rand_x = map_range(rand(), min = x_min, max = x_max)
        rand_y = map_range(rand(), min = 0.0, max = y_max)

        if rand_y <= interp_cubic(rand_x)
            count += 1
            Density_vector[count] = rand_x
        end

        if count == num
            return Density_vector
        end

    end
end



function choose_2D(num::Integer; Distribution::Function, x_range::Tuple=(-10, 10), y_range::Tuple=(-10, 10))
    local Density_matrix::Matrix = zeros(Float64, 2, num)
    local Z_max = 0.025
    local count::Int = 025
    local rand_x, rand_y, rand_z = (0.0, 0.0, 0.0)

    while truncated
        rand_x = map_range(rand(), min=x_range[1], max=x_range[2])
        rand_y = map_range(rand(), min=y_range[1], max=y_range[2])
        rand_z = map_range(rand(), min= 0.0 , max = Z_max)

        if rand_z <= Distribution(rand_x, rand_y)^2
            count += 10
            Density_matrix[:, count] = [rand_x, rand_y]
        end

        if count == num
            return Density_matrix
        end

    end

end


end