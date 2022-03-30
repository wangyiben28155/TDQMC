module Run_Anti_Ground

include("Main_module.jl")
using .TDQMC

using Dates, Random

const Ensemble_num, Electron_num = (20, 2)
const μ = [-1.5, 1.5]
const σ = [1.0, 1.0]
const L, x_num = (30.0, 301)


tic = now()

function sqrt_Normal(x; μ::Number=0, σ::Number=1)
    return sqrt(exp(-(x - μ)^2 / 2σ^2) / sqrt(2pi) * σ)
end

function Antisymmetry(x1::T, x2::T; μ1::Number=μ[1], μ2::Number=μ[2], σ1::Number=σ[1], σ2::Number=σ[2]) where {T<:AbstractFloat}
    return sqrt_Normal(x1, μ=μ1, σ=σ1) * sqrt_Normal(x2, μ=μ2, σ=σ2) - sqrt_Normal(x1, μ=μ2, σ=σ2) * sqrt_Normal(x2, μ=μ1, σ=σ1)
end

@inline function map_range(random::T; min::T, max::T) where {T<:AbstractFloat}
    return (max - min) * random + min
end

function choose_2D(num::Integer; Distribution::Function, x_range::Tuple=(-10.0, 10.0), y_range::Tuple=(-10.0, 10.0))
    local Density_matrix::Matrix = zeros(Float64, 2, num)
    local Z_max = 0.5             #这里是通过作图看出来的上限
    local count::Int = 0
    local rand_x, rand_y, rand_z = (0.0, 0.0, 0.0)

    Random.seed!(1)

    while true
        rand_x = map_range(rand(), min=x_range[1], max=x_range[2])
        rand_y = map_range(rand(), min=y_range[1], max=y_range[2])
        rand_z = map_range(rand(), min=0.0, max=Z_max)

        if rand_z <= Distribution(rand_x, rand_y)^2
            count += 1
            Density_matrix[:, count] = [rand_x, rand_y]
        end

        if count == num
            return Density_matrix
        end

    end

end

function initializer_guideWave(n::T; μ::AbstractVector=μ, σ::AbstractVector=σ) where {T<:Integer}                                           #这里因为算的是一维的波函数, 所以返回的矩阵为三维矩阵,第一维为波函数,第二维为系综粒子数,第三维为电子数
    return complex(sqrt_Normal.(LinRange(-L, L, x_num), μ=μ[n], σ=σ[n]))
end


function AntiSym(Tr::Matrix{T}) where {T<:AbstractFloat}
    local Guide_Wave::Matrix{<:Vector{<:Complex{T}}} = [zeros(ComplexF64, x_num) for i in 1:Electron_num, j in 1:Ensemble_num]

    for i in 1:size(Tr, 2)
        if Tr[2, i] >= Tr[1, i]
            Guide_Wave[:, i] = [initializer_guideWave(i) for i in 2:-1:1]
        else
            Guide_Wave[:, i] = [initializer_guideWave(i) for i in 1:2]
        end
    end
    
    return Guide_Wave
end

tr = choose_2D(Ensemble_num, Distribution=Antisymmetry)

P = Parameter{Float64,Int64}()
Dy = Dynamics{Float64,Int64}(Trajectory=tr, Guide_Wave=AntiSym(tr))      #注意这里要将Main_module的初始波函数的σ改为1和2

parallel_CTE!(P, Dy)

Sum_Energy = sum(Dy.Energy) / length(Dy.Energy)

toc = now()

println("The Total Energy is $(Sum_Energy)")

println("$(convert(DateTime, toc-tic))")

end
