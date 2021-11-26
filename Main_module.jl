module TDQMC

export Dynamics, Parameter

using Distributions, Random                   #用来初始化初始的波函数和系综粒子分布
import Base.@kwdef

const Ensemble_num, Electron_num = (50, 1)     #设定一些计算的参数
const L, x_num, step_t = (350.0, 1751, 13001) 
const μ, σ = (0, 30)


function generate_distribution(; x::Int = 123)
    Random.seed!(x)                                               #确定下一次使用rand函数的撒点是固定的, 每次调用程序也能保证下面的分布是确定的
    Initial_Distriution = rand(Normal(μ, σ), Ensemble_num)         #将确定的撒点分布储存起来, 用于下面实例化的赋值

    return Initial_Distriution
end

@kwdef mutable struct Dynamics{T<:AbstractFloat}
    Guide_Wave::Vector{Complex{T}} = complex(collect(pdf(Normal(μ, σ), LinRange(-L, L, x_num))))
    trajectory::Vector{T} = generate_distribution()                       #先给定一个初始的轨迹分布,目前只是计算氢原子的情况,所以先只使用向量
    Time::Union{T,Complex{T}} = 0.0
end


@kwdef struct Parameter{T<:AbstractFloat}              #用来控制计算参数的
    space_N::Int64 = x_num                                    #划分的格点的总数,后面做离散傅里叶变换的时候会用得到
    scope::T = L                                        #确定波函数的计算范围为-scope到+scope
    Δx::T = 2 * scope / (N - 1)                             #波函数的离散的空间间隔
    Δt::Union{T,Complex{T}} = 0.05                     #划分的时间间隔, 尝试时间迭代区间, 因为考虑到虚时演化, 所以类型设定为复数
    Step_t::Int64 = step_t
    sampling::Vector{T} = collect(LinRange(-scope, scope, N))
end

include("Numerical_Diff_DiscreteFunc.jl")      #用来做数值微分的函数, 对离散函数(不知道解析表达式的波函数进行求解)
include("Potential.jl")
include("Physical_quantity.jl")
# include("")
# include("")
# include("")



end