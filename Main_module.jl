module TDQMC

export Dynamics, Parameter, Movement!, CN_Evolution!

using Distributions, Random, SparseArrays                                   #用来初始化初始的波函数和系综粒子分布
import Base.@kwdef

const Ensemble_num, Electron_num = (50, 1)                                 #设定一些计算的参数
const L, x_num, step_t = (350.0, 1751, 13001)
const σ = 30


function initializer_guideWave(x::T) where {T<:AbstractFloat}                                           #这里因为算的是一维的波函数, 所以返回的矩阵为三维矩阵,第一维为波函数,第二维为系综粒子数,第三维为电子数

    return complex(sqrt.(pdf(Normal(x, σ), LinRange(-L, L, x_num))))
end

function generate_distribution(x::Int)
    Random.seed!(x)                                                       #确定下一次使用rand函数的撒点是固定的, 每次调用程序也能保证下面的分布是确定的
    Initial_Distriution = rand(Normal(μ, σ), Ensemble_num)                #将确定的撒点分布储存起来, 用于下面实例化的赋值

    return Initial_Distriution'
end

@kwdef mutable struct Dynamics{T<:AbstractFloat}
    Trajectory::Matrix{T} = vcat(generate_distribution.(1:Electron_num)...)   = 先给定一个初始的轨迹分布, 目前只是计算氢原子的情况,
     所以先只使用向量保存当前时刻的位置信息, 并且之后代不同电子之间粒子的位置之间的运算会有矢量化运算, 所以将其放在列坐标加快运算 =#
    Guide_Wave::Matrix{<:Vector{<:Complex{T}}} = initializer_guideWave.(Trajectory)     = 这里得到的分布是概率密度的开平方作为初始的波函数
    进行演化 =#
    Slater_Determinant::Matrix{T} = zeros(Complex{T}, (Electron_num, Electron_num))
    Slater_Determinant_Dericative::Matrix{T} = zeros(Complex{T}, (Electron_num, Electron_num))
    Energy::Vector{T} = zeros(T, Ensemble_num)
    Time::Union{T,Complex{T}} = 0.0
end


@kwdef struct Parameter{T1<:AbstractFloat,T2<:Integer}                        #用来控制计算参数的
    electron::T2 = Electron_num
    Group::T2 = Ensemble_num
    space_N::T2 = x_num                                                       #划分的格点的总数,后面做离散傅里叶变换的时候会用得到
    scope::T1 = L                                                             #确定波函数的计算范围为-scope到+scope
    Δx::T1 = 2 * scope / (N - 1)                                              #波函数的离散的空间间隔
    Square_Δx::T1 = 2 * (P.Δx)^2
    sampling::LinRange{T1} = LinRange(-scope, scope, space_N)
    Δt::Union{T1,Complex{T1}} = 0.05                                          #划分的时间间隔, 尝试时间迭代区间, 因为考虑到虚时演化, 所以类型设定为复数
    step_t::T2 = step_t
end

include("Numerical_Diff_DiscreteFunc.jl")                        #用来做数值微分的函数, 对离散函数(不知道解析表达式的波函数进行求解)
include("Find_nearest.jl")
include("Potential.jl")
include("Physical_quantity.jl")
include("Crank_Nicolson.jl")
include("Trajectory.jl")
include("visualization.jl")
include("Evolution.jl")
include("Calculation.jl")


using .Evolution
using .Trajectory

end