module Parallelize                             #这个模块是用来并行化的

export parallel_Evolution!, parallel_CTE!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Evolution
using ..TDQMC.Ground_state
using ..TDQMC.Quantity

using Base.Threads, SparseArrays, CSV, DataFrames


function record_Ground(Dy::Dynamics)           #这个函数虽然和function_1里的函数同名但是作用域是隔离的
    local df = DataFrame(Dy.Trajectory, :auto)

    CSV.write("Ground_Trajectory.csv", df)
end


# function Construct(P::Parameter)
#     local λ = 2im * (P.Twice_Δx² / real(P.Δt))      #这个为虚数
#     local λ₀ = 4 * (P.Twice_Δx² / -imag(P.Δt))       #这个为实数
#     local Constructure = ones(typeof(λ), P.space_N - 1)
#     local Section_A::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ₀ + 2.0) .- P.Twice_Δx² .* V_ne.(P.sampling))
#     local Section_B::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2.0) .- P.Twice_Δx² .* V_ne.(P.sampling))
#     local Section_C::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ₀ - 2.0) .+ P.Twice_Δx² .* V_ne.(P.sampling))
#     local Section_D::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2.0) .+ P.Twice_Δx² .* V_ne.(P.sampling))


#     local later_op::SparseMatrixCSC = spzeros(typeof(λ), P.space_N, P.space_N)
#     local former_op::SparseMatrixCSC = spzeros(typeof(λ), P.space_N, P.space_N)

#     later_op = Section_A * Section_B * Section_A
#     former_op = Section_C * Section_D * Section_C

#     return later_op, former_op

# end


function parallel_Evolution!(P::Parameter, Dy::Dynamics)
    local Threads_num::Integer = nthreads()
    local Thread_workload::Vector{<:Integer} = zeros(typeof(Threads_num), Threads_num)         #注意这里有数据竞争, 等下需要修改使用锁或者原子操作

    local λ = 2im * P.Twice_Δx² / P.Δt
    local Constructure = ones(typeof(λ), P.space_N - 1)
    local later_fix::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2.0) .- P.Twice_Δx² .* V_ne.(P.sampling))
    local former_fix::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2.0) .+ P.Twice_Δx² .* V_ne.(P.sampling))

    #这里对固定的构造参数使用传参, 而不在函数原来本身内部定义, 避免储存空间的浪费, 这些改动是在所有函数写好之后进行的更改
    @threads for i = 1:P.Group
        CN_Evolution!(P, Dy, i, later_fix = later_fix, former_fix = former_fix)
        Thread_workload[threadid()] += 1
        println(Thread_workload)
        Dy.Energy[i] = Group_Energy(P, Dy, i)
    end

    println("Caiculation is over!")

end

function parallel_CTE!(P::Parameter, Dy::Dynamics)
    local Threads_num::Integer = nthreads()
    local Thread_workload::Vector{<:Integer} = zeros(typeof(Threads_num), Threads_num)

    local λ = 2im * P.Twice_Δx² / P.Δt
    local Constructure = ones(typeof(λ), P.space_N - 1)
    local later_fix::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2.0) .- P.Twice_Δx² .* V_ne.(P.sampling))
    local former_fix::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2.0) .+ P.Twice_Δx² .* V_ne.(P.sampling))
    # local later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC = Construct(P)

    @threads for i = 1:P.Group
        CT_Evolution!(P, Dy, i, later_fix = later_fix, former_fix = former_fix)
        Dy.Energy[i] = Group_Energy(P, Dy, i)
        Thread_workload[threadid()] += 1
        println(Thread_workload)
    end
    record_Ground(Dy)
    println("Caiculation is over!")

end


end