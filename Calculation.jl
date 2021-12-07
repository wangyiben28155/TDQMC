module Parallelize                             #这个模块是用来并行化的

export parallel_Evolution!

using ..TDQMC
using ..TDQMC.Evolution
using Base.Threads


function parallel_Evolution!(P::Parameter, Dy::Dynamics)
    local Threads_num::Integer = nthreads()
    local Thread_workload::Vector{Atomic{<:Integer}} = [Atomic{Int}(0) for i = 1:Threads_num]         #注意这里有数据竞争, 等下需要修改使用锁或者原子操作

    @threads for i = 1:P.Group
        CN_Evolution!(P, Dy, i)
        @inbounds atomic_add!(Thread_workload[threadid()], 1)
    end

    println(Thread_workload)
    
end


end