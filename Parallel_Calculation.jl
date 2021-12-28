module Parallelize                             #这个模块是用来并行化的

export parallel_Evolution!

using ..TDQMC
using ..TDQMC.Evolution
using ..TDQMC.Quantity 

using Base.Threads


function parallel_Evolution!(P::Parameter, Dy::Dynamics)
    local Threads_num::Integer = nthreads()
    local Thread_workload::Vector{<:Integer} = zeros(typeof(Threads_num), Threads_num)         #注意这里有数据竞争, 等下需要修改使用锁或者原子操作

    @threads for i = 1:P.Group
        CN_Evolution!(P, Dy, i)
        @inbounds Thread_workload[threadid()] += 1
        println(Thread_workload)
    end

    println("Caiculation is over!")
    
end


end