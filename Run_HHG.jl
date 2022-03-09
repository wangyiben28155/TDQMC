module Run_HHG

include("Main_module.jl")
using .TDQMC

using DataFrames, CSV

using Dates

tic = now()

df = CSV.read("Ground_Guide_Wave.csv", DataFrame)
tr = Matrix(CSV.read("Ground_Trajectory.csv", DataFrame))

Electron_num = size(tr, 1)
Ensemble_num = size(tr, 2)
Total_num = Electron_num * Ensemble_num

initial_range = df[:, Total_num+1]                                                                     #最后一列保存的是空间的离散信息
Raw_Array = reshape(complex(Matrix(@. abs(parse(Complex{Float64}, df[:, 1:Total_num])))), (3001, Electron_num, Ensemble_num))  # 500代表的是系综*电子数的数目,因为是单电子所以不需要reshape

s_range = initial_range[end]
s_num = length(initial_range)

l_range = 200.0
l_num = extend_num(s_range, l_range, s_num)

P = Parameter{Float64,Int64}(space_N = l_num, scope = l_range, Δt = 0.05, step_t = 10000)

Raw_GuideWave = stitch_Matrix(Raw_Array, l_num)

Dy = Dynamics{Float64,Int64}(Trajectory = deepcopy(tr), Guide_Wave = deepcopy(Raw_GuideWave), Displace = zeros(Float64, (P.step_t + 1, Ensemble_num, Electron_num)),
    Time = zeros(typeof(P.Δt), Ensemble_num))

df = nothing
Raw_DuideWave = nothing
Raw_Array = nothing
tr = nothing
GC.gc()

parallel_Evolution!(P, Dy)

toc = now()

println("$(convert(DateTime, toc-tic))")

end