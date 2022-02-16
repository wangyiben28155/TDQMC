module Trajectory

export Movement!

using ..TDQMC
using ..TDQMC.Find_nearest_k

using Interpolations, LinearAlgebra   #然后对不在格点上的轨迹进行插值求得其导数等等

@inline function thermalization(t::Complex; cut_off = 15.0, amp = 1.0)            #衰减振幅包络
    return imag(t) >= cut_off ? 0.0 : amp * (1.0 - imag(t) / cut_off)
end


function find_lattice_wave(Particle_num::Integer, serial_num::Integer, P::Parameter, Dy::Dynamics; k::Integer = 5)    #这个对应的是粒子的轨迹矢量,用来得到相应的波函数插值的部分
    local localtion = Dy.Trajectory[Particle_num, serial_num]                          #得到某一个粒子的位置坐标
    local indexVec_k::UnitRange{<:Integer} = find_k_index(localtion, x = P.sampling, k = k)
    local Type_0 = eltype(eltype(Dy.Guide_Wave))
    local k_itp_Wave_Vector::Vector{<:Vector{<:Complex{<:AbstractFloat}}} = [zeros(Type_0, k) for i in 1:P.electron]     #选取五个点用来插值的波函数矢量,注意这种
    #赋值的方式每一个重复的部分在后面赋值的时#候不会产生关联.

    for i = 1:P.electron
        k_itp_Wave_Vector[i] = Dy.Guide_Wave[i, serial_num][indexVec_k]
    end

    return P.sampling[indexVec_k], k_itp_Wave_Vector                              #返回Particle_num数的电子轨迹位置处,每个组不同电子的波函数的五点取样的波函数值
end


function Interpolation_Wave(Particle_num::Integer, serial_num::Integer, P::Parameter, Dy::Dynamics)
    local localtion = Dy.Trajectory[Particle_num, serial_num]                                           #这里选出的是一个轨迹
    local xd::LinRange, yd::Vector = find_lattice_wave(Particle_num, serial_num, P, Dy)
    local Type_1 = eltype(eltype(yd))
    local Vec_Wave::Vector{<:Complex} = zeros(Type_1, P.electron)
    local Vec_Derivate::Vector{<:Complex} = zeros(Type_1, P.electron)
    local interp_cubic::Interpolations.Extrapolation

    for i = 1:P.electron
        interp_cubic = CubicSplineInterpolation(xd, yd[i])        #得到三次样条插值函数, i 代表每个电子的波函数的五点取样值
        Vec_Wave[i] = interp_cubic(localtion)                     #然后
        Vec_Derivate[i] = Interpolations.gradient(interp_cubic, localtion)[1]
    end
    return Vec_Wave, Vec_Derivate    #返回在particle_num的电子的轨迹处, 对所有电子的波函数和导数的取值
end


function Slater_determinant(P::Parameter, Dy::Dynamics, serial_num::Integer)           #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标,对应Electron_num维度的电子波函数
    local Type_2 = eltype(eltype(Dy.Guide_Wave))
    local Vec_Wave::Vector{<:Complex}, Vec_Derivate::Vector{<:Complex} = (zeros(Type_2, P.electron), zeros(Type_2, P.electron))
    local symmetric_determinate::Matrix{<:Complex} = zeros(Type_2, (P.electron, P.electron))
    local Derivate_eachcoodinate::Matrix{<:Complex} = zeros(Type_2, (P.electron, P.electron))

    for i = 1:P.electron
        Vec_Wave, Vec_Derivate = Interpolation_Wave(i, serial_num, P, Dy)
        symmetric_determinate[:, i] = Vec_Wave
        Derivate_eachcoodinate[:, i] = Vec_Derivate
    end

    return symmetric_determinate, Derivate_eachcoodinate
end


function Velocity(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local symmetric_determinate, Derivate_eachcoodinate = Slater_determinant(P, Dy, serial_num)
    local Derivate_WaveFunc::Vector{<:Matrix{<:Complex}} = [deepcopy(symmetric_determinate) for i = 1:P.electron]
    local symmetric_WaveFunc::Vector{<:Matrix{<:Complex}} = fill(symmetric_determinate, P.electron)
    local Vector_velocity::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), P.electron)

    for i = 1:P.electron
        Derivate_WaveFunc[i][:, i] = Derivate_eachcoodinate[:, i]
    end

    Vector_velocity[:] = @. imag(det(Derivate_WaveFunc) / det(symmetric_WaveFunc))

    return Vector_velocity
end


function Movement!(P::Parameter, Dy::Dynamics, serial_num::Integer; dt = P.Δt)                     # 这里我们使用欧拉法即可
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)
    local OutBoundary_index::Vector{<:Integer} = Int64[]
    local Total_time = P.step_t * imag(dt)

    Vec_Trajectory[:] .+= real(dt) * Velocity(P, Dy, serial_num)

    if imag(dt) != 0                          #这一部分是为了在寻找基态的过程中避免quantum dift,使用线性衰减的随机数
        Vec_Trajectory[:] += 0.1 * thermalization(Dy.Time[serial_num], cut_off = Total_time) * (rand(P.electron) .- 0.5)
    end

    OutBoundary_index = findall(x -> abs(x) > P.scope, Vec_Trajectory)       #下面的部分是为了将超过边界的粒子停留在边界

    if isempty(OutBoundary_index)
        return nothing
    else
        Vec_Trajectory[OutBoundary_index] = sign.(Vec_Trajectory[OutBoundary_index]) * P.scope
    end
end



end