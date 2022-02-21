module Trajectory

export Movement!

using ..TDQMC
using ..TDQMC.Find_nearest_k

using Interpolations, LinearAlgebra                                               #然后对不在格点上的轨迹进行插值求得其导数等等

@inline function thermalization(t::Complex; cut_off = 15.0, amp = 1.0)            #衰减振幅包络
    return real(t) <= cut_off ?  amp * (1.0 - real(t) / cut_off) : 0.0
end


function find_lattice_wave!(localtion::AbstractFloat, serial_num::Integer,
    k_itp_Wave_Vector::Vector{<:Vector{<:Complex}}, P::Parameter, Dy::Dynamics;
    k::Integer, In_num::Int64)
    #得到某一个粒子的位置坐标
    local Index = Dy.Index[serial_num]

    local indexVec_k::UnitRange{<:Integer} = find_k_index(localtion, x = P.sampling, k = k)
    #选取五个点用来插值的波函数矢量

    for i in 1:In_num
        k_itp_Wave_Vector[i] = Dy.Guide_Wave[Index[i], serial_num][indexVec_k]
    end

    return P.sampling[indexVec_k]           #返回Particle_num数的电子轨迹位置处,每个组不同电子的波函数的五点取样的波函数值
end


function Interpolation_Wave!(Particle_num::Integer, serial_num::Integer,
    Vec_Wave::Vector{<:Complex}, Vec_Derivate::Vector{<:Complex}, P::Parameter, Dy::Dynamics;
    k::Integer = 5, Type::DataType, In_num::Int64)

    local localtion::AbstractFloat = Dy.Trajectory[Particle_num, serial_num]                                                                   #小变量没必要传参,但读取需要时间,为方便下面使用,就使用变量读取保存
    local yd::Vector{<:Vector{<:Complex}} = [zeros(Type, k) for i = 1:In_num]
    local xd::LinRange = find_lattice_wave!(localtion, serial_num, yd, P, Dy, k = k, In_num = In_num)

    local interp_cubic::Vector{<:Interpolations.Extrapolation} = Vector{Interpolations.Extrapolation}(undef, In_num)

    for i in 1:In_num
        interp_cubic[i] = CubicSplineInterpolation(xd, yd[i])        #得到三次样条插值函数, i 代表每个电子的波函数的五点取样值 

        Vec_Wave[i] = interp_cubic[i](localtion)
        Vec_Derivate[i] = Interpolations.gradient(interp_cubic[i], localtion)[1]
    end

end


function Slater_determinant!(Symmetric::Matrix{<:Complex}, Derivate::Matrix{<:Complex},
    P::Parameter, Dy::Dynamics, serial_num::Integer; Type::DataType, In_num::Int64)   #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标
    local Index = Dy.Index[serial_num]

    local Vec_Wave::Vector{<:Complex}, Vec_Derivate::Vector{<:Complex} = (zeros(Type, In_num), zeros(Type, In_num))

    for i in 1:In_num                                                             #这里的for 循环只对边界内电子的索引进行,上面的循环都要这样
        Interpolation_Wave!(Index[i], serial_num, Vec_Wave, Vec_Derivate, P, Dy, Type = Type, In_num = In_num)             #这里应该改为在边界内对应的电子的行列式进行计算
    
        Symmetric[:, i] = Vec_Wave                                        #在Dy中加入一个Inbound来记录在边界内的索引
        Derivate[:, i] = Vec_Derivate
    end

end


function Velocity(P::Parameter, Dy::Dynamics, serial_num::Integer;
    Type::DataType = eltype(eltype(Dy.Guide_Wave)), In_num::Int)

    local symmetric_determinate, Derivate_eachcoodinate = (zeros(Type, (In_num, In_num)), zeros(Type, (In_num, In_num)))
    local Derivate_WaveFunc::Vector{<:Matrix{<:Complex}} = [deepcopy(symmetric_determinate) for i = 1:In_num]
    local symmetric_WaveFunc::Vector{<:Matrix{<:Complex}} = fill(symmetric_determinate, In_num)

    local Vector_velocity::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), In_num)


    Slater_determinant!(symmetric_determinate, Derivate_eachcoodinate, P, Dy, serial_num, Type = Type, In_num = In_num)

    Derivate_WaveFunc .= [deepcopy(symmetric_determinate) for i = 1:In_num]
    symmetric_WaveFunc .= fill(symmetric_determinate, In_num)

    for i in 1:In_num
        Derivate_WaveFunc[i][:, i] = Derivate_eachcoodinate[:, i]
    end

    Vector_velocity[:] = @. imag(det(Derivate_WaveFunc) / det(symmetric_WaveFunc))

    return Vector_velocity
end


function Movement!(P::Parameter, Dy::Dynamics, serial_num::Integer, Vec_Trajectory::SubArray;
    dt = P.Δt, In_num::Int64 = Dy.In_num[serial_num])                     # 这里我们使用欧拉法即可

    local Index = Dy.Index[serial_num]
    local OutBoundary_index::Vector{<:Integer} = Int64[]

    local Total_time = P.step_t * real(dt)


    Vec_Trajectory[Index] .+= real(dt) * Velocity(P, Dy, serial_num, In_num = In_num)

    if imag(dt) != 0.0                          #这一部分是为了在寻找基态的过程中避免quantum dift,使用线性衰减的随机数
        Vec_Trajectory[Index] .+= 0.1 * thermalization(Dy.Time[serial_num], cut_off = Total_time) * (rand(In_num) .- 0.5)
    end

    OutBoundary_index = findall(x -> abs(x) > P.scope, Vec_Trajectory)       #上面已经对粒子的位置进行了更新,
    #下面的部分是为了将超过边界的粒子停留在边界

    if isempty(OutBoundary_index)
        return nothing
    else
        Vec_Trajectory[OutBoundary_index] = sign.(Vec_Trajectory[OutBoundary_index]) * P.scope           #要改成Inf,等上面完工再说
    end
end



end