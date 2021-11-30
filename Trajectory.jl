module Trajectory

export Movement

using ..TDQMC
using ..TDQMC.Discrete_Func_diff               #使用此文件中定义的数值微分的函数
using NumericalAnalysis: Polynomial.Lagrange   #然后对不在格点上的轨迹进行插值求得其导数等等


function Slater_Determination(P::Parameter, Dy::Dynamics)                #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标,对应Electron_num维度的电子波函数




end


function Movement(P::Parameter, Dy::Dynamics)                     # 这里我们使用欧拉法即可
    local 







end