module visual

export plot_probability, plot_Ground

using ..TDQMC
using PyPlot

function plot_probability(P::Parameter)
    local df = CSV.read("Time_evolution_wave_function.csv", DataFrame)

    PyPlot.axes(yscale = "log")
    plot(df.x, df.Probability, color = "black")
    grid()
    title("E_0 = 0.1, ω=0.148, t=16T")
    ylabel("Probability")
    xticks(collect(LinRange(-P.scope, P.scope, 11)))
    xlabel("X(a.u)")

end

function plot_Ground()
    local df = CSV.read("Ground_state.csv", DataFrame)

    df.wave = @. real(parse(Complex{Float64}, df.wave))

    plot(df.x, df.wave, color = "black")
    grid()
    title("Ground state of wave function")
    xlabel("x")
    ylabel("amplitude")

end






end