module visual

export plot_probability, plot_Ground

using ..TDQMC

using PyPlot, CSV, DataFrames

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



function plot_HHG()
    local df = CSV.read("High_Harmonic_Generation.csv", DataFrame)
    local L = length(df.t)
    local xlocal = collect(0:T:16T)
    local Harmonic_order = 1 / (df.t[end] - df.t[1]) .* (LinRange(0, L - 1, L))

    subplot(2, 1, 1)
    plot(df.t, df.a_t, color = "maroon")
    xticks(xlocal, 1:length(xlocal))
    ylabel("a(t)")
    xlabel("cycles")
    title("The a(t) of the Wave packet")
    grid()
    subplot(2, 1, 2)
    semilogy(Harmonic_order, 2 / L * abs2.(fft(df.a_t)), color = "navy")
    ylabel(L"$log_{10}(a(\omega)^2)$")
    grid()
    xticks(ω_0 / (2pi) * collect(0:2:20), collect(0:2:20))
    xlabel("Hamonic order")
    xlim(0, 20 * ω_0 / (2pi))
    show()
end




end