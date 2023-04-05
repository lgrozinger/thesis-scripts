using Pyolin
using StatsPlots

opts = Dict(
    :fontfamily => "Computer Modern",
    :guidefontsize => 10,
    :tickfontsize => 10,
    :legendfontsize => 10,
    :linewidth => 2,
    :size => (327, 226)
)

experiments = RawExperiment(Pyolin.search("KT2440", "pSeva221", "1717"))
rawplt = plot(;xlabel="YFP (A.U.)", ylabel="Frequency Density", opts...)

color = 1
for i in eachindex(experiments)
    if experiments[i].iptg ∈ [0, 5, 10, 100, 500, 1000]
        density!(rawplt, events(experiments[i]); label="$(experiments[i].iptg)μM", color=color, opts...)
        vline!(rawplt, [median(events(experiments[i]))]; label=false, color=color, linestyle=:dash, linewidth=1)
        global color = color + 1
    end
end

plot!(rawplt, xlims=(0, 30))
