using Pyolin
using StatsBase
using StatsPlots
using LaTeXStrings
using LsqFit

opts = Dict(
    :fontfamily => "Computer Modern",
    :guidefontsize => 10,
    :tickfontsize => 10,
    :legendfontsize => 10,
    :linewidth => 2,
    :size => (327, 226)
)

ins = RPUExperiment(Pyolin.search("KT2440", "pSeva221", "1818"))
outs = RPUExperiment(Pyolin.search("KT2440", "pSeva221", "Lmra_n1"))

function fitting_plot(response)
    plt = plot(;xlabel=L"R", ylabel=L"f_{YFP}(R)", opts...)
    scatter!(plt, response.ins, response.outs; label="Data points")
    xs = 0.0:0.01:1.0
    plot!(plt, xs, response; label="Fitted hill", linewidth=2)
    σs = confidence_interval(response.fit, 0.11)
    @show stderror(response.fit)
    @show margin_error(response.fit, 0.11)
    @show estimate_covar(response.fit)
    @show response.fit.param
    return plt
end

aplt = plot(;xlabel="YFP (RPU)", ylabel="Probability density", opts..., legendfontsize=6)
for (i, iptg) in enumerate(unique(Pyolin.index.iptg))
    StatsPlots.density!(aplt, events(ins[i]), label="$(iptg)" * L"\mu M", color=i)
    vline!([median(ins[i])], color=i, linestyle=:dash, label="", alpha=0.7)
end
plot!(aplt, xlims=(0, 1.5))

bplt = plot(;xlabel="YFP (RPU)", ylabel="Probability density", opts..., legendfontsize=6)
for (i, iptg) in enumerate(unique(Pyolin.index.iptg))
    StatsPlots.density!(bplt, events(outs[i]), label="$(iptg)" * L"\mu M", color=i)
    vline!([median(outs[i])], color=i, linestyle=:dash, label="", alpha=0.7)
end
plot!(bplt, xlims=(0, 4))

cplt = plot(;xlabel="IPTG (μM)", ylabel="YFP (RPU)", opts...)
plot!(cplt, getproperty.(ins, :iptg), median.(ins); label=L"$R_i$")
scatter!(cplt, getproperty.(ins, :iptg), median.(ins); label=false, color=1)
plot!(cplt, getproperty.(outs, :iptg), median.(outs); label=L"$f_{YFP}$", color=2)
scatter!(cplt, getproperty.(outs, :iptg), median.(outs); label=false, color=2)

dplt = fitting_plot(Hill(ins, outs))

