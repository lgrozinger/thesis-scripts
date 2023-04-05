using Pyolin
using StatsBase
using StatsPlots
using LaTeXStrings

opts = Dict(
    :fontfamily => "Computer Modern",
    :guidefontsize => 10,
    :tickfontsize => 10,
    :legendfontsize => 10,
    :linewidth => 3,
    :size => (327, 226)
)

function lmra(strain, backbone)
    ins = RPUExperiment(Pyolin.search(strain, backbone, "1818"))
    outs = RPUExperiment(Pyolin.search(strain, backbone, "Lmra_n1"))
    return Hill(ins, outs)
end


function comparison_plot()
    plt = plot(;xlabel=L"R", ylabel=L"P_{R}(R)", opts..., legendfontsize=8)
    g1 = lmra("KT2440", "pSeva221")
    g2 = lmra("KT2440", "pSeva231")
    g3 = lmra("KT2440", "pSeva251")
    g4 = lmra("DH5alpha", "pAN")
    g5 = lmra("DH5alpha", "pSeva221")
    g6 = lmra("CC118Lpir", "pSeva221")
    g7 = lmra("CC118Lpir", "pSeva231")

    plot!(plt, 0.0:0.01:2.5, g1, label="KT2440 pSeva221", linewidth=2, color=1)
    # scatter!(plt, median.(g1.ins), median.(g1.outs), color=1, label=nothing, markersize=2)
    plot!(plt, 0.0:0.01:2.5, g2, label="KT2440 pSeva231", linewidth=2, color=2)
    # scatter!(plt, median.(g2.ins), median.(g2.outs), color=2, label=nothing, markersize=2)
    plot!(plt, 0.0:0.01:2.5, g3, label="KT2440 pSeva251", linewidth=2, color=3)
    # scatter!(plt, median.(g3.ins), median.(g3.outs), color=3, label=nothing, markersize=2)
    # plot!(plt, 0.0:0.01:2.5, g4, label=L"DH5$\alpha$ pAN", linewidth=2, color=4)
    # scatter!(plt, median.(g4.ins), median.(g4.outs), color=4, label=nothing, markersize=2)
    # plot!(plt, 0.0:0.01:2.5, g5, label=L"DH5$\alpha$ pSeva221", linewidth=2, color=5)
    # scatter!(plt, median.(g5.ins), median.(g5.outs), color=5, label=nothing, markersize=2)
    plot!(plt, 0.0:0.01:2.5, g6, label="CC118Lpir pSeva221", linewidth=2, color=6)
    # scatter!(plt, median.(g6.ins), median.(g6.outs), color=6, label=nothing, markersize=2)
    plot!(plt, 0.0:0.01:2.5, g7, label="CC118Lpir pSeva231", linewidth=2, color=7)
    # scatter!(plt, median.(g7.ins), median.(g7.outs), color=7, label=nothing, markersize=2)
    plot!(plt, ylims=(0, 3))
    return plt
end

# g = lmra("KT2440", "pSeva221")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("KT2440", "pSeva231")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("KT2440", "pSeva251")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("DH5alpha", "pSeva221")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("DH5alpha", "pAN")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("CC118Lpir", "pSeva221")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)

# g = lmra("CC118Lpir", "pSeva231")
# @show g.strain, g.backbone, Pyolin.y0(g), Pyolin.y1(g), Pyolin.K(g), Pyolin.N(g)
