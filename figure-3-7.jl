using Pyolin
using CSV
using Plots

opts = Dict(
    :fontfamily => "Computer Modern",
    :guidefontsize => 10,
    :tickfontsize => 10,
    :legendfontsize => 10,
    :linewidth => 2,
    :size => (500, 500),
    :showaxes => false,
)

function lmra(strain, backbone)
    ins = RPUExperiment(Pyolin.search(strain, backbone, "1818"))
    outs = RPUExperiment(Pyolin.search(strain, backbone, "Lmra_n1"))
    return Hill(ins, outs)
end

p1 = lmra("KT2440", "pSeva221")
p2 = lmra("KT2440", "pSeva231")
p3 = lmra("KT2440", "pSeva251")
p4 = lmra("DH5alpha", "pAN")
p5 = lmra("DH5alpha", "pSeva221")
p6 = lmra("CC118Lpir", "pSeva221")
p7 = lmra("CC118Lpir", "pSeva231")

function poly_curve(gate)
    return hcat(median.(gate.ins), median.(gate.outs))'
end

function frechet_distance(gateA, gateB)
    return Pyolin.frechet(poly_curve(gateA), poly_curve(gateB))
end

function frechet_matrix(gates)
    N = length(gates)
    M = Matrix{Float64}(undef, N, N)
    for j in eachindex(gates)
        for i in eachindex(gates)
            M[i, j] = frechet_distance(gates[j], gates[i])
        end
    end
    return M
end

function single_frechet_map(gates)
    zs = frechet_matrix(gates)
    zs = zs ./ maximum(zs)
    n = length(gates)

    xs = Vector{String}(undef, n)
    for i in 1:n
        xs[i] = gates[i].strain * "\n" * gates[i].backbone
    end
    plt = heatmap(collect(0:n), collect(0:n), zs; opts...)
    plot!(plt; xticks=(collect(0:n) .+ 0.5, xs), yticks=(collect(0:n) .+ 0.5, xs))
    plot!(plt; cbar=false, showaxes=false)
    plot!(plt; xlims=(0, n), ylims=(0, n))
    plot!(plt; aspectratio=1, xrotation=90)

    maketxt(m) = text(round(m; digits=2), "Computer Modern"; pointsize=9, color=(m < maximum(zs) / 2 ? :white : :black))

    anots = maketxt.(vec(zs))
    scatter!(
        plt,
        repeat(1:n, inner=n) .- 0.5, repeat(1:n, outer=n) .- 0.5;
        markerstrokecolor = RGBA(0, 0, 0, 0.0),
        seriesalpha = 0.0,
        label="",
        series_annotations = anots,
        opts...
    )
    return plt
end
