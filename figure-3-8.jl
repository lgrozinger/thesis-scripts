using Pyolin
using CSV
using Plots
using StatsPlots
using HypothesisTests
using LaTeXStrings


GATES = [
    "Amer_f1", "Amtr_a1", "Beti_e1", "Bm3r1_b2", "Bm3r1_b3",
    "Hiyiir_h1", "Lcara_i1", "Litr_l1", "Lmra_n1", "Phif_p1",
    "Phif_p2", "Psra_r1", "Qacr_q1", "Srpr_s1", "Srpr_s4",
]

CONTEXTS = [
    ("KT2440", "pSeva221"),
    ("KT2440", "pSeva231"),
    ("KT2440", "pSeva251"),
    ("DH5alpha", "pAN"),
    ("DH5alpha", "pSeva221"),
    ("CC118Lpir", "pSeva221"),
    ("CC118Lpir", "pSeva231"),
]

opts = Dict(
    :fontfamily => "Computer Modern",
    :guidefontsize => 10,
    :tickfontsize => 10,
    :legendfontsize => 10,
    :linewidth => 2,
    :size => (327, 226),
    :showaxes => false,
)

function contexts_in_gates(gate)
    results = Vector{Hill}(undef, length(CONTEXTS))
    for (i, context) in enumerate(CONTEXTS)
        ins = RPUExperiment(Pyolin.search(context..., "1818"))
        outs = RPUExperiment(Pyolin.search(context..., gate))
        results[i] = Hill(ins, outs)
    end
    return results
end

function gates_in_context(strain, backbone)
    ins = RPUExperiment(Pyolin.search(strain, backbone, "1818"))
    results = Vector{Hill}(undef, length(GATES))
    for (i, gate) in enumerate(GATES)
        outs = RPUExperiment(Pyolin.search(strain, backbone, gate))
        results[i] = Hill(ins, outs)
    end
    return results
end

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

function frechet_vector(gates)
    N = length(gates)
    distances = Float64[]
    for i in 1:N
        for j in i+1:N
            if gates[i].plasmid == gates[j].plasmid
                push!(distances, frechet_distance(gates[j], gates[i]))
            end
        end
    end
    return distances
end

function frechet_vector(gates, f)
    N = length(gates)
    distances = Float64[]
    for i in 1:N
        for j in i+1:N
            if f(gates[i], gates[j]) && gates[i].plasmid == gates[j].plasmid
                push!(distances, frechet_distance(gates[j], gates[i]))
            end
        end
    end
    return distances
end

# p1 = gates_in_context("KT2440", "pSeva221")
# p2 = gates_in_context("KT2440", "pSeva231")
# p3 = gates_in_context("KT2440", "pSeva251")
# p4 = gates_in_context("DH5alpha", "pAN")
# p5 = gates_in_context("DH5alpha", "pSeva221")
# p6 = gates_in_context("CC118Lpir", "pSeva221")
# p7 = gates_in_context("CC118Lpir", "pSeva231")

# everything = vcat(p1, p2, p3, p4, p5, p6, p7)
# all_distances = frechet_vector(everything)

# same_host_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain == y.strain
# )
# diff_host_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain != y.strain
# )

# putida_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain==y.strain=="KT2440"
# )

# coli_distances = frechet_vector(
#     everything,
#     (x, y) -> (x.strain∈["DH5alpha", "CC118Lpir"]) && (y.strain∈["DH5alpha", "CC118Lpir"])
# )

# putida_coli_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain=="KT2440" && (y.strain ∈ ["DH5alpha", "CC118Lpir"])
# )

# coli_coli_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain=="DH5alpha" && y.strain=="CC118Lpir"
# )

# dh5alpha_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain==y.strain=="DH5alpha",
# )

# cc118_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain==y.strain=="CC118Lpir",
# )

# same_coli_distances = frechet_vector(
#     everything,
#     (x, y) -> x.strain==y.strain && (x.strain ∈ ["DH5alpha", "CC118Lpir"])
# )

# same_backbone_distances = frechet_vector(
#     everything,
#     (x, y) -> x.backbone==y.backbone
# )

# diff_backbone_distances = frechet_vector(
#     everything,
#     (x, y) -> x.backbone!=y.backbone
# )

# function make_density_comparison(A, B, labels)
#     plt = density(A; label=first(labels), fill=(0, 0.5), linewidth=2, opts...)
#     density!(plt, B; label=last(labels), fill=(0, 0.5), linewidth=2, opts...)
#     plot!(plt, xlabel="Fréchet distance", ylabel="Probability density", xlims=(0, 10))
#     vline!(plt, [median(A), median(B)], color=[1, 2], label=false, linestyle=:dash)
#     return plt
# end

# println("NULL HYPOTHESIS: There is no difference in frechet distance whether the hosts are different or not") 
# print(HypothesisTests.MannWhitneyUTest(same_host_distances, diff_host_distances))

# println()
# println("NULL HYPOTHESIS: There is no difference in frechet distance whether the hosts the same bacteria or a mix of E. coli and P. putida") 
# print(HypothesisTests.MannWhitneyUTest(vcat(coli_distances, putida_distances), putida_coli_distances))

# println()
# println("NULL HYPOTHESIS: There is no difference in frechet distance whether the hosts are E. coli or P. putida")
# print(HypothesisTests.MannWhitneyUTest(putida_distances, coli_distances))
# println()
# print(HypothesisTests.ApproximateTwoSampleKSTest(putida_distances, coli_distances))

# println()
# println("NULL HYPOTHESIS: There is no difference in frechet distance whether the hosts are E. coli DH5 or E. coli CC118")
# print(HypothesisTests.MannWhitneyUTest(dh5alpha_distances, cc118_distances))
# println()
# print(HypothesisTests.ApproximateTwoSampleKSTest(dh5alpha_distances, cc118_distances))

# println()
# println("NULL HYPOTHESIS: There is no difference in frechet distance whether the backbones are the same or not")
# print(HypothesisTests.MannWhitneyUTest(same_backbone_distances, diff_backbone_distances))

# putidas = vcat(p1, p2, p3)
# colis = vcat(p4, p5, p6, p7)
# seva221 = vcat(p1, p5, p6)
# seva231 = vcat(p2, p7)
# seva251 = copy(p3)
# pan = copy(p4)

function frechet_map()
    zs = zeros(length(CONTEXTS), length(CONTEXTS))
    for (i, gate) in enumerate(GATES)
        gates = contexts_in_gates(gate)
        zs = zs .+ frechet_matrix(gates)
    end
    zs = zs ./ maximum(zs)

    n = length(CONTEXTS)

    xs = Vector{String}(undef, n)
    for i in 1:n
        xs[i] = first(CONTEXTS[i]) * "\n" * last(CONTEXTS[i])
    end

    plt = heatmap(collect(0:n), collect(0:n), zs; opts...)
    plot!(plt; xticks=(collect(0:n) .+ 0.5, xs), yticks=(collect(0:n) .+ 0.5, xs))
    plot!(plt; cbar=false, showaxes=false)
    plot!(plt; xlims=(0, n), ylims=(0, n))
    plot!(plt; aspectratio=1, xrotation=90)

    maketxt(m) = text(
        round(m; digits=2),
        "Computer Modern";
        pointsize=9,
        color=(m < maximum(zs) / 2 ? :white : :black)
    )
    
    anots = maketxt.(vec(zs))
    scatter!(
        plt,
        repeat(1:n, inner=n) .- 0.5, repeat(1:n, outer=n) .- 0.5;
        markerstrokecolor = RGBA(0, 0, 0, 0.0),
        seriesalpha = 0.0,
        label="",
        series_annotations = anots,
        opts...,
        size = (500, 500),
    )
    return plt
end
