##=================####==============##
## SCRIPT TO PRODUCE FIGURES
##=================####==============##

##=================####==============##

# using CSV
using Plots          # plotting
using LaTeXStrings   # latex labels in plots
using DelimitedFiles # latex labels in plots

##=================####==============##
# PATH TO SAVE IMGS (MODIFY WITH YOUR OWN!)

top_path = joinpath(homedir(), "Code", "Collaborations", "networkEpidemiology_MX", "julia_ABM")

data_path = joinpath(top_path, "data", "EOD")
img_out_path = joinpath(top_path, "figures", "EOD")

##=================####==============##
time_lockdown = 15
HC_capacity = 0.25

files = readdir(data_path)
pop_files = filter(x->occursin("populations", x), files)
time_files = filter(x->occursin("infTime", x), files)

cop_ag_vals = map(x->parse(Float64, match(r"_CA_(\d\.\d)_", x).captures[1]), pop_files)

num_agents = parse(Int, match(r"_N_(\d+)_", pop_files[1]).captures[1])
attack_rate = parse(Float64, match(r"_AR_(\d+\.\d+)_", pop_files[1]).captures[1])

##=================####==============##
# FRACTION OF INFECTED AGENTS AS A FUNTCION OF FRACTION OF LOCKDOWN AGENTS

plt = plot(
xlabel = L"t",
ylabel = L"N_I/N",
guidefontsize = 16,
tickfontsize = 10,
colorbar = :none,
legendtitle = L"p_C",
xlims = (0, 70)
)

max_I = Dict()

# for file in pop_files
for i in 1:length(pop_files)
    # data = readdlm(joinpath(data_path, file), ',')
    # cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", file).captures[1])
    # println(file, "|", cop_val)

    data = readdlm(joinpath(data_path, pop_files[i]), ',')
    cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", pop_files[i]).captures[1])
    println(pop_files[i], "|", cop_val)

    t_vals = data[2:end, 1]
    num_I = data[2:end, 3]

    max_I[cop_val] = maximum(num_I)

    plot!(plt,
        t_vals,
        num_I,
        line_z = cop_ag_vals[i],
        lab = cop_ag_vals[i],
        # line_z = cop_val,
        # lab = cop_val,
        lw = 1.25,
        color = :viridis,
        # ylims = (0.0, 0.2)
        )

end

annotate!(plt, 15, 0.61, text("lockdown", :red, :right, 10))
vline!(plt, [time_lockdown], ls = :dash, color = :red, lab = "")

# annotate!(plt, 50, 0.3, text("HC capacity", :blue, :left, 10))
# hline!(plt, [HC_capacity], ls = :dash, color = :blue, lab = "")

display(plt)

##=================####==============##

savefig(plt, joinpath(img_out_path, "SIR_N_$(num_agents)_AR_$(attack_rate)_CA.png"))

##=================####==============##
# FRACTION OF INFECTED AGENTS AS A FUNTCION OF FRACTION OF LOCKDOWN AGENTS

plt = plot(
xlabel = L"t",
ylabel = L"N_S/N",
guidefontsize = 16,
tickfontsize = 10,
colorbar = :none,
legendtitle = L"p_C",
xlims = (0, 70)
)

max_S = Dict()

# for file in pop_files
for i in 1:length(pop_files)
    # data = readdlm(joinpath(data_path, file), ',')
    # cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", file).captures[1])
    # println(file, "|", cop_val)

    data = readdlm(joinpath(data_path, pop_files[i]), ',')
    cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", pop_files[i]).captures[1])
    println(pop_files[i], "|", cop_val)

    t_vals = data[2:end, 1]
    num_I = data[2:end, 2]

    max_S[cop_val] = num_I[end]

    plot!(plt,
        t_vals,
        num_I,
        line_z = cop_ag_vals[i],
        lab = cop_ag_vals[i],
        # line_z = cop_val,
        # lab = cop_val,
        lw = 1.25,
        color = :viridis,
        # ylims = (0.0, 0.2)
        )

end

annotate!(plt, 15, 0.61, text("lockdown", :red, :right, 10))
vline!(plt, [time_lockdown], ls = :dash, color = :red, lab = "")

# annotate!(plt, 50, 0.3, text("HC capacity", :blue, :left, 10))
# hline!(plt, [HC_capacity], ls = :dash, color = :blue, lab = "")

display(plt)

##=================####==============##

savefig(plt, joinpath(img_out_path, "S_N_$(num_agents)_AR_$(attack_rate)_CA.png"))

##=================####==============##
# FRACTION OF INFECTED AGENTS AS A FUNTCION OF FRACTION OF LOCKDOWN AGENTS

plt = plot(
xlabel = L"\tau_I",
ylabel = L"N(\tau_I)",
guidefontsize = 16,
tickfontsize = 10,
colorbar = :none,
legendtitle = L"p_C",
xlims = (0, 70)
)

max_S = Dict()

# for file in pop_files
for i in 1:length(time_files)
    # data = readdlm(joinpath(data_path, file), ',')
    # cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", file).captures[1])
    # println(file, "|", cop_val)

    data = readdlm(joinpath(data_path, time_files[i]), ',')
    cop_val = parse(Float64, match(r"_CA_(\d\.\d)_", time_files[i]).captures[1])
    println(time_files[i], "|", cop_val)

    t_vals = data[2:end, 1]
    num_I = filter(x -> x > 0, data[2:end, 2])

    max_S[cop_val] = num_I[end]

    histogram!(plt,
        num_I,
        # line_z = cop_ag_vals[i],
        # lab = cop_ag_vals[i],
        fill_z = cop_val,
        lab = cop_val,
        # lw = 1.25,s
        color = :viridis,
        # normalized = true
        xlims = (15, 50)
        )

end

display(plt)

##=================####==============##

savefig(plt, joinpath(img_out_path, "infT_$(num_agents)_AR_$(attack_rate)_CA.png"))

##=================####==============##

# FRACTION OF INFECTED AGENTS AS A FUNTCION OF FRACTION OF LOCKDOWN AGENTS

plt = plot(
xlabel = L"p_C",
ylabel = L"\max(N_\star/N)",
guidefontsize = 16,
tickfontsize = 10,
leg = :bottomright
)

pC = sort(collect(keys(max_I)))

plot!(plt,
    pC,
    [max_I[x] for x in pC],
    m = :o,
    label = L"I"
    )

# plot!(plt,
#     pC,
#     [max_S[x] for x in pC],
#     m = :o,
#     label = L"S"
#     )

# annotate!(plt, 15, 0.61, text("lockdown", :red, :right, 10))
# vline!(plt, [time_lockdown], ls = :dash, color = :red, lab = "")

# annotate!(plt, 50, 0.3, text("HC capacity", :blue, :left, 10))
# hline!(plt, [HC_capacity], ls = :dash, color = :blue, lab = "")

display(plt)

##=================####==============##

savefig(plt, joinpath(img_out_path, "poblaciones_asint_N_$(num_agents)_AR_$(attack_rate)_CA.png"))

##=================####==============##
