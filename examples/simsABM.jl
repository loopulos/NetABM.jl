##=================####==============##
## TEST SIMULATION ABM NETWORK-EPIDEMICS
##=================####==============##

##=================####==============##

using StatsBase     # For sample function
using SparseArrays  # For Adyacency Matrix
using Distributions # For probability distributions

using Plots         # plotting
using LaTeXStrings  # latex labels in plots

##=================####==============##

include("NetABM.jl") # Include file with functions and definitions

##=================####==============##
# PATH TO SAVE IMGS (MODIFY WITH YOUR OWN!)

img_out_path = joinpath(homedir(), "Code", "Collaborations", "networkEpidemiology_MX", "julia_ABM", "figures")

##=================####==============##
## SYSTEM'S PARAMETERS

params = Params()

params.p_infected_t0 = 0.01 # Initial fraction of infected agents
params.attack_rate   = 0.25 # Virus attack rate
params.repetition    = 1 # Repetition ID for ensemble statistics

params.μ_cop_agents  = 0.15 # Mean of cooperation probability
params.σ_cop_agents  = 0.1 # Std deviation of cooperation probability
params.μ_recovT_agents = 7 # Mean recovery time

# Link probability (Erdos-Renyi)
params.p_link = ((1+5exp10(-2))*log(params.num_agents))/params.num_agents

##=================####==============##
# INITIALIZE ADJACENCY MATRIX (ERDOS-RENYI)

adj_mat, rows = initialize_erdos_renyi(params);

##=================####==============##
# INITIALIZE MEETS, COOPERATION AND RECOVERY TIMES DISTRIBUTION

rt_μ = params.μ_recovT_agents # MEAN RECOVERY TIME
pc_μ = params.μ_cop_agents # MEAN COOPERATION PROBABILITY
pc_σ = params.σ_cop_agents # STD_DEV COOPERATION PROBABILITY

# meets_dist = Exponential()
meets_dist  = Geometric()
coop_dist   = truncated(Normal(pc_μ, pc_σ), 0.0, 1.0)
recovt_dist = truncated(Poisson(rt_μ), 7, 35)

##=================####==============##

# DEFINE AGENTS POPULATION (DEFAULT VALUES)
agents = [Agent(id) for id in 1:params.num_agents]

# MODIFY AGENT'S CHARACTERISTICS ACCORDING TO ITS DEMOGRAPHICS
initialize_demographics!(agents, params, coop_dist, meets_dist, recovt_dist)

# ASSGINS CONTACTS TO agent FROM adj_mat
for ag in agents
    assign_contacts!(ag, agents, adj_mat, rows)
end

##=================####==============##

time_lockdown = 10

tau = []
num_S = []
num_I = []
num_R = []

# INITIALIZE TIME AND POPULATIONS
params.now_t = 0
populations = get_populations(agents, params)

while populations[2] > 0.0 # WHILE THERE ARE INFECTED AGENTS

    # RECORD DATA
    push!(tau, params.now_t)
    push!(num_S, populations[1])
    push!(num_I, populations[2])
    push!(num_R, populations[3])

    if params.now_t == time_lockdown
        set_coop_agents!(agents, params)
        println("start lockdown")
    end

   # FINDS NEXT STATE AND UPDATES ALL AGENTS
    update_all_agents!(agents, params)
    
    # GET CURRENT POPULATIONS AND UPDATE TIME
    global populations = get_populations(agents, params);
    params.now_t += 1
end

println("END")

##=================####==============##
# POPULATIONS GRAPH

plt = plot(
        xlabel = L"t",
        ylabel = L"N(s)/N",
        guidefontsize = 16,
        tickfontsize = 10
    )

plot!(plt, tau, num_S, color = :blue,  lab = "S", lw = 1.25)
plot!(plt, tau, num_I, color = :red,   lab = "I", lw = 1.25)
plot!(plt, tau, num_R, color = :green, lab = "R", lw = 1.25)
# plot!(plt, tau, 1.0 .- num_R, color = :green, lab = "R", lw = 1.25)

vline!(plt, [time_lockdown], ls = :dash, color = :red, lab = "")

display(plt)

##=================####==============##
savefig(plt, joinpath(img_out_path, "SIR_N_$(params.num_agents)_AR_$(params.attack_rate).png"))
##=================####==============##
# Distribución de Números de Encuentros

plt = plot(
        xlabel = L"N_m",
        ylabel = L"N (N_m)",
        title = "Encounters per agent per time step distribution",
        guidefontsize = 16,
        tickfontsize = 10
    )

vals = 1 .+ rand(meets_dist, params.num_agents)
histogram!(plt, vals, lab = false, bar_width = 0.85)

# vals = rand(meets_dist, params.num_agents)
# histogram!(plt, vals, normalized = true, lab = false, bar_width = 0.85)

xticks!(plt, 1:maximum(vals))
display(plt)

##=================####==============##
savefig(plt, joinpath(img_out_path, "meetsDist_N_$(params.num_agents).png"))
##=================####==============##
# Distribución de Tiempo de Recuperación

plt = plot(
        xlabel = L"R_t",
        ylabel = L"N(R_t)",
        title = "Agent's recovery times distribution",
        guidefontsize = 16,
        tickfontsize = 10
    )

vals = rand(recovt_dist, params.num_agents)
histogram!(plt, vals, lab = false)

# histogram!(plt, vals, normalized = true, lab = false)

# xticks!(plt, unique(vals))
display(plt)

##=================####==============##
savefig(plt, joinpath(img_out_path, "recovtDist_N_$(params.num_agents).png"))
##=================####==============##
# Tiempo de Infección

plt = plot(
        xlabel = L"t_I",
        ylabel = L"p(t_I)",
        title = "Time of infection distribution",
        guidefontsize = 16,
        tickfontsize = 10
    )

histogram!(plt, [ag.infection_t for ag in agents], normalized = true, lab = false)
display(plt)

##=================####==============##
savefig(plt, joinpath(img_out_path, "tiempoI_N_$(params.num_agents)_AR_$(params.attack_rate).png"))
##=================####==============##
# Distribución de Grado

plt = plot(
        xlabel = L"\kappa_i",
        ylabel = L"p\;(\kappa_i)",
        title = "Network's degree distribution ",
        guidefontsize = 16,
        tickfontsize = 10
    )

histogram!(plt, [ag.degree_t for ag in agents], normalized = true, lab = false)
display(plt)

##=================####==============##
savefig(plt, joinpath(img_out_path, "degDist_N_$(params.num_agents)_ER.png"))
##=================####==============##
