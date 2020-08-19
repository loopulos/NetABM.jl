##=================####==============##

using NetABM
using LightGraphs

using StatsBase
using Statistics
using Bootstrap

using Plots
using LaTeXStrings

using Distributions

##=================####==============##
# SYSYTEM PARAMETERS

params = Params()

params.repetition = 1      # Repetition ID for ensemble statistics
params.num_agents = 100000   # Number of agents (32768)
params.attack_rate = 0.5   # Virus attack rate
params.p_infected_t0 = 0.001# Initial fraction of infected agents

params.σ_cop_agents = 0.1  # Std deviation of cooperation probability
params.μ_cop_agents = 0.25  # Mean of cooperation probability
params.μ_recovT_agents = 15 # Mean recovery time

##=================####==============##
# INITIALIZE MEETS, COOPERATION AND RECOVERY TIMES DISTRIBUTION

pc_μ = params.μ_cop_agents    # MEAN COOPERATION PROBABILITY
pc_σ = params.σ_cop_agents    # STD_DEV COOPERATION PROBABILITY
rt_μ = params.μ_recovT_agents # MEAN RECOVERY TIME

# meets_dist = Exponential()
meets_dist = Geometric()
recovt_dist = truncated(Poisson(rt_μ), 7, 35)

coop_dist = truncated(Normal(pc_μ, pc_σ), 0.0, 1.0)

init_I = params.p_infected_t0
init_no_I = 1.0 - params.p_infected_t0

##=================####==============##

#  g = erdos_renyi(num_agents,0.4)
g = barabasi_albert(params.num_agents, 9)
agents = [Agent(i) for i = 1:(params.num_agents)];

NetABM.mod_init_demographics!(
    agents,
    coop_dist,
    meets_dist,
    recovt_dist;
    states = ["S", "I"],
    initial = [init_no_I, init_I],
    ϵ_thr = [0.25, 0.75]
)
map(x -> assign_contacts!(g, x), agents);

# TIME AT WHICH LOCKDOWN IS ANNOUNCED
#  time_lockdown = 10

tau = []
num_S = []
num_I = []
num_R = []
t_p_cops = Array{Float64}[]

# INITIALIZE TIME AND POPULATIONS
params.now_t = 0
populations = get_populations(agents, params)
global this_I = populations[2]

#  lockdown = false
#  time_lockdown = 0

push!(tau, params.now_t)
push!(num_S, populations[1])
push!(num_I, populations[2])
push!(num_R, populations[3])
push!(t_p_cops, [ag.p_cop for ag in agents])

println(params.now_t, "|", populations[2])
println("SIM STARTS")

while populations[2] > 0.0 # WHILE THERE ARE INFECTED AGENTS

    # RECORD DATA
    push!(tau, params.now_t)
    push!(num_S, populations[1])
    push!(num_I, populations[2])
    push!(num_R, populations[3])
    push!(t_p_cops, [ag.p_cop for ag in agents])

    # FINDS NEXT STATE AND UPDATES ALL AGENTS
    update_all_agents!(agents, params)

    # OPINION DYNAMICS
    for ag_id = 1:(params.num_agents)
        NetABM.HK_od_model!(agents, ag_id)
    end
    NetABM.set_athome_agents!(agents)

    #  EACH AGENTS DECIDES TO STAY AT HOME/ NOT TO STAY AT HOME
    #  if params.now_t >= time_lockdown
    #  #  if lockdown
    #      NetABM.set_athome_agents!(agents)
    #      # OPINION DYNAMICS
    #      for ag_id = 1:(params.num_agents)
    #          NetABM.HK_od_model!(agents, ag_id)
    #      end
    #  end

    # GET CURRENT POPULATIONS AND UPDATE TIME
    populations = get_populations(agents, params)
    params.now_t += 1

    prev_I = this_I
    this_I = populations[2]

    #  if !lockdown
    #      if this_I / prev_I >= 2.0
    #          lockdown = true
    #          time_lockdown = params.now_t
    #      end
    #  end

    #  println(params.now_t, "|", populations[2], "|", this_I / prev_I)
    println(params.now_t, populations)
end

println("NO MORE INFECTED!")

##=================####==============##

plt = plot()
plot!(plt, tau, num_I ./ params.num_agents, lab = "I")
#  plot!(plt, tau, num_S ./ params.num_agents, lab = "S")
#  plot!(plt, tau, num_R ./ params.num_agents, lab = "R")
#  vline!(plt, [time_lockdown], color = :red, ls = :dash, lab = "")

##=================####==============##

od_t = transpose(hcat(t_p_cops...));
plot(od_t, alpha = 0.45, lab = "")

##=================####==============##
