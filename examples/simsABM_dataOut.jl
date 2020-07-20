##=================####==============##
## TEST SIMULATION ABM NETWORK-EPIDEMICS
##=================####==============##

##=================####==============##

using ArgParse      # For CLI arguments
using StatsBase     # For sample function
using SparseArrays  # For Adyacency Matrix
using Distributions # For probability distributions
using DelimitedFiles # For data output

##=================####==============##

include("NetABM.jl") # Include file with functions and definitions

### ============== ### ============== ### ============== ###
### SYSTEM'S PARAMETERS FROM COMMAND LINE
### ============== ### ============== ### ============== ###

args_settings = ArgParseSettings()

@add_arg_table! args_settings begin
    "--numAgents", "-n"
        help = "Number of Agents"
        arg_type = Int
        default = 32768
    "--initialI", "-i"
        help = "Initial fraction of Infected Agents"
        arg_type = Float64
        default = 0.01
    "--attackR", "-a"
        help = "Attack Rate"
        arg_type = Float64
        default = 0.25
    "--meanRecT", "-t"
        help = "Mean recovery time"
        arg_type Int64
        default = 7
    "--meanCop", "-c"
        help = "Mean cooperation probability"
        arg_type = Float64
        default = 0.5
    "--stdCop", "-s"
        help = "Std Deviation cooperation probability"
        arg_type = Float64
        default = 0.1
    "--rep", "-r"
        help = "Repetition ID (for ensemble statistics)"
        arg_type = Int
        default = 1
end

args = parse_args(args_settings)

##=================####==============##
## ASIGN SYSTEM'S PARAMETERS FROM CLI

params = Params()

params.num_agents      = args["numAgents"] # Number of agents (32768)
params.p_infected_t0   = args["initialI"] # Initial fraction of infected agents
params.μ_recovT_agents = args["meanRecT"] # Std deviation of cooperation probability
params.μ_cop_agents    = args["meanCop"] # Mean of cooperation probability
params.σ_cop_agents    = args["stdCop"] # Std deviation of cooperation probability
params.attack_rate     = args["attackR"] # Virus attack rate
params.repetition      = args["rep"] # Repetition ID for ensemble statistics

# Link probability (Erdos-Renyi)
params.p_link = ((1+5exp10(-2))*log(params.num_agents))/params.num_agents

##=================####==============##
# TEST PARAMETERS

# params.num_agents    = 2^15 # Number of agents (32768)
# params.p_infected_t0 = 0.01 # Initial fraction of infected agents
# params.attack_rate   = 0.25 # Virus attack rate
# params.repetition    = 1 # Repetition ID for ensemble statistics

# Link probability (Erdos-Renyi)
# params.p_link = ((1+5exp10(-2))*log(params.num_agents))/params.num_agents

##=================####==============##
# PATH TO SAVE DATA (MODIFY WITH YOUR OWN!)

data_out_path = joinpath(homedir(),
    "Code", "Collaborations", "networkEpidemiology_MX", "julia_ABM", "data")

##=================####==============##
# INITIALIZE ADJACENCY MATRIX (ERDOS-RENYI)

println("BUILDING NETWORK...")
adj_mat, rows = initialize_erdos_renyi(params);
println("NETWORK DONE")

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

# TIME AT WHICH LOCKDOWN IS ANNOUNCED
time_lockdown = 10

tau = []
num_S = []
num_I = []
num_R = []

# INITIALIZE TIME AND POPULATIONS
params.now_t = 0
populations = get_populations(agents, params)

println("SIM STARTS")

while populations[2] > 0.0 # WHILE THERE ARE INFECTED AGENTS

    # RECORD DATA
    push!(tau, params.now_t)
    push!(num_S, populations[1])
    push!(num_I, populations[2])
    push!(num_R, populations[3])

    # IF LOCKDOWN IS ANNOUNCED, SELECT COOPERATING AGENTS FROM POPULATION
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

println("NO MORE INFECTED!")

##=================####==============##
# EXPORT PARAMETERS
export_parameters(params, data_out_path)

##=================####==============##
# EXPORT POPULATIONS DATA

params_filename =
    "_I0_$(params.p_infected_t0)_AR_$(params.attack_rate)_CA_$(params.p_coop_agents)_$(params.repetition).csv"

file_name = "populations_N_$(params.num_agents)"*params_filename

open(joinpath(data_out_path, file_name), "w") do io
    println(io, "t,S,I,R")
    writedlm(io, hcat(tau, num_S, num_I, num_R), ',')
end

##=================####==============##
# EXPORT INFECTION TIME OF AGENTS

file_name = "infTime_N_$(params.num_agents)"*params_filename

open(joinpath(data_out_path, file_name), "w") do io
    println(io, "agent_ID,time_of_infection")
    for ag in agents
        println(io, "$(ag.id),$(ag.infection_t)")
    end
end

##=================####==============##

println("TAN TAN...")

##=================####==============##
