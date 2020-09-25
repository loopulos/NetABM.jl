__precompile__()
module NetABM
using LightGraphs, Random, DelimitedFiles, StatsBase, SparseArrays, Distributions

export lectura_uw, lfr_network, Agent, Params, get_coop, set_adapt_agents!, update_coop!, set_coop_effect!,
    set_adapt_sign_agents!

export init_demographics!,
    get_coop!, next_state!, update_coop_infections!, update_coop_distance!
    #  , update_effect_given_distance_coop!, update_effect_given_distance_coop_free_dir!
    #  update_single_effect_distance_coop!, update_single_effect_distance_coop_free_dir!

export initialize_demographics!, set_fixed_coop_agents!, set_coop_agents!, SI_attitude!

export assign_contacts!,
    get_next_state!,
    update_state!,
    update_single_given_distance!,
    update_coop_given_distance!

export update_all_agents!,
    get_populations,
    export_parameters,
    update_single_effect_distance!,
    update_effect_given_distance!,
    update_single_effect_distance_coop!,
    update_single_effect_distance_coop_free_dir!,
    update_effect_given_distance_coop!,
    update_effect_given_distance_coop_free_dir!

"""
    Agent(id, neighs, p_cop, state)
Agent definition for the simulations with parameters:
`id          ::Int64 ` -> ID
`state       ::String ` -> Infection state (S, I, R)
`new_state   ::String ` -> Agent's state after meetings
`num_meets   ::Int64 ` -> Number of meetings in a time step
`recovery_t  ::Int64 ` -> Recovery time
`infection_t ::Int64 ` -> Time of infection
`contacts_t  ::Vector{Int64}` -> Contacts at time t
`age_group   ::Int64 ` -> Agent's age group (demographic)
`degree_t    ::Int64` -> Degree (Number of contacts) at time t
`ϵ           ::Float64` -> Bound of confidence
`p_cop       ::Float64 ` -> Cooperation probability
`attitude    ::String ` -> Risk aversion vs risk tolerant behavior: RA, RT
`coop_effect ::Float64 ` -> The mitigation efficacy of the complaince rules take by the agent
`at_home     ::Bool ` -> Flag to represent Agent is at home
`adapter     ::Bool ` -> Flag to represent Agent is willing to change behavior
`prior_bel   ::Float64 ` -> The mean value of the direction change of adaptive behavior
"""
mutable struct Agent
    id          ::Int64
    state       ::String
    new_state   ::String
    previous    ::Array{String}
    num_meets   ::Int64
    recovery_t  ::Int64
    infection_t ::Int64
    contacts_t  ::Vector{Int64}
    coopf       ::Vector{Int64}
    non_coopf   ::Vector{Int64}
    counter     ::Int64
    age_group   ::Int64
    degree_t    ::Int64
    ϵ           ::Float64
    p_cop       ::Float64
    attitude    ::String
    coop_effect ::Float64
    new_coop_effect ::Float64
    at_home     ::Bool
    adapter     ::Bool
    prior_bel   ::Float64
    # DEFAULT CONSTRUCTOR
    Agent(id) = new(
        id,
        "S",
        "S",
        Vector{String}(),
        1,
        5,
        0,
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        0,
        0,
        0,
        0.0,
        1.0,
        "rt",
        0.0,
        0.0,
        false,
        true,
        5.0
    )
end

##=================####==============##

"""
    Params(N, p_link, attack)
Global parameters for the simulation
`N::Int64 `  -> Number of agents in the population
`now_t::Int64 `  -> Current time (iteration)
`p_link::Float64` -> Link probability (Erdos-Renyi)
`p_infected_t0::Float64` -> Initial probability (fraction) of infected agents
`p_coop_agents::Float64` -> Initial probability (fraction) of lockdown agents
`attack_rate::Float64` -> Virus' attack rate
"""
mutable struct Params
    num_agents      ::Int64
    now_t           ::Int64
    μ_recovT_agents ::Int64
    p_link          ::Float64
    p_infected_t0   ::Float64
    μ_cop_agents    ::Float64
    σ_cop_agents    ::Float64
    attack_rate     ::Float64
    repetition      ::Int
    # DEFAULT CONSTRUCTOR
    Params() = new(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
end


#  include("agents.jl")
include("initialize.jl")
include("interactions.jl")
include("environmental.jl")
include("update_schemes.jl")
include("opinion_dynamics.jl")
include("data_out.jl")
end # module
