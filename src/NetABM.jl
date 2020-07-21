__precompile__()
module NetABM
    using LightGraphs, Random, DelimitedFiles, StatsBase, SparseArrays, Distributions
    export lectura_uw, lfr_network
    export init_demographics!
    export initialize_demographics!, set_fixed_coop_agents!, set_coop_agents!
    export assign_contacts!, get_next_state!, update_state!
    export update_all_agents!, get_populations, export_parameters

    """
        Agent(id, neighs, p_cop, state)
    Agent definition for the simulations with parameters:
    `id          ::Int64 ` -> ID
    `state       ::String ` -> Infection state (S, I, R)
    `new_state   ::String ` -> Agent's state after meetings
    `num_meets   ::Int64 ` -> Number of meetings in a time step
    `recovery_t  ::Int64 ` -> Recovery time
    `infection_t ::Int64 ` -> Time of infection
    `contacts_t  ::Vector{Agent}` -> Contacts at time t
    `age_group   ::Int64 ` -> Agent's age group (demographic)
    `degree_t    ::Int64` -> Degree (Number of contacts) at time t
    `p_cop       ::Float64 ` -> Cooperation probability
    `at_home     ::Bool ` -> Flag to represent Agent is at home
    """
    mutable struct Agent
        id          ::Int64
        state       ::String
        new_state   ::String
        num_meets   ::Int64
        recovery_t  ::Int64
        infection_t ::Int64
        contacts_t  ::Vector{Agent}
        age_group   ::Int64
        degree_t    ::Int64
        p_cop       ::Float64
        at_home     ::Bool
        # DEFAULT CONSTRUCTOR
        Agent(id) = new(id, "S", "S", 1, 5, 0, Vector{Agent}(), 1, 0, 1.0, false)
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
        num_agents    ::Int64
        now_t         ::Int64
        μ_recovT_agents::Int64
        p_link        ::Float64
        p_infected_t0 ::Float64
        μ_cop_agents ::Float64
        σ_cop_agents ::Float64
        attack_rate   ::Float64
        repetition    ::Int
        # DEFAULT CONSTRUCTOR
        Params() = new(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
    end


    include("agents.jl")
    include("environmental.jl")
end # module
