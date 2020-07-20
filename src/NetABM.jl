__precompile__()
module NetABM
    using LightGraphs, Random, DelimitedFiles, StatsBase, SparseArrays, Distributions
    export get_neigh, initialize_erdos_renyi, lectura_uw, get_network
    export initialize_demographics!, set_fixed_coop_agents!, set_coop_agents!
    export assign_contacts!, get_next_state!, update_state!
    export update_all_agents!, get_populations, export_parameters

    include("agents.jl")
    include("environmental.jl")
end # module