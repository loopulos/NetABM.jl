##=================####==============##

## Functions to initalize Agents and Networks

##=================####==============##

"""
    lectura_uw(red)
Creates a LightGraph's SimpleGraph object from an adjacency list array. Returns the SimpleGraph and the ordered node sequence,
this way node names can be strings.
...
# Arguments
- `red::Array`: Adjacency list.
...
"""
function lectura_uw(red)
    Nodes = union(unique(red[:, 1]), unique(red[:, 2]))
    g = SimpleGraph()
    last_node = Int64(length(Nodes))
    add_vertices!(g, last_node)
    for n = 1:size(red)[1]
        add_edge!(g, red[n, 1], red[n, 2])
        add_edge!(g, red[n, 2], red[n, 1])
    end
    return g, Nodes
end

##=================####==============##

"""
Generates a LFR network
...
# Arguments
- `N::Int`: Number of nodes.
- `k::Float`: Average degree.
- `maxk::Float`: Maximum degree.
- `mu::Float`: Mixing parameter.
- `t1::Float`:: Minus exponent for the degree sequence.
- `t2::Float`:: Minus exponent for the community size distribution.
- `minc::Int`: Minimum for the community sizes.
- `maxc::Int`: Maximum for the community sizes.
...
"""
function lfr_network(
    the_root,
    work_path;
    N = 100,
    k = 20,
    maxk = 40,
    mu = 0.1,
    t1 = 2,
    t2 = 1,
    minc = 10,
    maxc = 40,
)
    # the_root = pwd()
    run(`$the_root/benchmark -N $N -k $k -maxk $maxk -mu $mu -t1 $t1 -t2 $t2 -minc $minc -maxc $maxc`)

    out_files_path = joinpath(the_root, "..", "..", work_path)

    # the_net = readdlm(the_root*"/network.dat")
    # true_com = readdlm(the_root*"/community.dat")

    the_net = readdlm(joinpath(out_files_path, "network.dat"))
    true_com = readdlm(joinpath(out_files_path, "community.dat"))

    true_com = Int64.(true_com)
    return the_net, true_com
end

##=================####==============##
#
"""
    init_demographics!(agents; kwargs...)
Initialize agents' demographic attributes and
fraction of infected agents at t = 0
`coop_dist` -> Cooperation Distribution
`meets_dist` -> Distribution of number of meetings per time step
`recovt_dist` -> Agent's recovery time distribution
"""
function init_demographics!(
    agents;
    states::Array{String} = ["S", "I"],
    initial::Array{Float64} = [0.5, 0.5],
    coop_dist = [0],
    meets_dist = [0],
    recovt_dist = [0],
)
    for ag in agents
        ag.state = sample(states, Weights(initial))
        push!(ag.previous, ag.state)
        ag.p_cop = rand(coop_dist)
        ag.num_meets = 1 + rand(meets_dist)
        ag.recovery_t = ceil(rand(recovt_dist))
        if ag.state == "I"
            ag.counter = ag.counter + 1
        end
    end
end

##=================####==============##

# Distribucion uniforme de p_cop
function mod_init_demographics!(
    agents,
    coop_dist,
    meets_dist,
    recovt_dist;
    states::Array{String} = ["S", "I"],
    initial::Array{Float64} = [0.5, 0.5],
    ϵ_cm::Float64 = 0.11, # close minded agents
    ϵ_om::Float64 = 0.22, # open minded agents
    ϵ_thr::Vector{Float64} = [1.0/3.0, 2.0/3.0]
)
    for ag in agents
        ag.state = sample(states, Weights(initial))
        push!(ag.previous, ag.state)

        #  ag.p_cop = rand(coop_dist)
        ag.p_cop = rand()
        if ag.p_cop <= ϵ_thr[1] || ag.p_cop >= ϵ_thr[2]
            ag.ϵ = ϵ_cm
        else
            ag.ϵ = ϵ_om
        end

        ag.num_meets = 1 + rand(meets_dist)
        ag.recovery_t = ceil(rand(recovt_dist))

        if ag.state == "I"
            ag.counter = ag.counter + 1
        end
    end
end

##=================####==============##

"""
    set_fixed_coop_agents!(agents, params)
Set cooperating agents with probability `p_coop_agents`
"""
function set_fixed_coop_agents!(agents; p_coop_agents = 0.0)
    for ag in agents
        if rand() <= p_coop_agents
            ag.at_home = true
        end
    end
end

##=================####==============##

"""
    set_coop_agents!(agents, params)
Agent stays at home with probability `Agent.p_cop` at each time step
"""
function set_coop_agents!(agents; p_cop = 0.5)
    Threads.@threads for ag in agents
        if rand() <= p_cop
            ag.at_home = true
            ag.attitude = "ra"
            ag.coop_effect = 1.0
            ag.new_coop_effect = 1.0
        else
            ag.at_home = false
        end
    end
end

##=================####==============##

"""
    set_coop_agents!(agents, params)
Agent stays at home with probability `Agent.p_cop` at each time step
"""
function set_athome_agents!(agents)
    Threads.@threads for ag in agents
        if rand() <= ag.p_cop
            ag.at_home = true
        else
            ag.at_home = false
        end
    end
end

##=================####==============##

"""
    set_adapt_agents!(agents, params)
Agent stays at home with probability `Agent.p_cop` at each time step
"""
function set_adapt_agents!(agents; p_cop = 0.5)
    Threads.@threads for ag in agents
        if rand() <= p_cop
            ag.adapter = true
        else
            ag.adapter = false
        end
    end
end

##=================####==============##

"""
    assign_contacts!(agent, all_agents, adj_mat, row)
Assign contacts to `agent` from adjacency matrix
"""
function assign_contacts!(g, agent)
    agent.contacts_t = neighbors(g, agent.id)
    agent.degree_t = degree(g, agent.id)
end

##=================####==============##

function get_coop!(agents)
    Threads.@threads for agent in agents
        agent.coopf = findall(x -> agents[x].at_home == true, agent.contacts_t)
        agent.non_coopf = findall(x -> agents[x].at_home == false, agent.contacts_t)
    end
end

##=================####==============##
