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
    Nodes = union(unique(red[:,1]),unique(red[:,2]))
    g = SimpleGraph()
    last_node = Int64(length(Nodes))
    add_vertices!(g,last_node)
    for n in 1:size(red)[1]
        add_edge!(g,red[n,1],red[n,2])
        add_edge!(g,red[n,2],red[n,1])
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
function lfr_network(the_root, work_path;N = 100, k = 20, maxk = 40, mu = 0.1, t1 = 2, t2 = 1, minc = 10, maxc = 40)
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
function init_demographics!(agents; states::Array{String} = ["S","I"], initial::Array{Float64} = [0.5,0.5], coop_dist=[0], meets_dist=[0], recovt_dist=[0])
    for ag in agents
        ag.state = sample(states, Weights(initial))
        push!(ag.previous, ag.state)
        ag.p_cop      = rand(coop_dist)
        ag.num_meets  = 1 + rand(meets_dist)
        ag.recovery_t = ceil(rand(recovt_dist))
    end
end

##=================####==============##

"""
    set_fixed_coop_agents!(agents, params)
Set cooperating agents with probability `p_coop_agents`
"""
function set_fixed_coop_agents!(agents; p_coop_agents=0.0)
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
        else
            ag.at_home = false
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

function get_coop!(agents)
    for agent in agents
        agent.coopf = [agents[i].id for i in agent.contacts_t if agents[i].at_home]
        agent.non_coopf = [agents[i].id for i in agent.contacts_t if !agents[i].at_home]
    end
end
##=================####==============##

"""
This function looks at the infected neighbors of each agent to compute wether or not
she will get infected. If both agents adopt a cooperative behavior then the probability
of infection gets reduced by both terms. As the SIS and SIR models are basically the same
excep for reinfections then only a flag for SIR is needed.
"""
function SI_coop!(agent, agents;inf_prob=0.1, rec_prob=0.3, coop_red=0.7, R=false)
    infected_coop = [ag.state for ag in agents[agent.coopf] if ag.state == "I"] |> length
    infected_noncoop = [ag.state for ag in agents[agent.non_coopf] if ag.state == "I"] |> length
    if agent.state == "S"
        if agent.at_home
            inf_prob_co = inf_prob * (1-coop_red)^2
            inf_prob_no = inf_prob * (1-coop_red)
        else
            inf_prob_co = inf_prob * (1-coop_red)
            inf_prob_no = inf_prob
        end
        oddsc = [sample([true,false],Weights([inf_prob_co,1-inf_prob_co])) for i in 1:infected_coop] |> sum
        oddsn = [sample([true,false],Weights([inf_prob_no,1-inf_prob_no])) for i in 1:infected_noncoop] |> sum
        odds = oddsc + oddsn
        if odds > 0
            agent.new_state = "I"
        else
            agent.new_state = "S"
        end
    elseif agent.state == "I"
        if !R
            if rand() <= rec_prob
                agent.new_state = "S"
            else
                agent.new_state = "I"
            end
        else
            if rand() <= rec_prob
                agent.new_state = "R"
            else
                agent.new_state = "I"
            end
        end
    end
end

function next_state!(agents;fun=SI_coop!,kwargs...)
    for agent in agents
        fun(agent,agents;kwargs...)
    end
end





function SI_next!(agents; inf_prob=0.1, rec_prob=0.3)
    for agent in agents
        if agent.state == "S"
            infected = [ag.state for ag in agents[agent.contacts_t] if ag.state == "I"] |> length
            odds = [sample([true,false],Weights([inf_prob,1-inf_prob])) for i in 1:infected] |> sum
            if odds > 0
                agent.new_state = "I"
            else
                agent.new_state = "S"
            end
        elseif agent.state == "I"
            if rand() <= rec_prob
                agent.new_state = "S"
            else
                agent.new_state = "I"
            end
        end
    end
end


function SIR_next!(agents; inf_prob=0.1, rec_prob=0.3)
    for agent in agents
        if agent.state == "S"
            infected = [ag.state for ag in agents[agent.contacts_t] if ag.state == "I"] |> length
            odds = [sample([true,false],Weights([inf_prob,1-inf_prob])) for i in 1:infected] |> sum
            if odds > 0
                agent.new_state = "I"
            else
                agent.new_state = "S"
            end
        elseif agent.state == "I"
            if rand() <= rec_prob
                agent.new_state = "R"
            else
                agent.new_state = "I"
            end
        end
    end
end


"""
    update_state!(agents)
Finds agent's next state and updates it, it just packages
`get_next_state!` and `update_state!`
"""
function update_state!(agents)
    Threads.@threads for ag in agents
        push!(ag.previous,ag.state)
        ag.state = ag.new_state
    end
end

function update_coop!(agents,threshoold)
    Threads.@threads for ag in agents
        if changer
            if length(ag.coopf)/ag.degree_t > threshoold
                ag.at_home = true
            else
                ag.at_home = false
            end
        end
    end
end





"""
    get_next_state!(agent, params)
Finds `agent.new_state` according to its contacts and meetings
"""
function get_next_state!(agent;now_t)

    # IF AGENT IS SUCEPTIBLE AND IS *NOT* AT HOME...
    if agent.state == "S" && agent.at_home == false

        if agent.num_meets <= agent.degree_t
            meets = sample(agent.contacts_t, agent.num_meets, replace=false)
        else
            meets = agent.contacts_t
        end

        # println(agent.id, "|meets:", [x.id for x in meets])

        for c_agent in meets # LOOP THROUGH AGENTS MET THAT TIMESTEP
            coll_state = agent.state * c_agent.state
            # SUCEPTIBLE <- INFECTED
            if coll_state == "SI" && c_agent.at_home == false
                # println(coll_state)
                if rand() < params.attack_rate
                    agent.new_state = "I"
                    agent.infection_t = now_t
                    break
                end
            end
        end
    elseif agent.state == "I" # INFECTED AGENT TO RECOVER

        time_from_infection = now_t - agent.infection_t
        if time_from_infection == agent.recovery_t
            agent.new_state = "R"
            # println(agent.id, "|", time_from_infection, "|RECOVERED!")
        end

    end
end

##=================####==============##

#  """
#      update_state!(agent)
#  Updates `agent.state` with the one computed in `get_next_state!`
#  """
#  function update_state!(agent)
#      agent.state = agent.new_state
#  end

##=================####==============##

"""
    update_all_agents!(agents, params)
Finds agent's next state and updates it, it just packages
`get_next_state!` and `update_state!`
"""
function update_all_agents!(agents, params)
    # FIND AGENTS' NEXT STATE FROM INTERACTIONS
    Threads.@threads for ag in agents
        get_next_state!(ag, params)
    end

    # UPDATE AGENTS' STATE
    Threads.@threads for ag in agents
        update_state!(ag)
    end
end

##=================####==============##

"""
    get_populations(agents, params)
Get number of agents in each state:
    - S: Suceptible
    - I: Infected
    - R: Recovered
    - Q: Quarantained (not implemented)
    - D: Decesead (not implemented)
"""
function get_populations(agents, params)
    # species = num_S, num_I, num_R, num_Q
    species = zeros(4)

    for ag in agents
        if ag.state == "S"
            species[1] += 1.0
        elseif ag.state == "I"
            species[2] += 1.0
        elseif ag.state == "R"
            species[3] += 1.0
        # else
        #     species[4] += 1.0
        end
    end
    return species ./ params.num_agents
end
##=================####==============##

"""
    export_parameters(params, out_path)
Export parameters to file
"""
function export_parameters(params, out_path)

    filename = "parameters_N_$(params.num_agents)_rep_$(params.repetition).txt"

    open(joinpath(out_path, filename), "w") do io
    println(io, "num_agents   :\t", params.num_agents   )
    println(io, "p_link       :\t", params.p_link       )
    println(io, "p_infected_t0:\t", params.p_infected_t0)
    println(io, "attack_rate  :\t", params.attack_rate  )
    println(io, "repetition   :\t", params.repetition   )
    end

end
##=================####==============##
