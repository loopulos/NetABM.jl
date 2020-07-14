module NetABM
    ##=================####==============##
    ## ABM NETWORK-EPIDEMICS
    ##=================####==============##

    ##=================####==============##

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

    ##=================####==============##

    """
        get_neigh(id, adj_mat, rows)
    Get contacts `id` from adjacency matrix
    `adj_mat::SparseArray ` -> Adjacency Matrix
    `rows::Array{Int,1}` -> `SparseArray` rows with non-zero values
    `agent_id::Int64` -> Agent `id` to find its contacts
    """
    function get_neigh(agent_id, adj_mat, rows)
        contacts = Int[]
        for j in nzrange(adj_mat, agent_id)
            push!(contacts, rows[j])
        end
        return contacts
    end

    ##=================####==============##

    """
        initialize_erdos_renyi(params)
    Generate adjacency matrix for a random (undirected) Erdos-Renyi graph
    with `params.num_agents` nodes and `params.p_link` edge probability
    Returns:
    `adj_mat::SparseArray ` -> Adjacency Matrix
    `rows::Array{Int,1}` -> `SparseArray` rows with non-zero values
    """
    function initialize_erdos_renyi(params)
        adj_mat = spzeros(Int8, params.num_agents, params.num_agents)

        for i in 1:params.num_agents, j in i+1:params.num_agents
            if rand() < params.p_link
                adj_mat[i,j] = 1
                adj_mat[j,i] = 1
            end
        end

        rows = rowvals(adj_mat)
        return adj_mat, rows
    end

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
    function get_network(the_root, work_path;N = 100, k = 20, maxk = 40, mu = 0.1, t1 = 2, t2 = 1, minc = 10, maxc = 40)
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

    """
        initialize_demographics!(agents, params)
    Initialize agents' demographic attributes and
    fraction of infected agents at t = 0
    `coop_dist` -> Distribution of cooperation
    `meets_dist` -> Distribution of number of meetings per time step
    `recovt_dist` -> Agent's recovery time distribution
    """
    function initialize_demographics!(agents, params, coop_dist, meets_dist, recovt_dist)

        for ag in agents
            if rand() <= params.p_infected_t0
                ag.state = "I"
            end
            ag.p_cop      = rand(coop_dist)
            ag.num_meets  = 1 + rand(meets_dist)
            ag.recovery_t = ceil(rand(recovt_dist))
        end
    end

    ##=================####==============##

    """
        set_fixed_coop_agents!(agents, params)
    Set cooperating agents with probability `params.p_coop_agents`
    """
    function set_fixed_coop_agents!(agents, params)

        for ag in agents
            if rand() <= params.p_coop_agents
                ag.at_home = true
            end
        end

    end

    ##=================####==============##

    """
        set_coop_agents!(agents, params)
    Agent stays at home with probability `Agent.p_cop` at each time step
    """
    function set_coop_agents!(agents, params)

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
        assign_contacts!(agent, all_agents, adj_mat, row)
    Assign contacts to `agent` from adjacency matrix
    """
    function assign_contacts!(agent, all_agents, adj_mat, rows)

        agent.contacts_t = Vector{Agent}()

        neigh = get_neigh(agent.id, adj_mat, rows)
        map(x->push!(agent.contacts_t, all_agents[x]), neigh)
        agent.degree_t = length(agent.contacts_t)
        # println(neigh)
    end

    ##=================####==============##

    """
        get_next_state!(agent, params)
    Finds `agent.new_state` according to its contacts and meetings
    """
    function get_next_state!(agent, params)

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
                        agent.infection_t = params.now_t
                        break
                    end
                end
            end
        elseif agent.state == "I" # INFECTED AGENT TO RECOVER

            time_from_infection = params.now_t - agent.infection_t
            if time_from_infection == agent.recovery_t
                agent.new_state = "R"
                # println(agent.id, "|", time_from_infection, "|RECOVERED!")
            end

        end
    end

    ##=================####==============##

    """
        update_state!(agent)
    Updates `agent.state` with the one computed in `get_next_state!`
    """
    function update_state!(agent)
        agent.state = agent.new_state
    end

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
end # module
