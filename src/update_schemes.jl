##=================####==============##

# Different update schemes

##=================####==============##

"""
    update_state!(agents)
Finds agent's next state and updates it, it just packages
`get_next_state!` and `update_state!`
"""
function update_single_state!(ag)
    push!(ag.previous, ag.state)
    if (ag.state == "S" && ag.new_state == "I")
        ag.counter = ag.counter + 1
    end
    ag.state = ag.new_state
end

function update_state!(agents)
    Threads.@threads for ag in agents
        update_single_state!(ag)
    end
end

##=================####==============##

function update_coop!(agents,threshold;lrt=false)
    Threads.@threads for ag in agents
        if ag.adapter
            if length(ag.coopf)/ag.degree_t > threshold
                ag.at_home = true
            else
                if lrt
                    ag.at_home = false
                end
            end
        end
    end
end

##=================####==============##

function update_coop_distance!(agents,g,d,threshold;lrt=false)
    Threads.@threads for ag in agents
        if ag.adapter
            status = agents[neighborhood(g,ag.id,d)|> unique] |> f -> map(x -> x.at_home,f)
            thecoop = sum(status)
            tot = length(status)
            if thecoop/tot >= threshold
                ag.at_home = true
            else
                if lrt
                    ag.at_home = false
                end
            end
        end
    end
end

##=================####==============##

function prev_infected(agent,threshold)
    historic = agent.previous[end-threshold:end]
    Infected = filter(x->x=="I",historic)
    Susceptible = filter(x->x=="S",historic)
    Recovered = filter(x->x=="R",historic)
    Infected, Susceptible, Recovered
end

##=================####==============##

function update_coop_infections!(agents,threshold;lrt=false)
    Threads.@threads for ag in agents
        if ag.adapter
            Infected, Susceptible, Recovered = prev_infected(ag, threshold)
            if length(Infected) >= 1
                ag.at_home = true
                #  ag.at_home = !(ag.at_home)
            else
                if lrt
                    ag.at_home = false
                    #  ag.at_home = !(ag.at_home)
                end
            end
        end
    end
end

##=================####==============##

function update_coop_distance_inf!(agents,g,d,threshold;lrt=false)
    Threads.@threads for ag in agents
        if ag.adapter
            status = agents[neighborhood(g,ag.id,d)|> unique] |> f -> map(x -> prev_infected(x,threshold)[1],f)
            theinfected = sum(map(x -> "I" in x, status))
            tot = length(status)
            if theinfected >= 1
                ag.at_home = true
            else
                if lrt
                    ag.at_home = false
                end
            end
        end
    end
end

##=================####==============##
##### Check infections and change behavior depending on distance

function update_single_given_distance!(agents,g,v,d,threshold,probs)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)
    distances = last.(neigh)
    changed = false
    for dist in 0:d
        current = findall(x->x==dist,distances)
        infected = findall(x->x.counter >= threshold, agents[current])
        for time in 1:length(infected)
            if rand() < probs[dist+1]
                agents[v].at_home = true
                changed = true
                break
            end
        end
        if changed
            break
        end
    end
end

##=================####==============##

function update_coop_given_distance!(agents,g,d,threshold,probs;lrt=false)
    Threads.@threads for ag in agents
        if ag.adapter
            update_single_given_distance!(agents,g,ag.id,d,threshold,probs)
        end
    end
end

##=================####==============##

function update_single_effect_distance!(agents,g,d,threshold,step;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)
    distances = last.(neigh)
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        infected = findall(x->x.counter >= threshold, agents[current_nodes])
        for time in 1:length(infected)
            if agents[v].attitude == "ra"
                agents[v].coop_effect = max(agents[v].coop_effect-step[dist+1],0)
            elseif agents[v].attitude == "rt"
                agents[v].coop_effect = min(agents[v].coop_effect+step[dist+1],1)
            end
        end
    end
end

##=================####==============##

function update_effect_given_distance!(agents,g,d,threshold,step)
    Threads.@threads for ag in agents
        if ag.adapter
            update_single_effect_distance!(agents,g,d,threshold,step;v=ag.id)
        end
    end
end

##=================####==============##

function update_single_effect_distance_coop!(agents,g,d,threshold,steps;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].coop_effect
    the_final = 0
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        mean_effect = mean([x.coop_effect for x in agents[current_nodes]])
        the_diff = abs(curr_effect - mean_effect)
        the_change = the_diff*steps[dist]
        if rand() <= 1/dist
            if agents[v].attitude == "ra"
                the_final = curr_effect - the_change
            elseif agents[v].attitude == "rt"
                the_final = curr_effect + the_change
            end
        end
    end
    if agents[v].attitude == "ra"
        agents[v].new_coop_effect = max(the_final,0)
    elseif agents[v].attitude == "rt"
        agents[v].new_coop_effect = min(the_final,1)
    end
end

##=================####==============##

function update_single_effect_distance_coop_free_dir!(agents,g,d,threshold,steps;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].coop_effect
    the_final = 0
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        mean_effect = mean([x.coop_effect for x in agents[current_nodes]])
        the_diff = (curr_effect - mean_effect)
        the_change = the_diff*steps[dist]
        if rand() <= 1/dist
            the_final = curr_effect - the_change
        end
    end
    if the_final > 1
        the_final = 1
    elseif the_final < 0
        the_final = 0
    end
    agents[v].new_coop_effect = the_final
end
##=================####==============##

function update_effect_given_distance_coop!(agents,g,d,threshold,step)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_coop!(agents,g,d,threshold,step;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end

##=================####==============##

function update_effect_given_distance_coop_free_dir!(agents,g,d,threshold,step)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_coop_free_dir!(agents,g,d,threshold,step;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
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
    update_all_agents!(agents, now_t)
Finds agent's next state and updates it, it just packages
`get_next_state!` and `update_state!`
"""
function update_all_agents!(agents, params)
    # FIND AGENTS' NEXT STATE FROM INTERACTIONS
    Threads.@threads for ag_id in 1:params.num_agents
        get_next_state!(agents, ag_id, params)
    end

    #  Threads.@threads for ag in agents
    #      get_next_state!(ag, now_t)
    #  end

    # UPDATE AGENTS' STATE
    Threads.@threads for ag in agents
        update_state!(ag)
    end
end

##=================####==============##
