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
    #  if (ag.state == "S" && ag.new_state == "I")
    if ag.new_state == "I"
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

function update_single_effect_distance_coop_stochastic!(agents,g,d,threshold,steps;v)
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
        the_sign = rand(di) |> rand_sign
#         the_sign = rand() |> rand_sign
        if rand() <= 1/dist
            if agents[v].attitude == "ra"
                the_final = curr_effect - the_change * the_sign
            elseif agents[v].attitude == "rt"
                the_final = curr_effect + the_change * the_sign
            end
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

function update_effect_given_distance_coop_stochastic!(agents,g,d,threshold,step)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_coop_stochastic!(agents,g,d,threshold,step;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end

##=================####==============##

function update_single_effect_distance_coop_prior!(agents,g,d,threshold,steps,sd_change;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].coop_effect
    the_final = 0
    di = truncated(Normal(agents[v].prior_bel,sd_change),0,1)
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        mean_effect = mean([x.coop_effect for x in agents[current_nodes]])
        the_diff = abs(curr_effect - mean_effect)
        the_change = the_diff*steps[dist]
        the_sign = rand(di) |> rand_sign
        if rand() <= 1/dist
            if agents[v].attitude == "ra"
                the_final = curr_effect - (the_change * the_sign)
            elseif agents[v].attitude == "rt"
                the_final = curr_effect + (the_change * the_sign)
            end
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

function update_effect_given_distance_coop_prior!(agents,g,d,threshold,step,sd_change)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_coop_prior!(agents,g,d,threshold,step,sd_change;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
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

function update_single_effect_distance_alpha!(agents,g,d,threshold,steps,sd_change,sd_alpha;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].coop_effect
    the_final = 0
    di = truncated(Normal(agents[v].prior_bel,sd_change),0,1)
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        mean_effect = mean([x.coop_effect for x in agents[current_nodes]])
        the_diff = abs(curr_effect - mean_effect)
        if agents[v].attitude == "ra"
            side = -1
        elseif agents[v].attitude == "rt"
            side = 1
        end
        dialpha = truncated(Normal(steps[dist]*side,sd_alpha),-1,1)
        the_change = the_diff*rand(dialpha)
        the_sign = rand(di) |> rand_sign
        if rand() <= 1/dist
            the_final = curr_effect + (the_change * the_sign)
        end
    end
    if the_final > 1
        the_final = 1
    elseif the_final < 0
        the_final = 0
    end
    if the_final != agents[v].coop_effect
        agents[v].new_coop_effect = the_final
    else
        agents[v].new_coop_effect = agents[v].coop_effect
    end
end

##=================####==============##

function update_effect_given_distance_coop_alpha!(agents,g,d,threshold,step,sd_change,sd_alpha)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_alpha!(agents,g,d,threshold,step,sd_change,sd_alpha;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end

##=================####==============##

#  function risk_assessment!(agents, g, d; v)
#      neigh = neighborhood_dists(g,v,d)
#      nodes = last.(neigh)[2:end]
#      distances = last.(neigh)[2:end]
#      curr_effect = agents[v].coop_effect
#      the_final = 0
#      for dist in unique(distances)
#          current = findall(x->x==dist, distances)
#          current_nodes = nodes[current]
#          populations = get_populations(agents[current_nodes]);
#          I = get(populations,"I",0)
#          if I > 0
#              the_final = curr_effect + (1/dist)*(I/length(current))*0.01#(length(current)/nv(g))
#          else
#              the_final = curr_effect - (1/dist)*(length(current)/length(nodes))*0.01#(length(current)/nv(g))
#          end
#      end
#      if the_final > 1
#          the_final = 1
#      elseif the_final < 0
#          the_final = 0
#      end
#      if the_final != agents[v].coop_effect
#          agents[v].new_coop_effect = the_final
#      else
#          agents[v].new_coop_effect = agents[v].coop_effect
#      end
#  end

##=================####==============##

function update_risk_assessment!(agents,g,d)
    Threads.@threads for ag in agents#[the_adaps]
        risk_assessment(agents, g, d; v=ag.id)
    end
    Threads.@threads for ag in agents#[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end

##=================####==============##

function update_past!(agents,g,d,threshold,step,sd_change,sd_alpha)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_alpha!(agents,g,d,threshold,step,sd_change,sd_alpha;v=ag.id)
        end
    end
end

##=================####==============##

function risk_up!(agents, g, d; v)
    present = agents[v].coop_effect
    future = agents[v].new_coop_effect
    neigh = neighborhood_dists(g,v,d)
    nodes = last.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    the_prob = 0
    for dist in unique(distances)
        current = findall(x->x==dist, distances)
        current_nodes = nodes[current]
        populations = get_populations(agents[current_nodes]);
        I = get(populations,"I",0)
        if I > 0
            the_prob = (the_prob + (1/dist)*(I/length(current)))
        end
    end
    the_diff = abs(present - future)
    if rand() < the_prob
        updated = agents[v].coop_effect + the_diff
    if updated > 1
        updated = 1
    elseif updated < 0
        updated = 0
    end
    agents[v].new_coop_effect = updated
    end
end

##=================####==============##

function update_single_neighbors!(agents,g,d,threshold,steps,sd_change,sd_alpha;v)
    neigh = neighborhood_dists(g,v,d)
    nodes = first.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].coop_effect
    curr_bel = agents[v].prior_bel
    the_final = 0
    final_bel = 0
    di = truncated(Normal(agents[v].prior_bel,sd_change),0,1)
    for dist in unique(distances)
        current = findall(x->x==dist,distances)
        current_nodes = nodes[current]
        mean_effect = mean([x.coop_effect for x in agents[current_nodes]])
        the_diff = abs(curr_effect - mean_effect)
        populations = get_populations(agents[current_nodes]);
        I = get(populations,"I",0)
        if agents[v].attitude == "ra"
            side = -1
        elseif agents[v].attitude == "rt"
            side = 1
        end
        if I > 0
            final_bel = curr_bel + (steps[dist])*(I/length(current))
        else
            final_bel = curr_bel - (steps[dist])*(1/length(current))
        end
        dialpha = truncated(Normal(steps[dist]*side,sd_change),-1,1)
        the_change = the_diff*rand(dialpha)
        the_sign = rand(di) |> rand_sign
        if rand() <= 1/dist
            the_final = curr_effect + (the_change * the_sign)
        end
    end
    if the_final > 1
        the_final = 1
    elseif the_final < 0
        the_final = 0
    end
    if final_bel > 1
        final_bel = 1
    elseif final_bel < 0
        final_bel = 0
    end
    agents[v].new_coop_effect = the_final
    agents[v].prior_bel = final_bel
end

##=================####==============##

function update_all_neighbors!(agents,g,d,threshold,step,sd_change,sd_alpha)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            #  update_single_effect_distance_alpha!(agents,g,d,threshold,step,sd_change,sd_alpha;v=ag.id)
            update_single_neighbors!(agents,g,d,threshold,step,sd_change,sd_alpha;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end

##=================####==============##

function risk_assessment!(agents, g, d, risk, sd_risk; v)
    neigh = neighborhood_dists(g,v,d)
    nodes = last.(neigh)[2:end]
    distances = last.(neigh)[2:end]
    curr_effect = agents[v].new_coop_effect
    the_final = curr_effect
    for dist in unique(distances)
        current = findall(x->x==dist, distances)
        current_nodes = nodes[current]
        populations = get_populations(agents[current_nodes]);
        I = get(populations,"I",0)
        num = risk*(1/dist)*(I/length(current))
        di = truncated(Normal(num,sd_risk),-1,1)
        the_final = the_final + rand(di)
    end
    if the_final > 1
        the_final = 1
    elseif the_final < 0
        the_final = 0
    end
    if the_final != agents[v].new_coop_effect
        agents[v].new_coop_effect = the_final
    else
        agents[v].new_coop_effect = agents[v].new_coop_effect
    end
end

##=================####==============##

function update_effect_given_distance_coop_alpha_risk!(agents,g,d,threshold,step,sd_change,sd_alpha,risk, sd_risk)
    the_adaps = findall(x->x.adapter,agents)
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            update_single_effect_distance_alpha!(agents,g,d,threshold,step,sd_change,sd_alpha;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        if ag.adapter == true
            risk_assessment!(agents,g,d, risk, sd_risk;v=ag.id)
        end
    end
    Threads.@threads for ag in agents[the_adaps]
        ag.coop_effect = ag.new_coop_effect
    end
end


