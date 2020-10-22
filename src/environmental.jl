##=================####==============##

# OD Agents Interactions

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

##=================####==============##
function SI_attitude_old!(agent, agents;inf_prob=0.1, rec_prob=0.3, R=false)
    #  sus = findall(x -> x.state == "S",agents)
    if agent.state == "S"
        infec = [ag.coop_effect for ag in agents[agent.contacts_t] if ag.state == "I"]
        the_odds = @. inf_prob * (1-agent.coop_effect) * (1-infec)
        odds = the_odds |> f -> map(x -> sample([true,false], Weights([x,1-x])),f) |> sum
        dis = truncated(Normal(7,5),1,Inf)
        if odds > 0
            agent.new_state = "I"
            agent.days = Int(ceil(rand(dis)))
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

##=================####==============##

function SI_attitude!(agent, agents;inf_prob=0.3, rec_prob=0.03, R=false)
    #  sus = findall(x -> x.state == "S",agents)
    if agent.state == "S"
        infec = [ag.coop_effect for ag in agents[agent.contacts_t] if ag.state == "I"]
        the_odds = @. inf_prob * (1-agent.coop_effect) * (1-infec)
        odds = the_odds |> f -> map(x -> sample([true,false], Weights([x,1-x])),f) |> sum
        if odds > 0
            agent.new_state = "I"
        else
            agent.new_state = "S"
        end
    elseif agent.state == "I"
        if !R
            if agent.counter >= agent.recovery_t
                agent.new_state = "S"
            else
                agent.new_state = "I"
            end
        else
            if agent.counter >= agent.recovery_t
                agent.new_state = "R"
            else
                agent.new_state = "I"
            end
        end
    end
end

##=================####==============##

function next_state!(agents;fun=SI_coop!,kwargs...)
    for agent in agents
        fun(agent,agents;kwargs...)
    end
end

##=================####==============##

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

##=================####==============##

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

##=================####==============##
