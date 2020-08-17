##=================####==============##

# MZH Agents Interactions

##=================####==============##

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
