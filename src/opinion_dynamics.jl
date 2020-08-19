##=================####==============##
# Opinion dynamics
##=================####==============##

"""
    HK_od_model!(agent;kwargs...)
Hegselmann-Krause model
Updates `agent` cooperation opinion as a function of its neighbors
ϵ -> Bounded confidence
"""
function HK_od_model!(agents, ag_id)

    # Agent to be updated
    ag = agents[ag_id]

    # Agent's confidence set
    ag_conf_set = [
        agents[neigh]
        for neigh in ag.contacts_t if abs(ag.p_cop - agents[neigh].p_cop) <= ag.ϵ
    ]

    if !isempty(ag_conf_set)
        # Mean opinion of agent's confidence set
        conf_set_op = mean([neigh.p_cop for neigh in ag_conf_set])

        # Fraction of infected agents in the confidence set
        ϕ_I =
            count(x -> x == "I", [neigh.state for neigh in ag_conf_set]) / length(ag_conf_set)

        # New agents opninion
        ag.p_cop = conf_set_op + 0.1 * ϕ_I * (1.0 - conf_set_op)
        #  ag.p_cop = conf_set_op
    end
end

##=================####==============##
