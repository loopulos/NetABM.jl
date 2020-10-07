##=================####==============##

# Data analysis, export SIR populations and parameters

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
function get_populations_old(agents, params)
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
    #  return species ./ params.num_agents
    return species
end
##=================####==============##

"""
    get_populations(agents)
It returns a dictionary with the number
of agents in each state. Might be:
    - S: Suceptible
    - I: Infected
    - R: Recovered
    - Q: Quarantained (not implemented)
    - D: Decesead (not implemented)
"""
function get_populations_old(agents)
    map(x -> x.state, agents) |> countmap
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
