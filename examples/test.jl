#  using NetABM, LightGraphs
using NetABM
using LightGraphs
using StatsBase
using Plots

g = erdos_renyi(1000,0.4)
g = barabasi_albert(1000,9)
agents = [Agent(i) for i = 1:1000]
init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
set_coop_agents!(agents;p_cop=0.1)
map(x -> assign_contacts!(g,x), agents);
get_coop!(agents)


S = Array{Int64}(undef,0)
I = Array{Int64}(undef,0)
R = Array{Int64}(undef,0)

while [ag.state for ag in agents if ag.state=="I"] |> length > 0
    next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.7, R=true)
    update_state!(agents)
    push!(S,[ag.state for ag in agents if ag.state=="S"] |> length)
    push!(I,[ag.state for ag in agents if ag.state=="I"] |> length)
    push!(R,[ag.state for ag in agents if ag.state=="R"] |> length)
end

plot(1:length(S),S,label="S")
plot!(1:length(S),I,label="I")
plot!(1:length(S),R,label="R")


