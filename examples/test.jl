#  using NetABM, LightGraphs
using LightGraphs
using StatsBase

g = erdos_renyi(100,0.7)
g = barabasi_albert(1000,9)
agents = [Agent(i) for i = 1:1000]
init_demographics!(agents;states=["S","I"],initial=[0.7,0.3])
set_coop_agents!(agents;p_cop=0.1)
map(x -> assign_contacts!(g,x), agents);
get_coop!(agents)

next_state!(agents;inf_prob=0.01, rec_prob=0.1, coop_red=0.5, R=true)
update_state!(agents)
[ag.state for ag in agents if ag.state=="S"] |> length
[ag.state for ag in agents if ag.state=="I"] |> length
[ag.state for ag in agents if ag.state=="R"] |> length


