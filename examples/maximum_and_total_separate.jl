using Pkg
#  If NetABM is not available then it should be installed via
#  Uncomment if not installed
] add https://github.com/loopulos/NetABM.jl#develop
#  ] add https://github.com/ollin18/NetABM.jl

function useit(list::Array{Symbol})
    installed = [key for key in keys(Pkg.installed())]
    strpackages = @. string(list)
    uninstalled = setdiff(strpackages,installed)

    map(Pkg.add,uninstalled)
    for package ∈ list
        @eval using $package
    end
end

thep = [:LightGraphs, :StatsBase, :Plots, :GraphPlot, :Statistics, :LaTeXStrings, :Bootstrap, :NetABM]
useit(thep)


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

function get_network(;N = 100, k = 20, maxk = 40, mu = 0.1, t1 = 2, t2 = 1, minc = 10, maxc = 40, c=0.3)
        the_root = pwd()
        try
            run(`$the_root/src/binary_networks/benchmark -N $N -k $k -maxk $maxk -mu $mu -C $c`)
            #  run(`$the_root/src/binary_networks/benchmark -N $N -k $k -maxk $maxk -mu $mu`)
        catch
            N
        end
end
the_root = pwd()

get_network(N = 100, mu=0.1);

################################################################
################################################################
################################################################
################################################################
################################################################
# The same as before but instead of looking at neighbor infections
# we look use a mean cooperation effect by distance

global num_agents = 100
global n_boot = 1000
global cil = 0.95
global j
#  global muladol = [(0,[0.1]), (1,repeat([0.1],2)),(2,repeat([0.1],3)),(3,repeat([0.1],4))]
#  global muladol = [(i,repeat([0.6],i+1)) for i in 0:10]
#  global muladol = [(i,repeat([0.1],i+1)) for i in 0:10]
global muladol = [(i,repeat([.1],i)) for i in 1:6]
global muladol = [(i,LinRange(0.01,0.001,6)|>collect) for i in 1:6]

global adoptstep = (1,[0.1, 0.1])
global maxI = Array{Float64}(undef,0)
global sumI = []

@time let
    #  for adoptstep in muladol
    #  for numnei in 6:4:26
    #  for numnei in Int.(floor.(LinRange(4,80,6)))
    for numnei in Int.(ceil.(exp.(range(1,stop=4,length=11))))
    #  for mixing in [0.9, 0.1]#0.1:0.1:0.6
    for copera in 0.0:0.1:1
        for adap in 0.0:0.1:1
            global meantotI = Array{Float64}(undef,0)
            global Imean = Array{Float64}(undef,0)
            global Ilow = Array{Float64}(undef,0)
            global Ihigh = Array{Float64}(undef,0)
            global Itotall = Array{Float64}(undef,0)
            global Iall = Array{Array}(undef,0)
            global Sall = Array{Array}(undef,0)
            global Rall = Array{Array}(undef,0)
            for a in 1:15
                g = barabasi_albert(num_agents,numnei)
                        #  get_network(N = num_agents, mu=0.5, c=mixing);
                        #  the_net = readdlm(the_root*"/network.dat")
                        #  g, Nodes = lectura_uw(the_net)
                agents = [Agent(i) for i = 1:num_agents];
                init_demographics!(agents;states=["S","I"],initial=[1.0,0.0]);
                pat_0 = rand(1:nv(g))
                agents[pat_0].state = "I"
                #  init_demographics!(agents;states=["S","I"],initial=[0.9,0.1]);
                set_coop_agents!(agents;p_cop=copera);
                set_adapt_agents!(agents;p_cop=adap);
                set_coop_effect!(agents)
                map(x -> assign_contacts!(g,x), agents);
                num_S = []
                num_I = []
                num_R = []
                global oldS = ([ag.state for ag in agents if ag.state=="S"] |> length)/num_agents
                global oldI = ([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents
                global totI = ([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents
                for i in 1:100
                    push!(num_S, ([ag.state for ag in agents if ag.state=="S"] |> length)/num_agents)
                    push!(num_I, ([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
                    push!(num_R, ([ag.state for ag in agents if ag.state=="R"] |> length)/num_agents)

                    next_state!(agents;fun=SI_attitude!,inf_prob=0.1, rec_prob=0.03, R=true);
                    update_state!(agents);
                    update_effect_given_distance_coop!(agents,g,adoptstep[1],1,adoptstep[2]);

                    newI = oldS - num_S[end]
                    global totI = totI + newI
                    oldS = num_S[end]

                end
                push!(Iall,num_I)
                push!(Itotall,totI)
                GC.gc()
            end
            for step in eachindex(Iall[1])
                bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(Imean,bci1[1][1])
                push!(Ilow,bci2[1][2])
                push!(Ihigh,bci2[1][3])
            end
            bs1 = bootstrap(mean, Itotall,BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            push!(sumI,bci1[1][1])
            push!(maxI,maximum(Imean))
        end
        GC.gc()
    end
    end
end

theme(:default)
using DelimitedFiles
numnei = Int.(ceil.(exp.(range(1,stop=4,length=11))))

plt = heatmap(
        xlabel = "Compliant Agents",
        ylabel = "Adaptive Agents",
        title = "Maximum Infections",
        guidefontsize = 16,
        tickfontsize = 10,
        xlims = (0,11)
    )
xs = [string(i) for i = 0:0.1:1]
ys = [string(i) for i = 0:0.1:1]
z = reshape(maxI, (11, 11, 11))
    z1 = z[:,:,1]
mean(z1)
plt = heatmap!(xs, ys, z1, aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
png(plt,"figs/heatmap_scaled_maximum_coop.png")
savefig(plt,"figs/heatmap_scaled_maximum_coop.pdf")
writedlm("data/max_inf_coop_neigh_num.csv",maxI)
writedlm("data/max_inf_coop_neigh_exp.csv",maxI)

for i in 1:11
    j = numnei[i]
    z1 = z[:,:,i]
    #  plt = heatmap!(xs, ys, z1, title=L"\langle k\rangle =$i", aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
    plt = heatmap!(xs, ys, z1, aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
    png(plt,"figs/heatmap_scaled_maximum_coop_neigh_$i.png")
    savefig(plt,"figs/heatmap_scaled_maximum_coop_neigh_$i.pdf")
end


plt = heatmap(
        xlabel = "Compliant Agents",
        ylabel = "Adaptive Agents",
        title = "Maximum Infections",
        guidefontsize = 16,
        tickfontsize = 10,
        xlims = (0,11)
    )
xs = [string(i) for i = 0:0.1:1]
ys = [string(i) for i = 0:0.1:1]
z = reshape(sumI, (11, 11, 11))
    z1 = z[:,:,1]
mean(z1)
plt = heatmap!(xs, ys, z1, aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
png(plt,"figs/heatmap_scaled_maximum_coop.png")
savefig(plt,"figs/heatmap_scaled_maximum_coop.pdf")
writedlm("data/max_inf_coop_neigh_num.csv",maxI)
writedlm("data/max_inf_coop_neigh_exp.csv",maxI)
