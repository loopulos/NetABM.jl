#  using NetABM, LightGraphs
using NetABM
using LightGraphs
using StatsBase
using Plots
using Statistics
using Bootstrap
using LaTeXStrings

num_agents = 1000
#  g = erdos_renyi(num_agents,0.4)
g = barabasi_albert(num_agents,9)
agents = [Agent(i) for i = 1:num_agents]
init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
set_coop_agents!(agents;p_cop=0.1)
map(x -> assign_contacts!(g,x), agents);
get_coop!(agents)


S = Array{Float64}(undef,0)
I = Array{Float64}(undef,0)
R = Array{Float64}(undef,0)

for i in 1:100
    next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.7, R=true)
    update_state!(agents)
    push!(S,([ag.state for ag in agents if ag.state=="S"] |> length)/num_agents)
    push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
    push!(R,([ag.state for ag in agents if ag.state=="R"] |> length)/num_agents)
end

plot(1:length(S),S,label="S")
plot!(1:length(S),I,label="I")
plot!(1:length(S),R,label="R")


#  What about levels of cooperation?
#


n_boot = num_agents
some_data = randn(100)
bs1 = bootstrap(std, some_data, BasicSampling(n_boot))
cil = 0.95
bci1 = confint(bs1, BasicConfInt(cil))
II = [zeros(100) for i in 1:100]

g = barabasi_albert(5000,10)
global p = plot()
begin
    Threads.@threads for j in [0.0,0.3,0.6,0.9]
        agents = [Agent(i) for i = 1:5000]
        init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
        set_coop_agents!(agents;p_cop=j)
        map(x -> assign_contacts!(g,x), agents);
        get_coop!(agents)
        I = Array{Float64}(undef,0)
        Threads.@threads for i in 1:100
            next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.7, R=true)
            update_state!(agents)
            push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/5000)
        end
        p = plot!(1:length(I),I,label=string(j))
    end
end
p

#  Changing the average degree of the network

num_agents = 1000
global Imean = Array{Float64}(undef,0)
global Ilow = Array{Float64}(undef,0)
global Ihigh = Array{Float64}(undef,0)
global Iall = Array{Array}(undef,0)
global n_boot = 10
global cil = 0.95
global j
        #  Threads.@threads
global p = plot()
begin
    for j in [2, 6, 10, 20]
        global Imean = Array{Float64}(undef,0)
        global Ilow = Array{Float64}(undef,0)
        global Ihigh = Array{Float64}(undef,0)
        global Iall = Array{Array}(undef,0)
        Threads.@threads for a in 1:100
            g = barabasi_albert(num_agents,j)
            agents = [Agent(i) for i = 1:num_agents]
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
            set_coop_agents!(agents;p_cop=0.9)
            map(x -> assign_contacts!(g,x), agents);
            get_coop!(agents)
            I = Array{Float64}(undef,0)
            for i in 1:100
                next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.5, R=true)
                update_state!(agents)
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
            end
            push!(Iall,I)
            #  p = plot!(1:length(I),I,label=string(j))
        end
        #  end
        for step in eachindex(Iall[1])
            bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
            bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            bci2 = confint(bs2, BasicConfInt(cil))
            push!(Imean,bci1[1][1])
            push!(Ilow,bci2[1][2])
            push!(Ihigh,bci2[1][3])
        end
        #  p=plot!(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.35,label=string(j))
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
end
p

#  step = 1
#  for step in eachindex(Iall[1])
#      bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
#      bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
#      bci1 = confint(bs1, BasicConfInt(cil))
#      bci2 = confint(bs2, BasicConfInt(cil))
#      push!(Imean,bci1[1][1])
#      push!(Ilow,bci2[1][2])
#      push!(Ihigh,bci2[1][3])
#  end
#
#  Imean
#  Ilow
#  Ihigh
#  plot(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.5)
#  plot(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.5)




#  Changing the cooperation according to neighboors

num_agents = 1000
global Imean = Array{Float64}(undef,0)
global Ilow = Array{Float64}(undef,0)
global Ihigh = Array{Float64}(undef,0)
global Iall = Array{Array}(undef,0)
global n_boot = 10
global cil = 0.95
global j
        #  Threads.@threads
global p = plot()
begin
    for j in [0.1, 0.4, 0.7]
        global Imean = Array{Float64}(undef,0)
        global Ilow = Array{Float64}(undef,0)
        global Ihigh = Array{Float64}(undef,0)
        global Iall = Array{Array}(undef,0)
        Threads.@threads for a in 1:100
            g = barabasi_albert(num_agents,10)
            agents = [Agent(i) for i = 1:num_agents]
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
            set_coop_agents!(agents;p_cop=0.4)
            set_adapt_agents!(agents;p_cop=j)
            map(x -> assign_contacts!(g,x), agents);
            get_coop!(agents)
            I = Array{Float64}(undef,0)
            for i in 1:100
                next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.8, R=true)
                update_state!(agents)
                update_coop!(agents,0.3)
                get_coop!(agents)
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
            end
            push!(Iall,I)
            #  p = plot!(1:length(I),I,label=string(j))
        end
        #  end
        for step in eachindex(Iall[1])
            bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
            bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            bci2 = confint(bs2, BasicConfInt(cil))
            push!(Imean,bci1[1][1])
            push!(Ilow,bci2[1][2])
            push!(Ihigh,bci2[1][3])
        end
        #  p=plot!(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.35,label=string(j))
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
end
p


#  Changing the cooperation according to previous infections

num_agents = 1000
global Imean = Array{Float64}(undef,0)
global Ilow = Array{Float64}(undef,0)
global Ihigh = Array{Float64}(undef,0)
global Iall = Array{Array}(undef,0)
global n_boot = 10
global cil = 0.95
global j
        #  Threads.@threads
global p = plot()
begin
    for j in [0.1, 0.4, 0.7, 1]
        global Imean = Array{Float64}(undef,0)
        global Ilow = Array{Float64}(undef,0)
        global Ihigh = Array{Float64}(undef,0)
        global Iall = Array{Array}(undef,0)
        Threads.@threads for a in 1:100
            g = barabasi_albert(num_agents,10)
            agents = [Agent(i) for i = 1:num_agents]
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
            set_coop_agents!(agents;p_cop=0.4)
            set_adapt_agents!(agents;p_cop=j)
            map(x -> assign_contacts!(g,x), agents);
            get_coop!(agents)
            I = Array{Float64}(undef,0)
            for i in 1:100
                next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.8, R=false)
                update_state!(agents)
                update_coop_infections!(agents,i;lrt=false)
                get_coop!(agents)
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
            end
            push!(Iall,I)
            #  p = plot!(1:length(I),I,label=string(j))
        end
        #  end
        for step in eachindex(Iall[1])
            bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
            bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            bci2 = confint(bs2, BasicConfInt(cil))
            push!(Imean,bci1[1][1])
            push!(Ilow,bci2[1][2])
            push!(Ihigh,bci2[1][3])
        end
        #  p=plot!(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.35,label=string(j))
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
end
p

p=plot!(xlims=(0,100),ylims=(0,1))
p=plot!(xlabel = L"t_i", ylabel = L"I/N")
savefig(p,"../figs/sis_adapting_infections.png")



#  Changing the cooperation according to previous infections including lrt

num_agents = 1000
global Imean = Array{Float64}(undef,0)
global Ilow = Array{Float64}(undef,0)
global Ihigh = Array{Float64}(undef,0)
global Iall = Array{Array}(undef,0)
global n_boot = 10
global cil = 0.95
global j
        #  Threads.@threads
global p = plot()
begin
    for j in [0.1, 0.4, 0.7, 1]
        global Imean = Array{Float64}(undef,0)
        global Ilow = Array{Float64}(undef,0)
        global Ihigh = Array{Float64}(undef,0)
        global Iall = Array{Array}(undef,0)
        Threads.@threads for a in 1:100
            g = barabasi_albert(num_agents,10)
            agents = [Agent(i) for i = 1:num_agents]
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
            set_coop_agents!(agents;p_cop=0.4)
            set_adapt_agents!(agents;p_cop=j)
            map(x -> assign_contacts!(g,x), agents);
            get_coop!(agents)
            I = Array{Float64}(undef,0)
            for i in 1:100
                next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.8, R=false)
                update_state!(agents)
                update_coop_infections!(agents,i;lrt=true)
                get_coop!(agents)
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
            end
            push!(Iall,I)
            #  p = plot!(1:length(I),I,label=string(j))
        end
        #  end
        for step in eachindex(Iall[1])
            bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
            bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            bci2 = confint(bs2, BasicConfInt(cil))
            push!(Imean,bci1[1][1])
            push!(Ilow,bci2[1][2])
            push!(Ihigh,bci2[1][3])
        end
        #  p=plot!(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.35,label=string(j))
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
end
p

p=plot!(xlims=(0,100),ylims=(0,1))
p=plot!(xlabel = L"t_i", ylabel = L"I/N")
savefig(p,"../figs/sir_adapting_infections_incl_lrt.png")



#  Changing the cooperation according to neighboors given distance

num_agents = 1000
global Imean = Array{Float64}(undef,0)
global Ilow = Array{Float64}(undef,0)
global Ihigh = Array{Float64}(undef,0)
global Iall = Array{Array}(undef,0)
global n_boot = 10
global cil = 0.95
global j
        #  Threads.@threads
global p = plot()
begin
    for j in [1,2,3]
        global Imean = Array{Float64}(undef,0)
        global Ilow = Array{Float64}(undef,0)
        global Ihigh = Array{Float64}(undef,0)
        global Iall = Array{Array}(undef,0)
        Threads.@threads for a in 1:100
            g = barabasi_albert(num_agents,10)
            agents = [Agent(i) for i = 1:num_agents]
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
            set_coop_agents!(agents;p_cop=0.4)
            set_adapt_agents!(agents;p_cop=0.4)
            map(x -> assign_contacts!(g,x), agents);
            get_coop!(agents)
            I = Array{Float64}(undef,0)
            for i in 1:100
                next_state!(agents;inf_prob=0.1, rec_prob=0.1, coop_red=0.8, R=true)
                update_state!(agents)
                update_coop_distance!(agents,g,j,0.3;lrt=true)
                get_coop!(agents)
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
            end
            push!(Iall,I)
            #  p = plot!(1:length(I),I,label=string(j))
        end
        #  end
        for step in eachindex(Iall[1])
            bs1 = bootstrap(mean,getindex.(Iall,step),BasicSampling(n_boot))
            bs2 = bootstrap(std,getindex.(Iall,step),BasicSampling(n_boot))
            bci1 = confint(bs1, BasicConfInt(cil))
            bci2 = confint(bs2, BasicConfInt(cil))
            push!(Imean,bci1[1][1])
            push!(Ilow,bci2[1][2])
            push!(Ihigh,bci2[1][3])
        end
        #  p=plot!(1:length(Imean), Imean, ribbon = (Imean .- Ilow, Ihigh .- Imean), fillalpha=0.35,label=string(j))
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
end
p

p=plot!(xlims=(0,100),ylims=(0,1))
p=plot!(xlabel = L"t_i", ylabel = L"I/N")
savefig(p,"../figs/sir_adapting_distance.png")






