### A Pluto.jl notebook ###
# v0.11.13

using Markdown
using InteractiveUtils

# ╔═╡ 28761bd0-f1fb-11ea-1f82-717a48c220a1
using Pkg

# ╔═╡ 6292b148-f1fb-11ea-2c32-91cd70cb3bd9
function useit(list::Array{Symbol})
    installed = [key for key in keys(Pkg.installed())]
    strpackages = @. string(list)
    uninstalled = setdiff(strpackages,installed)

    map(Pkg.add,uninstalled)
    for package ∈ list
        @eval using $package
    end
end

# ╔═╡ 036d9cfe-f1fc-11ea-2b6b-4fae4cf7dc02
Pkg.add(url="https://github.com/loopulos/NetABM.jl")

# ╔═╡ 62c3033e-f1fb-11ea-32e1-e94bbd4058f6
thep2 = [:LightGraphs, :StatsBase, :Plots, :GraphPlot, :Statistics, :LaTeXStrings, :Bootstrap]

# ╔═╡ 62c3d066-f1fb-11ea-3320-abbb1d14df89
useit(thep2)

# ╔═╡ 5eaa7916-f1fc-11ea-3354-fb5b33917bc7
global rtall = Array{Array}(undef,0);

# ╔═╡ 75ae33b4-f1fc-11ea-0ddb-9f912fed94e9
global lrtall = Array{Array}(undef,0);

# ╔═╡ 75b0a64e-f1fc-11ea-3245-ade5236eabf4
global raall = Array{Array}(undef,0);

# ╔═╡ 75c10a8e-f1fc-11ea-2949-5784fe7301d2
global lraall = Array{Array}(undef,0);

# ╔═╡ 75ebbe1e-f1fc-11ea-2d02-9fed550f2214
global n_boot = 100

# ╔═╡ 76001a58-f1fc-11ea-277c-e9e68b50fb9a
global cil = 0.95

# ╔═╡ 7600de2a-f1fc-11ea-07d8-c32270c2e0bd
global j

# ╔═╡ 76353418-f1fc-11ea-0fa4-73369a86b186
@time let
    #  for j in [0.0, 0.1, 0.4, 0.7, 1.0]
    for j in [0.5]
        Imean = Array{Float64}(undef,0)
        Ilow = Array{Float64}(undef,0)
        Ihigh = Array{Float64}(undef,0)
        Iall = Array{Array}(undef,0)
        for a in 1:100
            g = barabasi_albert(num_agents,3)
            agents = [Agent(i) for i = 1:num_agents];
            init_demographics!(agents;states=["S","I"],initial=[0.9,0.1]);
            set_coop_agents!(agents;p_cop=0.5);
            set_adapt_agents!(agents;p_cop=j);
            map(x -> assign_contacts!(g,x), agents);
            #  get_coop!(agents);
            I = Array{Float64}(undef,0);
            rt = Array{Float64}(undef,0);
            ra = Array{Float64}(undef,0);
            lrt = Array{Float64}(undef,0);
            lra = Array{Float64}(undef,0);
            for i in 1:100
                next_state!(agents;fun=SI_attitude!,inf_prob=0.1, rec_prob=0.03, R=false);
                update_state!(agents);
                update_effect_given_distance!(agents,g,adoptstep[1],1,adoptstep[2]);
                #  get_coop!(agents);
                push!(I,([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents);
                push!(rt,mean([ag.coop_effect for ag in agents if ag.attitude == "rt" && !ag.adapter]));
                push!(lrt,mean([ag.coop_effect for ag in agents if ag.attitude == "rt" && ag.adapter]));
                push!(ra,mean([ag.coop_effect for ag in agents if ag.attitude == "ra" && !ag.adapter]));
                push!(lra,mean([ag.coop_effect for ag in agents if ag.attitude == "ra" && ag.adapter]));
            end
            push!(Iall,I)
            push!(rtall,rt)
            push!(lrtall,lrt)
            push!(raall,ra)
            push!(lraall,lra)
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
            #  GC.gc()
        end
        p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label=string(j))
    end
    GC.gc()
end

# ╔═╡ 764b35e2-f1fc-11ea-0d88-630a537082e9
Imean = Array{Float64}(undef,0)

# ╔═╡ 764c297a-f1fc-11ea-110c-8181fb96b08c
Ilow = Array{Float64}(undef,0)

# ╔═╡ 76686946-f1fc-11ea-0fa2-b7edcca389b4
Ihigh = Array{Float64}(undef,0)

# ╔═╡ 7682b67a-f1fc-11ea-290e-3be1c6b72b2c
Iall = Array{Array}(undef,0)

# ╔═╡ 76840afc-f1fc-11ea-33d8-e7a199ac97f6
for step in eachindex(lraall[1])
    bs1 = bootstrap(mean,getindex.(lraall,step),BasicSampling(n_boot))
    bs2 = bootstrap(std,getindex.(lraall,step),BasicSampling(n_boot))
    bci1 = confint(bs1, BasicConfInt(cil))
    bci2 = confint(bs2, BasicConfInt(cil))
    push!(Imean,bci1[1][1])
    push!(Ilow,bci2[1][2])
    push!(Ihigh,bci2[1][3])
end

# ╔═╡ 76a007c0-f1fc-11ea-04cc-0d91a50209ff
pyplot()

# ╔═╡ 76cb3094-f1fc-11ea-347d-c74f8474d987
p

# ╔═╡ 8642ec88-f1fc-11ea-0b4c-1bea3b95596e
begin
		num_agents = 1000
		#  g = erdos_renyi(num_agents,0.4)
		g = barabasi_albert(num_agents,9)
		agents = [Agent(i) for i = 1:num_agents]
		init_demographics!(agents;states=["S","I"],initial=[0.9,0.1])
		set_coop_agents!(agents;p_cop=0.1)
		map(x -> assign_contacts!(g,x), agents);
		get_coop!(agents)
end

# ╔═╡ 76041c34-f1fc-11ea-16a6-3987f38e1d5b
global adoptstep = (0,[0.1])

# ╔═╡ 761bee22-f1fc-11ea-0a83-31849f4d3661
global p = plot(xlims=(0,100),ylims=(0,1),xlabel = L"t_i", ylabel = L"I/N")

# ╔═╡ 76b68ab8-f1fc-11ea-2aaa-711f0ef72d4f
p=plot!(1:length(Imean), Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35,label="lrt")

# ╔═╡ 75d8b670-f1fc-11ea-3505-996f9c062d6c
global num_agents = 1000

# ╔═╡ 761ab8f4-f1fc-11ea-303b-4d26188526b2
global p = plot(xlims=(0,100),ylims=(0,1),xlabel = L"t_i", ylabel = L"I/N")

# ╔═╡ 75c18024-f1fc-11ea-2fcd-a91bcc345fc0
global adoptstep = (1,[0.1,0.05])

# ╔═╡ Cell order:
# ╠═28761bd0-f1fb-11ea-1f82-717a48c220a1
# ╠═036d9cfe-f1fc-11ea-2b6b-4fae4cf7dc02
# ╠═6292b148-f1fb-11ea-2c32-91cd70cb3bd9
# ╠═62c3033e-f1fb-11ea-32e1-e94bbd4058f6
# ╠═62c3d066-f1fb-11ea-3320-abbb1d14df89
# ╠═8642ec88-f1fc-11ea-0b4c-1bea3b95596e
# ╠═5eaa7916-f1fc-11ea-3354-fb5b33917bc7
# ╠═75ae33b4-f1fc-11ea-0ddb-9f912fed94e9
# ╠═75b0a64e-f1fc-11ea-3245-ade5236eabf4
# ╠═75c10a8e-f1fc-11ea-2949-5784fe7301d2
# ╠═75c18024-f1fc-11ea-2fcd-a91bcc345fc0
# ╠═75d8b670-f1fc-11ea-3505-996f9c062d6c
# ╠═75ebbe1e-f1fc-11ea-2d02-9fed550f2214
# ╠═76001a58-f1fc-11ea-277c-e9e68b50fb9a
# ╠═7600de2a-f1fc-11ea-07d8-c32270c2e0bd
# ╠═76041c34-f1fc-11ea-16a6-3987f38e1d5b
# ╠═761ab8f4-f1fc-11ea-303b-4d26188526b2
# ╠═761bee22-f1fc-11ea-0a83-31849f4d3661
# ╠═76353418-f1fc-11ea-0fa4-73369a86b186
# ╠═764b35e2-f1fc-11ea-0d88-630a537082e9
# ╠═764c297a-f1fc-11ea-110c-8181fb96b08c
# ╠═76686946-f1fc-11ea-0fa2-b7edcca389b4
# ╠═7682b67a-f1fc-11ea-290e-3be1c6b72b2c
# ╠═76840afc-f1fc-11ea-33d8-e7a199ac97f6
# ╠═76a007c0-f1fc-11ea-04cc-0d91a50209ff
# ╠═76b68ab8-f1fc-11ea-2aaa-711f0ef72d4f
# ╠═76cb3094-f1fc-11ea-347d-c74f8474d987
