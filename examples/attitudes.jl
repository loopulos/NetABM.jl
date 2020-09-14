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

thep = [:LightGraphs, :StatsBase, :Plots, :GraphPlot, :Statistics, :LaTeXStrings, :Bootstrap, :Distributions, :NetABM]
useit(thep)

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

global adoptstep = (1,[0.1, 0.1])
global maxI = Array{Float64}(undef,0)
global sumI = []

numnei = Int.(ceil.(exp.(range(1,stop=4,length=11))))
numnei = numnei[5]
icol="#E83F6F"
scol="#05668D"
rcol="#90A959"

@time let
global plt = plot(
        xlabel = L"t",
        ylabel = L"N",
        guidefontsize = 16,
        tickfontsize = 10,
        dpi=300
    )

    #  for copera in 0.0:0.1:1
    for copera in [0.5]
        #  for adap in 0.0:0.1:1
        for adap in [0.5]
            global Smean = Array{Float64}(undef,0)
            global Slow = Array{Float64}(undef,0)
            global Shigh = Array{Float64}(undef,0)
            global Imean = Array{Float64}(undef,0)
            global Ilow = Array{Float64}(undef,0)
            global Ihigh = Array{Float64}(undef,0)
            global Rmean = Array{Float64}(undef,0)
            global Rlow = Array{Float64}(undef,0)
            global Rhigh = Array{Float64}(undef,0)
            global Iall = Array{Array}(undef,0)
            global Sall = Array{Array}(undef,0)
            global Rall = Array{Array}(undef,0)

            global RAmean = Array{Float64}(undef,0)
            global RAlow = Array{Float64}(undef,0)
            global RAhigh = Array{Float64}(undef,0)
            global LRAmean = Array{Float64}(undef,0)
            global LRAlow = Array{Float64}(undef,0)
            global LRAhigh = Array{Float64}(undef,0)
            global LRTmean = Array{Float64}(undef,0)
            global LRTlow = Array{Float64}(undef,0)
            global LRThigh = Array{Float64}(undef,0)
            global RTmean = Array{Float64}(undef,0)
            global RTlow = Array{Float64}(undef,0)
            global RThigh = Array{Float64}(undef,0)
            global RAall = Array{Array}(undef,0)
            global LRAall = Array{Array}(undef,0)
            global LRTall = Array{Array}(undef,0)
            global RTall = Array{Array}(undef,0)

            for a in 1:10
                g = barabasi_albert(num_agents,numnei)
                agents = [Agent(i) for i = 1:num_agents];
                init_demographics!(agents;states=["S","I"],initial=[1.0,0.0]);
                pat_0 = rand(1:nv(g))
                agents[pat_0].state = "I"
                set_coop_agents!(agents;p_cop=copera);
                set_adapt_agents!(agents;p_cop=adap);
                set_coop_effet!(agents)
                map(x -> assign_contacts!(g,x), agents);
                global num_S = []
                global num_I = []
                global num_R = []

                global num_ra = []
                global num_lra = []
                global num_lrt = []
                global num_rt = []
                for i in 1:100
                    push!(num_S, ([ag.state for ag in agents if ag.state=="S"] |> length)/num_agents)
                    push!(num_I, ([ag.state for ag in agents if ag.state=="I"] |> length)/num_agents)
                    push!(num_R, ([ag.state for ag in agents if ag.state=="R"] |> length)/num_agents)

                    push!(num_ra,  mean([ag.coop_effect for ag in agents if ag.attitude=="ra" && !ag.adapter]))
                    push!(num_lra,  mean([ag.coop_effect for ag in agents if (ag.attitude=="rt" && ag.adapter)]))
                    push!(num_lrt,  mean([ag.coop_effect for ag in agents if (ag.attitude=="ra" && ag.adapter)]))
                    push!(num_rt,  mean([ag.coop_effect for ag in agents if ag.attitude=="rt" && !ag.adapter]))

                    next_state!(agents;fun=SI_attitude!,inf_prob=0.1, rec_prob=0.03, R=true);
                    update_state!(agents);
                    update_effect_given_distance_coop!(agents,g,adoptstep[1],1,adoptstep[2]);
                end
                if maximum(num_I) > 0.04
                    plot!(plt, 1:100, num_I, alpha=0.15, color = icol, lw = 0.3, label="")
                    plot!(plt, 1:100, num_S, alpha=0.15, color = scol, lw = 0.3, label="")
                    plot!(plt, 1:100, num_R, alpha=0.15, color = rcol, lw = 0.3, label="")
                    push!(Iall,num_I)
                    push!(Sall,num_S)
                    push!(Rall,num_R)

                    push!(RAall,num_ra)
                    push!(LRAall,num_lra)
                    push!(LRTall,num_lrt)
                    push!(RTall,num_rt)
                end
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

                bs1 = bootstrap(mean,getindex.(Sall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(Sall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(Smean,bci1[1][1])
                push!(Slow,bci2[1][2])
                push!(Shigh,bci2[1][3])

                bs1 = bootstrap(mean,getindex.(Rall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(Rall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(Rmean,bci1[1][1])
                push!(Rlow,bci2[1][2])
                push!(Rhigh,bci2[1][3])

                bs1 = bootstrap(mean,getindex.(RAall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(RAall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(RAmean,bci1[1][1])
                push!(RAlow,bci2[1][2])
                push!(RAhigh,bci2[1][3])

                bs1 = bootstrap(mean,getindex.(LRAall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(LRAall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(LRAmean,bci1[1][1])
                push!(LRAlow,bci2[1][2])
                push!(LRAhigh,bci2[1][3])

                bs1 = bootstrap(mean,getindex.(LRTall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(LRTall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(LRTmean,bci1[1][1])
                push!(LRTlow,bci2[1][2])
                push!(LRThigh,bci2[1][3])

                bs1 = bootstrap(mean,getindex.(RTall,step),BasicSampling(n_boot))
                bs2 = bootstrap(std,getindex.(RTall,step),BasicSampling(n_boot))
                bci1 = confint(bs1, BasicConfInt(cil))
                bci2 = confint(bs2, BasicConfInt(cil))
                push!(RTmean,bci1[1][1])
                push!(RTlow,bci2[1][2])
                push!(RThigh,bci2[1][3])
            end
        end
        GC.gc()
    end
end

plot!(plt, 1:100, Imean, color = icol,   lab = "I", lw = 1.25)
plot!(plt, 1:100, Smean, color = scol,   lab = "S", lw = 1.25)
plot!(plt, 1:100, Rmean, color = rcol,   lab = "R", lw = 1.25)
plt

png(plt,"figs/cop05_adap05_k10.png")
savefig(plt,"figs/cop05_adap05_k10.pdf")



plt2 = plot(
        xlabel = L"t",
        ylabel = L"N",
        guidefontsize = 16,
        tickfontsize = 10,
        dpi=300
    )

plot!(plt2, 1:100, Imean, ribbon = (Ilow, Ihigh), fillalpha=0.35, color = "#E83F6F",   lab = "I", lw = 1.25)
plot!(plt2, 1:100, Smean, ribbon = (Slow, Shigh), fillalpha=0.35, color = "#05668D",   lab = "S", lw = 1.25)
plot!(plt2, 1:100, Rmean, ribbon = (Rlow, Rhigh), fillalpha=0.35, color = "#90A959",   lab = "R", lw = 1.25)

png(plt2,"figs/cop05_adap05_k10_95.png")
savefig(plt2,"figs/cop05_adap05_k10_95.pdf")
# plot!(plt, tau, num_R, color = :green, lab = "R", lw = 1.25)

plt3 = plot(
        xlabel = L"t",
        ylabel = L"N",
        guidefontsize = 16,
        tickfontsize = 10,
        dpi=300,
        ylim=[0,1]
    )

plot!(plt3, 1:100, RAmean, ribbon = (RAlow, RAhigh), fillalpha=0.35,    lab = "RA", lw = 1.25)
plot!(plt3, 1:100, LRAmean, ribbon = (LRAlow, LRAhigh), fillalpha=0.35,   lab = "LRA", lw = 1.25)
plot!(plt3, 1:100, LRTmean, ribbon = (LRTlow, LRThigh), fillalpha=0.35,   lab = "LRT", lw = 1.25)
plot!(plt3, 1:100, RTmean, ribbon = (RTlow, RThigh), fillalpha=0.35,   lab = "RT", lw = 1.25)



for i in 1:11
    j = numnei[i]
    z1 = z[:,:,i]
    #  plt = heatmap!(xs, ys, z1, title=L"\langle k\rangle =$i", aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
    plt = heatmap!(xs, ys, z1, title= h, aspect_ratio = 1, clim=(0,1), c = :algae,dpi=300)
    png(plt,"figs/heatmap_scaled_maximum_coop_neigh_$i.png")
    savefig(plt,"figs/heatmap_scaled_maximum_coop_neigh_$i.pdf")
end
