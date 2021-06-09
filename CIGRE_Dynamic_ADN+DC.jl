import Base: @__doc__
using PowerDynamics
dir = @__DIR__
###### knoten ######

include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
include("$dir/power_flow.jl")

include("$dir/CIGRE_static_2ADN.jl")
include("$dir/ACDCACLineU.jl")

busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
# busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
busses_static2 = [PQAlgebraic(P = -0.5,Q = -0.5)]
# busses_static2 = [PQAlgebraic(P = t->-0.5,Q = t->-0.5)]
# busses_static2 = [VSIMinimal(;τ_P=0.632,τ_Q=0.632,K_P=1.396,K_Q=1.198,V_r=0.5,P=-0.5,Q=-0.5)]
DCline = ACDCACLine()
    # C_dc = ω * c * transmission_length / 4
    # C_conv = ω * 20E-6
    # Z_dcconv =  r * transmission_length + 1im * (ω * l * transmission_length - 1 / (2 * C_dc + 2 * C_conv))
    # C1 = 0.1511749E-6 
    # ldata =  2.82
    # Ysh = 1im .* ω .* C1 .* ldata
    # linetest = [PiModelLine(;from=1, to = 13 ,y=Z_dcconv ,y_shunt_km = Ysh / 2.0 / base_admittance,
    # y_shunt_mk = Ysh / 2.0 / base_admittance, )]  #instead of the DC line

    # DCline = [PiModelLine(;from=1, to = 13 ,y=(5.45385-1im*18.54203)*100 ,y_shunt_km = 8.63937E-5,
    # y_shunt_mk =8.63937E-5, )]  #DC line of 100km long move into the ACDCACLine.jl



    # R1 = 0.501 # Ω/km
    # X1 = 0.716 # Ω/km
    # C1 = 0.1511749E-6 # F/km

    # Z = complex(R1, X1) .* 3 
    # Ysh = 1im .* ω .* C1 .* 2
    # linetest = [PiModelLine(;from=1, to = 13 ,
    # y = inv(Z) / base_admittance,
    # y_shunt_km = Ysh / 2.0 / base_admittance,
    # y_shunt_mk = Ysh / 2.0 / base_admittance, )]  #instead of the DC line
    
# function CIGRE_Dynamic_2ADN()


###### construct 2 ADNs and 1 DC bus line######
# this part does not differ much as before, but change the P_limit of the 2nd ADN
    busses_static = []
    append!(busses_static, busses_static1)
    append!(busses_static, busses_static2) # bus from 1 to 24
    lines = []
    append!(lines, lines1)
    # append!(lines, lines2)
    append!(lines, DCline) # all the lines
    # append!(lines, linetest)

    pg_static = PowerGrid(busses_static, lines) 
    # power_flow = pf_sol(pg_static, ones(48), nothing)
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    busses = copy(busses_static)

    DG_locs1 = collect(2:12) 
    # DG_locs2 =  collect(14:24)
    DG_locs = []
    append!(DG_locs, DG_locs1) # located, 1 for slack, so 2 for MV1
    # append!(DG_locs, DG_locs2)
    quad = 0
    # two loops for two ADNs(set the DG unit)
    for i in DG_locs1 # this is a loop 

        S_bkgrnd = zero(im)

    try

        S_bkgrnd = busses_static[i].S

    catch

        S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)

    end

        busses[i] = DGUnit(;

        K_pll=1632.993, # Hz/pu

        T_int=2.0,

        K_PT1=1.0,

        T_PT1=1e-8,

        K_FRT=2.0,

        I_max=1.0, # pu

        P_limit=1.0, # pu   the value here is changed

        Q_limit=1.0, # pu

        Vref=1.0, # pu

        Vdead=0.1, # pu

        S_pq=V -> S_bkgrnd, #* (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

        Y_n=0.0,

    )


    end
    # for i in DG_locs2 # this is a loop 

    #     S_bkgrnd = zero(im)

    # try

    #     S_bkgrnd = busses_static[i].S

    # catch

    #     S_bkgrnd = complex(busses_static[i].P, busses_static[i].Q)

    # end

    #     busses[i] = DGUnit(;

    #     K_pll=1632.993, # Hz/pu

    #     T_int=2.0,

    #     K_PT1=1.0,

    #     T_PT1=1e-8,

    #     K_FRT=2.0,

    #     I_max=1.0, # pu

    #     P_limit=20.0, # pu  the value here is changed 

    #     Q_limit=1.0, # pu

    #     Vref=1.0, # pu

    #     Vdead=0.1, # pu

    #     S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

    #     Y_n=0.0,

    # )


    # end
   
    pg = PowerGrid(busses, lines)

    # S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power  # the value to be checked
    S_total = complex(24.309026442693472, 6.116128028641012) * 1e6 / base_power  # the largest value the system can reach after added the 13th node
    # S_total = complex(22.37314104875412, 6.115974820490637) * 1e6 / base_power 
    P_ref(t) = t > 1.5 ?  0.75*real(S_total) : real(S_total) # try 0.1 or something small, vary the 1.1 here to see the change of P_g or P_int
    Q_ref(t) = t > 1.5 ?  0.75*imag(S_total) : imag(S_total) 


    cpg = ADN(pg, DGUnit, P_ref, Q_ref) 
    fp = initial_guess(cpg, power_flow[:, :u])
    
    op = find_steady_state(cpg, fp) 
    verbose ? check_operationpoint(cpg, op) : nothing   #how the verbose is defined?

    # return pg, cpg
#end


tspan = (0, 5.0)
ode = rhs(cpg)
prob = ODEProblem(
    ode,
    op.vec, #op
    tspan)
# sol = solve(ode, op, tspan)
sol = solve(prob, Rodas5())

using Plots

###### active power flow at the slack ######
m = sol.prob.f.f.cpg.controller(sol(sol.t[1]), nothing, sol.t[1]) |> first
Δ = repeat([m],inner=length(sol.t))
# Δ =Any[]
for i in 1:length(sol.t)
    Δ[i] = sol.prob.f.f.cpg.controller(sol(sol.t[i]), nothing, sol.t[i]) |> first
end
Perr=[d[1][] for d in Δ]   # the index j is removed in the second blanket of array d
P_slack = zeros(length(sol.t),1)
# P_slack[1,1] = P_ref(sol.t[1])
for (k,P) in enumerate(Perr)
    P_slack[k,1] = P_ref(sol.t[k]).-Perr[k] 
println(P_slack)
end
plot(sol.t, P_slack, title = "Active Power flow", 
    label = ["the real power"],xlabel=("time step"), ylabel=("P/W"))
plot!(sol.t, t->P_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_PerrSingle+DC0.75.png")
# plot!(sol.t, Perr)

###### reactive power flow at the slack ######
Qerr=[d[2][] for d in Δ]
# plot(sol.t, Qerr)
Q_slack = zeros(length(sol.t),1)
for (k,Q) in enumerate(Qerr)
    Q_slack[k,1] = Q_ref(sol.t[k]).-Qerr[k]
    println(Q_slack)
end

# savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_QerrSingle.png")
plot(sol.t, Q_slack, title = "Reactive Power flow", 
    label = ["the imag power"],xlabel=("time step"), ylabel=("Q/Var"))
plot!(sol.t, t->Q_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_QerrSingle+DC0.75.png")
# Δ = sol.prob.f.f.cpg.controller(sol(sol.t), nothing, sol.t) |> first

###### absolute value of voltage ######
array_length = length([u[1] for u in sol.u])
node_length = length(busses)
U_real_n = zeros(node_length,array_length)
U_imag_n = zeros(node_length,array_length) 
voltage_n = zeros(node_length,array_length)
for (j,k) in (enumerate([1,3,10,17,24,31,38,45,52,59,66,73,80]))
    
    U_real_n[j,:]=[u[k] for u in sol.u].^2
    U_imag_n[j,:]= [u[k+1] for u in sol.u].^2
    for i in 1:array_length
        voltage_n[j,i] = sqrt.(U_real_n[j,i] .+ U_imag_n[j,i])
    end
end

result = DataFrame(   #data table
    t = sol.t,
    voltage_slack = voltage_n[1,:],
    voltage_node2 = voltage_n[2,:],
    voltage_node3 = voltage_n[3,:],
    voltage_node4 = voltage_n[4,:],
    voltage_node5 = voltage_n[5,:],
    voltage_node6 = voltage_n[6,:],
    voltage_node7 = voltage_n[7,:],
    voltage_node8 = voltage_n[8,:],
    voltage_node9 = voltage_n[9,:],
    voltage_node10 = voltage_n[10,:],
    voltage_node11 = voltage_n[11,:],
    voltage_node12 = voltage_n[12,:],
    voltage_node13 = voltage_n[13,:],
)

plot(result.t, result.voltage_slack,ylims=(0.93,1.05),title = "Voltage", 
label = ["voltage_slack"],xlabel=("time t"), ylabel=("U/V"))
plot!(result.t, result.voltage_node2,label = ["voltage_node2"])
plot!(result.t, result.voltage_node3,label = ["voltage_node3"])
plot!(result.t, result.voltage_node4,label = ["voltage_node4"])
plot!(result.t, result.voltage_node5,label = ["voltage_node5"])
plot!(result.t, result.voltage_node6,label = ["voltage_node6"])
plot!(result.t, result.voltage_node7,label = ["voltage_node7"])
plot!(result.t, result.voltage_node8,label = ["voltage_node8"])
plot!(result.t, result.voltage_node9,label = ["voltage_node9"])
plot!(result.t, result.voltage_node10,label = ["voltage_node10"])
plot!(result.t, result.voltage_node11,label = ["voltage_node11"])
plot!(result.t, result.voltage_node12,label = ["voltage_node12"])
plot!(result.t, result.voltage_node13,label = ["voltage_node13"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_VoltageAbsoluteSingle+DC0.75.png")