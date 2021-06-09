import Base: @__doc__
using PowerDynamics
dir = @__DIR__
###### knoten ######

include("$dir/control.jl")
# custom types
include("$dir/DGUnit.jl")
include("$dir/OLTC.jl")
# load actual data
#include("$dir/cigre_static.jl")
# failure model
#include("$dir/short_circuit.jl")  #there is problem when running through this julia file
# static solution
include("$dir/power_flow.jl")

include("$dir/CIGRE_static_2ADN.jl")
include("$dir/ACDCACLineU.jl")

busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
# busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
# DCline = ACDCACLine()

# function CIGRE_Dynamic_2ADN()


###### construct 2 ADNs and 1 DC bus line######
# this part does not differ much as before, but change the P_limit of the 2nd ADN
    busses_static = []
    append!(busses_static, busses_static1)
    # append!(busses_static, busses_static2) # bus from 1 to 24
    lines = []
    append!(lines, lines1)
    # append!(lines, lines2)
    # append!(lines, DCline) # all the lines

    pg_static = PowerGrid(busses_static, lines) 
    # power_flow = pf_sol(pg_static, ones(48), nothing)
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    busses = copy(busses_static)

    DG_locs1 = collect(2:12) 
    # DG_locs2 =  collect(14:24)
    DG_locs = []
    append!(DG_locs, DG_locs1) # located, 1 for slack, so 2 for MV1
    # append!(DG_locs, DG_locs2)

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

        S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

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

    S_total = complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power  # the value to be checked
    P_ref(t) = t > 0.25 ?  0.9*real(S_total) : real(S_total)  # try 0.1 or something small, vary the 1.1 here to see the change of P_g or P_int
    Q_ref(t) = t > 0.25 ? 0.9*imag(S_total) : imag(S_total) 

    # P_ref = 50 #ADN provide power to the higher grid
    # Q_ref = 0 # DC part has no reactive power
    cpg = ADN(pg, DGUnit, P_ref, Q_ref) 
    #cpg = ADN(pg, DGUnit,  P_ref(t), Q_ref(t))
    fp = initial_guess(cpg, power_flow[:, :u])
    #power_flow_modified = power_flow
    #power_flow_modified[:, :u][end] = power_flow[:,:u][12]
    #fp_mod = initial_guess(cpg, power_flow_modified[:, :u])
    
    op = find_steady_state(cpg, fp) 
    # initial_guess(pg)  return State(pg, sol.u)
    verbose ? check_operationpoint(cpg, op) : nothing   #how the verbose is defined?

    # return pg, cpg
#end


#fp = find_valid_initial_condition(pg, ones(163))
#fp = initial_guess(cpg)
tspan = (0, 2.0)
ode = rhs(cpg)
prob = ODEProblem(
    ode,
    op.vec, #op
    tspan)
# sol = solve(ode, op, tspan)
sol = solve(prob, Rodas5())

using Plots
# plot(sol,fmt=:png)
# sol = CIGRE_Dynamic_2ADN()
# state = sol[end]
# Δ = sol.prob.f.f.cpg.controller(state, nothing, last(sol.t)) |> first
# ΔP = Δ[1]
# ΔQ = Δ[2]  # part of the result

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
    P_slack[k,1] = P_ref(sol.t[k]).-Perr[k] # calculation here! deos the first real power have an error?
    println(P_slack)
end
plot(sol.t, P_slack, title = "Active Power flow", 
    label = ["the real power"],xlabel=("time step"), ylabel=("P/W"))
plot!(sol.t, t->P_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_PerrSingle0.9.png")
# plot!(sol.t, Perr)

###### reactive power flow at the slack ######
Qerr=[d[2][] for d in Δ]
# plot(sol.t, Qerr)
Q_slack = zeros(length(sol.t),1)
for (k,Q) in enumerate(Qerr)
    Q_slack[k,1] = Q_ref(sol.t[k]).-Qerr[k]
    println(Q_slack)
end

# plot(sol.t,Q_slack)
# plot!(sol.t, t->Q_ref(t))
# savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_QerrSingle.png")
plot(sol.t, Q_slack, title = "Reactive Power flow", 
    label = ["the imag power"],xlabel=("time step"), ylabel=("Q/Var"))
plot!(sol.t, t->Q_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_QerrSingle0.9.png")
# Δ = sol.prob.f.f.cpg.controller(sol(sol.t), nothing, sol.t) |> first

###### absolute value of voltage ######
# U_real_slack=[u[1] for u in sol.u].^2
# U_imag_slack= [u[2] for u in sol.u].^2
# voltage_slack = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_slack[i] = sqrt.(U_real_slack[i] .+ U_imag_slack[i])
# end

# U_real_n2=[u[3] for u in sol.u].^2
# U_imag_n2= [u[4] for u in sol.u].^2
# voltage_n2 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n2[i] = sqrt.(U_real_n2[i] .+ U_imag_n2[i])
# end

# U_real_n3=[u[10] for u in sol.u].^2
# U_imag_n3= [u[11] for u in sol.u].^2
# voltage_n3 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n3[i] = sqrt.(U_real_n3[i] .+ U_imag_n3[i])
# end

# U_real_n4=[u[17] for u in sol.u].^2
# U_imag_n4= [u[18] for u in sol.u].^2
# voltage_n4 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n4[i] = sqrt.(U_real_n4[i] .+ U_imag_n4[i])
# end

# U_real_n5=[u[24] for u in sol.u].^2
# U_imag_n5= [u[25] for u in sol.u].^2
# voltage_n5 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n5[i] = sqrt.(U_real_n5[i] .+ U_imag_n5[i])
# end

# U_real_n6=[u[31] for u in sol.u].^2
# U_imag_n6= [u[32] for u in sol.u].^2
# voltage_n6 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n6[i] = sqrt.(U_real_n6[i] .+ U_imag_n6[i])
# end

# U_real_n7=[u[38] for u in sol.u].^2
# U_imag_n7= [u[39] for u in sol.u].^2
# voltage_n7 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n7[i] = sqrt.(U_real_n7[i] .+ U_imag_n7[i])
# end

# U_real_n8=[u[45] for u in sol.u].^2
# U_imag_n8= [u[46] for u in sol.u].^2
# voltage_n8 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n8[i] = sqrt.(U_real_n8[i] .+ U_imag_n8[i])
# end

# U_real_n9=[u[52] for u in sol.u].^2
# U_imag_n9= [u[53] for u in sol.u].^2
# voltage_n9 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n9[i] = sqrt.(U_real_n9[i] .+ U_imag_n9[i])
# end

# U_real_n10=[u[59] for u in sol.u].^2
# U_imag_n10= [u[60] for u in sol.u].^2
# voltage_n10 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n10[i] = sqrt.(U_real_n10[i] .+ U_imag_n10[i])
# end

# U_real_n11=[u[66] for u in sol.u].^2
# U_imag_n11= [u[67] for u in sol.u].^2
# voltage_n11 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n11[i] = sqrt.(U_real_n11[i] .+ U_imag_n11[i])
# end

# U_real_n12=[u[73] for u in sol.u].^2
# U_imag_n12= [u[74] for u in sol.u].^2
# voltage_n12 = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_n12[i] = sqrt.(U_real_n12[i] .+ U_imag_n12[i])
# end

# U_real_slack=[u[1] for u in sol.u].^2
# U_imag_slack= [u[2] for u in sol.u].^2
# voltage_slack = zeros(length(U_real_slack))
# for i in 1:length(U_real_slack)
#     voltage_slack[i] = sqrt.(U_real_slack[i] .+ U_imag_slack[i])
# end
array_length = length([u[1] for u in sol.u])
node_length = length(busses)
U_real_n = zeros(node_length,array_length)
U_imag_n = zeros(node_length,array_length) 
voltage_n = zeros(node_length,array_length)
# U_real_n[1,:] = U_real_slack
# U_imag_n[1,:] = U_imag_slack
# voltage_n[1,:] = voltage_slack
for (j,k) in (enumerate([1,3,10,17,24,31,38,45,52,59,66,73]))
    
    U_real_n[j,:]=[u[k] for u in sol.u].^2
    U_imag_n[j,:]= [u[k+1] for u in sol.u].^2
    # println(U_real_n)
    # println(U_imag_n)
    # voltage_n[j+1,length(U_real_slack)] = zeros(length(U_real_slack))
    for i in 1:array_length
        voltage_n[j,i] = sqrt.(U_real_n[j,i] .+ U_imag_n[j,i])
    end
    # println(voltage_n)
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
    # voltage_node2 = voltage_n2,
    # voltage_node3 = voltage_n3,
    # voltage_node4 = voltage_n4,
    # voltage_node5 = voltage_n5,
    # voltage_node6 = voltage_n6,
    # voltage_node7 = voltage_n7,
    # voltage_node8 = voltage_n8,
    # voltage_node9 = voltage_n9,
    # voltage_node10 = voltage_n10,
    # voltage_node11 = voltage_n11,
    # voltage_node12 = voltage_n12,
)

plot(result.t, result.voltage_slack,ylims=(0.93,1.05),title = "Voltage", 
label = ["voltage_slack"],xlabel=("time step"), ylabel=("U/V"))
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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_VoltageAbsoluteSingle0.9.png")
