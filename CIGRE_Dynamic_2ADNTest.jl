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
busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
DCline = ACDCACLine()

# function CIGRE_Dynamic_2ADN()


###### construct 2 ADNs and 1 DC bus line######
# this part does not differ much as before, but change the P_limit of the 2nd ADN
    busses_static = []
    append!(busses_static, busses_static1)
    append!(busses_static, busses_static2) # bus from 1 to 24
    lines = []
    append!(lines, lines1)
    append!(lines, lines2)
    append!(lines, DCline) # all the lines

    pg_static = PowerGrid(busses_static, lines) 
    # power_flow = pf_sol(pg_static, ones(48), nothing)
    power_flow = pf_sol(pg_static, initial_guess(pg_static), nothing)

    busses = copy(busses_static)

    DG_locs1 = collect(2:12) 
    DG_locs2 =  collect(14:24)
    DG_locs = []
    append!(DG_locs, DG_locs1) # located, 1 for slack, so 2 for MV1
    append!(DG_locs, DG_locs2)

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
    for i in DG_locs2 # this is a loop 

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

        P_limit=20.0, # pu  the value here is changed 

        Q_limit=1.0, # pu

        Vref=1.0, # pu

        Vdead=0.1, # pu

        S_pq=V -> S_bkgrnd * (quad * V^2 + 1 - quad), # quad can be chosen 1 or 0, voltage dependent or not

        Y_n=0.0,

    )


    end
   
    pg = PowerGrid(busses, lines)

    S_total = 2*complex(24.373141067954748, 6.115974820490637) * 1e6 / base_power  # the value to be checked
    P_ref(t) = t > 0.25 ?  0.7*real(S_total) : real(S_total)  # try 0.1 or something small, vary the 1.1 here to see the change of P_g or P_int
    Q_ref(t) = t > 0.25 ?  imag(S_total)-4 : imag(S_total) 

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
m = sol.prob.f.f.cpg.controller(sol(sol.t[1]), nothing, sol.t[1]) |> first
Δ = repeat([m],inner=length(sol.t))
# Δ =Any[]
for i in 1:length(sol.t)
    Δ[i] = sol.prob.f.f.cpg.controller(sol(sol.t[i]), nothing, sol.t[i]) |> first
end
Perr=[d[1][] for d in Δ]   # the index j is removed in the second blanket of array d
P_slack = zeros(length(sol.t),1)
for (k,P) in enumerate(Perr)
    P_slack[k,1] = P_ref(sol.t[k]).-Perr[k]
end
plot(sol.t, P_slack, title = "Power flow", 
    label = ["the real power"],xlabel=("time t"), ylabel=("P/u"))
plot!(sol.t, t->P_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Perr.png")
# plot!(sol.t, Perr)

Qerr=[d[2][] for d in Δ]
# plot(sol.t, Qerr)
Q_slack = zeros(length(sol.t),1)
for (k,Q) in enumerate(Qerr)
    Q_slack[k,1] = Q_ref(sol.t[k]).-Qerr[k]
end
# plot(sol.t,Q_slack)
# plot!(sol.t, t->Q_ref(t))
# savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_QerrSingle.png")
plot(sol.t, Q_slack, title = "Reactive Power flow", 
    label = ["the imag power"],xlabel=("time t"), ylabel=("Q/u"))
plot!(sol.t, t->Q_ref(t), label = ["the reference"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Qerr.png")
# Δ = sol.prob.f.f.cpg.controller(sol(sol.t), nothing, sol.t) |> first

# absolute value of voltage
# 1st ADN
U_real_slack=[u[1] for u in sol.u].^2
U_imag_slack= [u[2] for u in sol.u].^2
voltage_slack = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_slack[i] = sqrt.(U_real_slack[i] .+ U_imag_slack[i])
end

U_real_n2=[u[3] for u in sol.u].^2
U_imag_n2= [u[4] for u in sol.u].^2
voltage_n2 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n2[i] = sqrt.(U_real_n2[i] .+ U_imag_n2[i])
end

U_real_n3=[u[10] for u in sol.u].^2
U_imag_n3= [u[11] for u in sol.u].^2
voltage_n3 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n3[i] = sqrt.(U_real_n3[i] .+ U_imag_n3[i])
end

U_real_n4=[u[17] for u in sol.u].^2
U_imag_n4= [u[18] for u in sol.u].^2
voltage_n4 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n4[i] = sqrt.(U_real_n4[i] .+ U_imag_n4[i])
end

U_real_n5=[u[24] for u in sol.u].^2
U_imag_n5= [u[25] for u in sol.u].^2
voltage_n5 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n5[i] = sqrt.(U_real_n5[i] .+ U_imag_n5[i])
end

U_real_n6=[u[31] for u in sol.u].^2
U_imag_n6= [u[32] for u in sol.u].^2
voltage_n6 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n6[i] = sqrt.(U_real_n6[i] .+ U_imag_n6[i])
end

U_real_n7=[u[38] for u in sol.u].^2
U_imag_n7= [u[39] for u in sol.u].^2
voltage_n7 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n7[i] = sqrt.(U_real_n7[i] .+ U_imag_n7[i])
end

U_real_n8=[u[45] for u in sol.u].^2
U_imag_n8= [u[46] for u in sol.u].^2
voltage_n8 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n8[i] = sqrt.(U_real_n8[i] .+ U_imag_n8[i])
end

U_real_n9=[u[52] for u in sol.u].^2
U_imag_n9= [u[53] for u in sol.u].^2
voltage_n9 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n9[i] = sqrt.(U_real_n9[i] .+ U_imag_n9[i])
end

U_real_n10=[u[59] for u in sol.u].^2
U_imag_n10= [u[60] for u in sol.u].^2
voltage_n10 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n10[i] = sqrt.(U_real_n10[i] .+ U_imag_n10[i])
end

U_real_n11=[u[66] for u in sol.u].^2
U_imag_n11= [u[67] for u in sol.u].^2
voltage_n11 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n11[i] = sqrt.(U_real_n11[i] .+ U_imag_n11[i])
end

U_real_n12=[u[73] for u in sol.u].^2
U_imag_n12= [u[74] for u in sol.u].^2
voltage_n12 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n12[i] = sqrt.(U_real_n12[i] .+ U_imag_n12[i])
end

# 2nd ADN
U_real_PQ=[u[80] for u in sol.u].^2
U_imag_PQ= [u[81] for u in sol.u].^2
voltage_PQ = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_PQ[i] = sqrt.(U_real_PQ[i] .+ U_imag_PQ[i])
end

U_real_n14=[u[82] for u in sol.u].^2
U_imag_n14= [u[83] for u in sol.u].^2
voltage_n14 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n14[i] = sqrt.(U_real_n14[i] .+ U_imag_n14[i])
end

U_real_n15=[u[89] for u in sol.u].^2
U_imag_n15= [u[90] for u in sol.u].^2
voltage_n15 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n15[i] = sqrt.(U_real_n15[i] .+ U_imag_n15[i])
end

U_real_n16=[u[96] for u in sol.u].^2
U_imag_n16= [u[97] for u in sol.u].^2
voltage_n16 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n16[i] = sqrt.(U_real_n16[i] .+ U_imag_n16[i])
end

U_real_n17=[u[103] for u in sol.u].^2
U_imag_n17= [u[104] for u in sol.u].^2
voltage_n17 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n17[i] = sqrt.(U_real_n17[i] .+ U_imag_n17[i])
end

U_real_n18=[u[110] for u in sol.u].^2
U_imag_n18= [u[111] for u in sol.u].^2
voltage_n18 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n18[i] = sqrt.(U_real_n18[i] .+ U_imag_n18[i])
end

U_real_n19=[u[117] for u in sol.u].^2
U_imag_n19= [u[118] for u in sol.u].^2
voltage_n19 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n19[i] = sqrt.(U_real_n19[i] .+ U_imag_n19[i])
end

U_real_n20=[u[124] for u in sol.u].^2
U_imag_n20= [u[125] for u in sol.u].^2
voltage_n20 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n20[i] = sqrt.(U_real_n20[i] .+ U_imag_n20[i])
end

U_real_n21=[u[131] for u in sol.u].^2
U_imag_n21= [u[132] for u in sol.u].^2
voltage_n21 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n21[i] = sqrt.(U_real_n21[i] .+ U_imag_n21[i])
end

U_real_n22=[u[138] for u in sol.u].^2
U_imag_n22= [u[139] for u in sol.u].^2
voltage_n22 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n22[i] = sqrt.(U_real_n22[i] .+ U_imag_n22[i])
end

U_real_n23=[u[145] for u in sol.u].^2
U_imag_n23= [u[146] for u in sol.u].^2
voltage_n23 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n23[i] = sqrt.(U_real_n23[i] .+ U_imag_n23[i])
end

U_real_n24=[u[152] for u in sol.u].^2
U_imag_n24= [u[153] for u in sol.u].^2
voltage_n24 = zeros(length(U_real_slack))
for i in 1:length(U_real_slack)
    voltage_n24[i] = sqrt.(U_real_n24[i] .+ U_imag_n24[i])
end

# for (j,k) in ([3,10,17,24,31,38,45,52,59,66,73],[2:12])
#     U_real_n = zeros(length(U_real_slack))
#     U_imag_n = zeros(length(U_real_slack))
#     U_real_n[k]=[u[j] for u in sol.u].^2
#     U_imag_n[k]= [u[j+1] for u in sol.u].^2
#     voltage_n[k] = zeros(length(U_real_slack))
#     for i in 1:length(U_real_slack)
#         voltage_n[k][i] = sqrt.(U_real_slack[i] .+ U_imag_slack[i])
#     end
# end

result = DataFrame(   #data table
    t = sol.t,
    voltage_slack = voltage_slack,
    voltage_node2 = voltage_n2,
    voltage_node3 = voltage_n3,
    voltage_node4 = voltage_n4,
    voltage_node5 = voltage_n5,
    voltage_node6 = voltage_n6,
    voltage_node7 = voltage_n7,
    voltage_node8 = voltage_n8,
    voltage_node9 = voltage_n9,
    voltage_node10 = voltage_n10,
    voltage_node11 = voltage_n11,
    voltage_node12 = voltage_n12,
    voltage_PQ = voltage_PQ,
    voltage_node14 = voltage_n14,
    voltage_node15 = voltage_n15,
    voltage_node16 = voltage_n16,
    voltage_node17 = voltage_n17,
    voltage_node18 = voltage_n18,
    voltage_node19 = voltage_n19,
    voltage_node20 = voltage_n20,
    voltage_node21 = voltage_n21,
    voltage_node22 = voltage_n22,
    voltage_node23 = voltage_n23,
    voltage_node24 = voltage_n24,
)

plot(result.t, result.voltage_slack,
# ylims=(0.98,1.06),
title = "Voltage", 
label = ["voltage_slack"],xlabel=("time t"), ylabel=("V/u"))
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
plot!(result.t, result.voltage_PQ,label = ["voltage_nodePQ"])
plot!(result.t, result.voltage_node14,label = ["voltage_node14"])
plot!(result.t, result.voltage_node15,label = ["voltage_node15"])
plot!(result.t, result.voltage_node16,label = ["voltage_node16"])
plot!(result.t, result.voltage_node17,label = ["voltage_node17"])
plot!(result.t, result.voltage_node18,label = ["voltage_node18"])
plot!(result.t, result.voltage_node19,label = ["voltage_node19"])
plot!(result.t, result.voltage_node20,label = ["voltage_node20"])
plot!(result.t, result.voltage_node21,label = ["voltage_node21"])
plot!(result.t, result.voltage_node22,label = ["voltage_node22"])
plot!(result.t, result.voltage_node23,label = ["voltage_node23"])
plot!(result.t, result.voltage_node24,label = ["voltage_node24"])
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_VoltageAbsolute.png")

# the figure of voltage
plot(sol, vars=1)# slack bus 
plot!(sol, vars=3)#DG unit at node 2 and so on 
plot!(sol, vars=10)
plot!(sol, vars=17)
plot!(sol, vars=24)
plot!(sol, vars=31)
plot!(sol, vars=38)
plot!(sol, vars=45)
plot!(sol, vars=52)
plot!(sol, vars=59)
plot!(sol, vars=66)
plot!(sol, vars=73)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_voltageSingleReal.png")

plot(sol, vars=1:2)# slack bus 
plot!(sol, vars=3:4)#DG unit at node 2 and so on 
plot!(sol, vars=10:11)
plot!(sol, vars=17:18)
plot!(sol, vars=24:25)
plot!(sol, vars=31:32)
plot!(sol, vars=38:39)
plot!(sol, vars=45:46)
plot!(sol, vars=52:53)
plot!(sol, vars=59:60)
plot!(sol, vars=66:67)
plot!(sol, vars=73:74)
plot!(sol, vars=80:81)
plot!(sol, vars=82:83)
plot!(sol, vars=89:90)
plot!(sol, vars=96:97)
plot!(sol, vars=103:104)
plot!(sol, vars=110:111)
plot!(sol, vars=117:118)
plot!(sol, vars=124:125)
plot!(sol, vars=131:132)
plot!(sol, vars=138:139)
plot!(sol, vars=145:146)
plot!(sol, vars=152:153)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_voltageSingle.png")

# the figure of angle
plot(sol, vars=5)
plot!(sol, vars=12)
plot!(sol, vars=19)
plot!(sol, vars=26)
plot!(sol, vars=33)
plot!(sol, vars=40)
plot!(sol, vars=47)
plot!(sol, vars=54)
plot!(sol, vars=61)
plot!(sol, vars=68)
plot!(sol, vars=75)
plot!(sol, vars=84)
plot!(sol, vars=91)
plot!(sol, vars=98)
plot!(sol, vars=105)
plot!(sol, vars=112)
plot!(sol, vars=119)
plot!(sol, vars=126)
plot!(sol, vars=133)
plot!(sol, vars=140)
plot!(sol, vars=147)
plot!(sol, vars=154)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_angleSingle.png")

# the figure of P_int (= angle + 1)
plot(sol, vars=6, ylims=(-2,0))
plot!(sol, vars=13)
plot!(sol, vars=20)
plot!(sol, vars=27)
plot!(sol, vars=34)
plot!(sol, vars=41)
plot!(sol, vars=48)
plot!(sol, vars=55)
plot!(sol, vars=62)
plot!(sol, vars=69)
plot!(sol, vars=76)
plot!(sol, vars=85)
plot!(sol, vars=92)
plot!(sol, vars=99)
plot!(sol, vars=106)
plot!(sol, vars=113)
plot!(sol, vars=120)
plot!(sol, vars=127)
plot!(sol, vars=134)
plot!(sol, vars=141)
plot!(sol, vars=148)
plot!(sol, vars=155)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_P_intSingle.png")

# figure of Q_int (= P_int + 1)
plot(sol, vars=7, ylims=(-2,0))
plot!(sol, vars=14)
plot!(sol, vars=21)
plot!(sol, vars=28)
plot!(sol, vars=35)
plot!(sol, vars=42)
plot!(sol, vars=49)
plot!(sol, vars=56)
plot!(sol, vars=63)
plot!(sol, vars=70)
plot!(sol, vars=77)
plot!(sol, vars=86)
plot!(sol, vars=93)
plot!(sol, vars=100)
plot!(sol, vars=107)
plot!(sol, vars=114)
plot!(sol, vars=121)
plot!(sol, vars=128)
plot!(sol, vars=135)
plot!(sol, vars=142)
plot!(sol, vars=149)
plot!(sol, vars=156)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Q_intSingle.png")

# the figure of P_g (= Q_int + 1)
plot(sol, vars=8, ylims=(-2,0))
plot!(sol, vars=15)
plot!(sol, vars=22)
plot!(sol, vars=29)
plot!(sol, vars=36)
plot!(sol, vars=43)
plot!(sol, vars=50)
plot!(sol, vars=57)
plot!(sol, vars=64)
plot!(sol, vars=71)
plot!(sol, vars=78)
plot!(sol, vars=87)
plot!(sol, vars=94)
plot!(sol, vars=101)
plot!(sol, vars=108)
plot!(sol, vars=115)
plot!(sol, vars=122)
plot!(sol, vars=129)
plot!(sol, vars=136)
plot!(sol, vars=143)
plot!(sol, vars=150)
plot!(sol, vars=157)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_P_gSingle.png")

# the figure of Q_g (= P_g + 1)
plot(sol, vars=9)
plot!(sol, vars=16)
plot!(sol, vars=23)
plot!(sol, vars=30)
plot!(sol, vars=37)
plot!(sol, vars=44)
plot!(sol, vars=51)
plot!(sol, vars=58)
plot!(sol, vars=65)
plot!(sol, vars=72)
plot!(sol, vars=79)
plot!(sol, vars=88)
plot!(sol, vars=95)
plot!(sol, vars=102)
plot!(sol, vars=109)
plot!(sol, vars=116)
plot!(sol, vars=123)
plot!(sol, vars=130)
plot!(sol, vars=137)
plot!(sol, vars=144)
plot!(sol, vars=151)
plot!(sol, vars=158)
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Q_gSingle.png")