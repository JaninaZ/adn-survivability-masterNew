import Base: @__doc__
using PowerDynamics
dir = @__DIR__
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
include("$dir/CIGRE_Dynamic_2ADN.jl")
using Plots
# plot(sol,fmt=:png)
sol = CIGRE_Dynamic_2ADN()
# state = sol[end]
# Δ = sol.prob.f.f.cpg.controller(state, nothing, last(sol.t)) |> first
# ΔP = Δ[1]
# ΔQ = Δ[2]  # part of the result

# the figure of voltage
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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_voltage.png")

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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_angle.png")

# the figure of P_int (= angle + 1)
plot(sol, vars=6)
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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_P_int.png")

# figure of Q_int (= P_int + 1)
plot(sol, vars=7)
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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Q_int.png")

# the figure of P_g (= Q_int + 1)
plot(sol, vars=8)
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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_P_g.png")

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
savefig("$dir/simulationplots/CIGRE_Dynamic_2ADN_Q_g.png")



# result = DataFrame(   #data table
# 	R = R,
# 	T = T,
# 	P = valid_P,
# 	Q = valid_Q,
# 	nsc_surv = first.(sim.u),
# 	nsc_surv_t = [u[2] for u in sim.u], #sim is sol
# 	nsc_fin_v = [u[3] for u in sim.u],
# 	nsc_min_v = [u[4] for u in sim.u],
# 	nsc_success = categorical(last.(sim.u)),
#   
# )

# now = unix2datetime(time())
# CSV.write("$dir/mc_res_$(Dates.format(now, "e_d_u_yy-H_M_S"))_step_nsc_node_4_quad0_error.csv", result)

# resultPQ = CSV.read("$dir/results/mc_res_Fri_21_Aug_20-16_43_13_step_nsc_node_4_quad0_error.csv")

# plot(resultPQ.t, resultPQ.P) #plot the power over the time
