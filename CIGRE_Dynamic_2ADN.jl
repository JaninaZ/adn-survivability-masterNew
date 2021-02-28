import Base: @__doc__

dir = @__DIR__
###### knoten ######

busses_static1, lines1, T1, elist1, Zs1, Yshs1 = CIGRE_static_ADN1()#1-12
busses_static2, lines2, T2, elist2, Zs2, Yshs2 = CIGRE_static_ADN2()#13-24
DCline = ACDCACLine()

function CIGRE_Dynamic_2ADN()


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
    P_ref(t) = t > 0.25 ? 0.5 * real(S_total) : real(S_total)  # try 0.1 or something small, vary the 1.1 here to see the change of P_g or P_int
    Q_ref(t) = t > 0.25 ? 0.5 * imag(S_total) : imag(S_total) 

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
tspan = (0, 10.0)
ode = rhs(cpg)
prob = ODEProblem(
    ode,
    op.vec, #op
    tspan)
# sol = solve(ode, op, tspan)
sol = solve(prob, Rodas5())
sol
end
