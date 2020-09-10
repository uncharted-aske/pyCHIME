using DifferentialEquations
using Plots

function sir_ode(du,u,p,t)
    S,I,R,B,G,D = u
    g,c,k,n,m,b,q = p
    du[1] = -B*S*I
    du[2] = B*S*I-g*I
    du[3] = g*I
    du[4] = b*(G + g) / n*(1 - c)
    du[5] = -q*G*D
    du[6] = k*D*(1 - (D/m))
end

params = [0.07,0.05,0.1,1000.0,10.0,0.01,0.0005]
init = [999.0,1.0,0.0,6.8e-5,0.14,0.1]
tspan = (0.0,200.0)
sir_prob = ODEProblem(sir_ode,init,tspan,params)
sir_sol = solve(sir_prob,saveat = 0.1,alg_hints=[:stiff])

plt = plot(sir_sol,vars=[(0,1), (0,2), (0,3)])
#plt = plot(sir_sol,vars=[(0,4)])
gui(plt)