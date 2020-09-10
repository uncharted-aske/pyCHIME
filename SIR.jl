using DifferentialEquations
using Plots

function sir_ode(du,u,p,t)
    S,I,R = u
    b,g = p
    du[1] = -b*S*I
    du[2] = b*S*I-g*I
    du[3] = g*I
end

params = [0.1,0.05]
init = [0.999,0.001,0.0]
tspan = (0.0,200.0)
sir_prob = ODEProblem(sir_ode,init,tspan,params)
sir_sol = solve(sir_prob,saveat = 0.1)

plt = plot(sir_sol)
gui(plt)