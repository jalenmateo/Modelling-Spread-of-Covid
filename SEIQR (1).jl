# With reference to DOI: 10.31579/CRCT.2020/006 Differential Equation Analysis on COVID-19
using Plots

si = 1.3
gamma = 0.99
theta = 0.1

B = 0.5 #contagion rate
A = 0.5 #removal rate for quarantine
nu = 0.5 #removal rate for latent
N = 54500000 #Population
v = 5 #incubation period
phi = 1.25 #(1/phi) $$average delayed reporting period
w = 1 #infection reucing factors in the exposed infectious
lamb = 0.8 #infection reducing factors in the latent infections

p_mask = 0.84     #0.84

avg_new = 14
avg_old = 30

eps = 0.1
tau = 0.09

itr = 1000

Sus = zeros(itr)
Sus[1] = N-100

Exp = zeros(itr)
Exp[1] = 100

Infc = zeros(itr)
Infc[1] = 5

Q = zeros(itr)
Q[1] = 3000

R = zeros(itr)
R[1] = 1000

L = zeros(itr)
L[1] = 0

p = [w, B, lamb, N, tau, si, v, phi, eps, A, nu, gamma, p_mask, theta, avg_old, avg_new]
# print(p)


for i in 2:itr
    Sus[i] = Sus[i-1] - B * Sus[i-1] * (w * Exp[i-1] + Infc[i-1] + Q[i-1]) / N
    Exp[i] = Exp[i-1] + (- eps * tau * Exp[i-1] +  B * Sus[i-1] * (w * Exp[i-1] + Infc[i-1] + Q[i-1]) / N)*(1+theta*((-avg_old+avg_new)/avg_old)-si*(p_mask-0.5))
    Infc[i] = Infc[i-1] + eps * tau * Exp[i-1] - A * (1-v) * Infc[i-1] - v * Infc[i-1] / phi
    Q[i] = Q[i-1] + (v * Infc[i-1] / phi - A * Q[i-1])*gamma
    R[i] = R[i-1] + A * Q[i-1] + A * (1-v) * Infc[i-1]
end
# print(Infc[1:100])
plot(1:itr, Infc[1:itr], title = "Simulation Using an SEIQR Model", linewidth = 2 ,linecolour = "turquoise", yaxis = "Number of people infected", xaxis = "Time (Days)", legend=false)

#Alternatively, use the SEIQR package in julia
using DifferentialEquations
u0 = [N-30, 25, 5, 0, 0]

function seiqr!(du, u, p,t) 
    du[1] =  - p[2] * u[1] * (p[1] * u[2] + u[3] + u[4]) / p[4]
    du[2] = (1 + p[14]*(p[15] - p[16]) + (p[13] - 0.5)*p[6]) - p[9] * p[5] * u[2] +  p[2] * u[1] * (p[1] * u[2] + u[3] + u[4]) / p[4]
    du[3] = p[9] * p[5] * u[2] - p[10] * (1-p[7]) * u[3] - p[7] * u[3] / p[8]
    du[4] = p[12] * (p[7] * u[3] / p[8] - p[10] * u[4])
    du[5] = (p[10] * u[4] + p[10] * (1-p[7]) * u[3])
end

p = [w, B, lamb, N, tau, si, v, phi, eps, A, nu, gamma, p_mask, theta, avg_old, avg_new]
tspan = (0,320)
prob = ODEProblem(seiqr!,u0,tspan,p)
sol = solve(prob)
# sol1 = sol[1]
# plot(sol1)
plot(sol, layout =(5,1), label = ["S" "E" "I" "Q" "R"], yaxis = "", linecolour = "turquoise", size = (700,700), linewidth = 3, ylabel = "Ppl", xlabel = "time (days)")
