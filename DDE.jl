#Delayed Differential Equations

#With reference from ISSN: 2248-9622 , A Delay Differential Equation Model for the Spread of COVID-19

using Plots
#rate of emergence = spread of existing * prob of susceptibility * existing cases

#To modify this model, we can add on to it by tweaking m0 to have components related to social distancing and mask wearing.
#I focused on social distancing and mask wearing for this project

maxiter = 1000  #time
y = zeros(maxiter)  #number infected

y[1] = 3  #initial conditions for infected people
y[2] = 3
y[3] = 5
y[4] = 5
y[5] = 5
y[6] = 9

N =  6000000 #susceptible people 
# t = time
mu1 = 0.1  #fraction of asymptomatic cases
t1 = 3  #asymptomatic infectious time
t2 = 5  #incubation time/latency time
t2_ = 3  #time remained at large/not caught(t2/2)
t3 = 6 # recovery time
mu3 = 0.2  #fraction of people escaped from contact tracing
m0 = 0.5 #rate at which the at large cases transmits the disease to other people, assuming all are susceptible
# 3 classes at large - 1-mu3 [contact traced cases], mu3(1-mu1) [untraced symptomatic] and mu3mu1 [untraced asymptomatic]

n_iter = 100  #range to plot
n_iter = maxiter-n_iter + 1

#### To plot number of infected against time
for i in 7:maxiter
    y[i] = y[i-1] + m0*(1-y[i-1]/N)*((y[i-1])-(1-mu3)*y[i-1-t2_]-(1-mu1)*mu3*y[i-1-t2]-mu1*mu3*y[i-1-t1])
end

new_y = zeros(maxiter)
for i in 1:6
    new_y[i] = y[i]
end
for i in 7:maxiter
    new_y[i] = y[i] - y[i-t3]
end
plot(1:maxiter, new_y[1:maxiter], title = "Simulation using a DDE Model", linecolour = "turquoise", linewidth = 2, yaxis = "Number of people infected", xaxis = "Time (Days)", legend = false)
#### To plot total number of infected against time
plot(1:maxiter, y[1:maxiter])
####

