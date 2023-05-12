#Random walk attempt
using Random
using Plot
#Initialise the variance of the drift
phi_sigma = 3
iter = 1000
monte = 50
phi_sigmama = randn!(zeros(iter,monte))

# u can never be negative so if u is <= 1, and sigmama is <0, flip sigmama
    
for i in 1:iter
    for j in 1:monte
       phi_sigmama[i,j] = phi_sigma * phi_sigmama[i,j] 
    end
end

n_sigma = 2
n_sigmama = randn!(zeros(iter, monte))
for i in 1:iter
    for j in 1:monte
        n_sigmama[i,j] = n_sigma * n_sigmama[i,j]
    end
end

beta = zeros(iter, monte)
#initialise

for i in 1:iter-1
    beta[i,1] = 0.09
    for j in 1:monte
        if beta[i,j] < 0 && phi_sigmama[i+1,j] < 0
            beta[i+1,j] = beta[i,j] - phi_sigmama[i+1,j]
        else
            beta[i+1,j] = beta[i,j] + phi_sigmama[i+1,j]
        end
    end
end
    
u = zeros(iter, monte)
#initialise
u[1,1:monte] .= 27
for i in 1:iter-1
    for j in 1:monte
        u[i+1,j] = u[i,j] + beta[i,j] + n_sigmama[i,j]
    end
end
p1 = plot()
for i in 1:monte
    plot!(p1, 1:iter, u[1:iter,i], legend = false, linecolour = "black")
end
# p1 = plot()
# plot!(1:iter, beta[1:iter])
# plot!(1:iter, u[1:iter])
display(p1)
