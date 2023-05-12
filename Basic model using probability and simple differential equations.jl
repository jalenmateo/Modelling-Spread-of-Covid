#Imad A. Moosa (2020) The effectiveness of social distancing in containing Covid-19, Applied Economics, 52:58, 6292-6305, DOI: 10.1080/00036846.2020.1789061

using Plots

#plot R against f
#Initialising values

R_0 = 3
a = [0.25, 0.5, 0.75, 1]
f = [x for x in range(0.0, stop = 1, step = 0.01)]
max_iter = length(f)
numberof = length(a)

function r(a,f,R_0)
    ( 1 - ( 1 - a^2 ) * f ) * R_0
end

R = zeros(max_iter, numberof)

for i in 1:numberof
    for j in 1:max_iter
        R[j,i] = r(a[i], f[j], R_0)
    end
end

p1 = plot(xlims = (0,1), ylims=(0,3))
for i in 1:numberof
    plot!(p1, f[1:max_iter], R[:,i])
end
plot!()

#plot of R against a

R_0 = 3
f = [0.25, 0.5, 0.75, 1]
a = [x for x in range(0.0, stop = 1, step = 0.01)]
max_iter = length(a)
numberof = length(f)

function r(a,f,R_0)
    ( 1 - ( 1 - a^2 ) * f ) * R_0
end
R = zeros(max_iter, numberof)

for i in 1:numberof
    for j in 1:max_iter
        R[j,i] = r(a[j], f[i], R_0)
    end
end

p1 = plot(xlims = (0,1), ylims=(0,3))
for i in 1:numberof
    plot!(p1, a[1:max_iter], R[:,i])
end
plot!()

#partial derivative of R wrt f

R_0 = 3
f = [0.25, 0.5, 0.75, 1]
a = [x for x in range(0.0, stop = 1, step = 0.01)]
max_iter = length(a)
numberof = length(f)

function rf(R_0, a)
    ( a^2 - 1 ) * R_0
end
Rf = zeros(max_iter)

for j in 1:max_iter
        Rf[j] = rf(R_0, a[j])
end

plot(a[1:max_iter] ,Rf[1:max_iter])
