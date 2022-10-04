#LabourInequality PS4 Krzakala
#Task 2
using Plots

const sigma = 1.5
const beta = 0.993422
const q=0.998
global pi=[0.925 0.075;0.5 0.5]

function u(c)
    if sigma==1
        log(c)
    else
        c^(1-sigma)/(1-sigma)
    end
end



grid_A = LinRange(-4, 4, 100)
grid_e = [1,0.1]


# Initiate the value function guess
global V = zeros(length(grid_A),length(grid_e))
global g = zeros(length(grid_A),length(grid_e))

max_iter = 2000
tolerance = 0.00000001

for iter = 1:max_iter
    global TV = zeros(length(grid_A),length(grid_e))

    for (idx_e, e) in enumerate(grid_e)
        for (idx_a, a) in enumerate(grid_A)
            U = zeros(length(grid_A))
            for (idx_ap, ap) in enumerate(grid_A)
                c = a + e - q*ap
                if c>0
                    U[idx_ap] = u(c)+pi[idx_e,1]*beta*V[idx_ap, 1]+pi[idx_e,2]*beta*V[idx_ap,2]
                else
                    U[idx_ap] = -9e8
                end
            end

            ax, bx = findmax(U)
            TV[idx_a,idx_e]=ax
            g[idx_a, idx_e]=bx
        end
    end
    dev = maximum(abs.(TV-V))
    global V = deepcopy(TV)

    if dev <= tolerance
        break
    end
end

plot(grid_A, g[:,1])
plot!(grid_A, g[:,2])
plot!(grid_A, grid_A)