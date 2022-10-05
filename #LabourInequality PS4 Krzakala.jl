#LabourInequality PS4 Krzakala
#Task 2
using Plots
using Statistics
using DataFrames
using CSV

cd("C:\\Users\\mpkrz\\OneDrive - SGH\\Dokumenty\\GitHub\\pslabour4")

const sigma = 1.5
const beta = 0.993422
const q = 0.998
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
global p = zeros(length(grid_A),length(grid_e))

max_iter = 2000
tolerance = 1e-5

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
            p[idx_a, idx_e]=grid_A[bx]
        end
    end
    dev = maximum(abs.(TV-V))
    global V = deepcopy(TV)

    if dev <= tolerance
        break
    end
end

plot(grid_A, p[:,1])
plot!(grid_A, p[:,2])
plot!(grid_A, grid_A)
savefig("polfunc.png")

#Task3

Transition=zeros(200,200)
ae_dist=zeros(1,200)


for (idx_e_r, e) in enumerate(grid_e)
    for (idx_e_c, e) in enumerate(grid_e)
        for i in 1:100
            Transition[floor(Int, i+100*(idx_e_r-1)),floor(Int, 100*(idx_e_c-1)+g[i,idx_e_c])]=pi[idx_e_r,idx_e_c]
        end
    end
end
df=DataFrame(Transition, :auto)
CSV.write("transition_m.csv", df)

#= Transition=zeros(200,200)
ae_dist=ones(1,200)
for (idx_s, s) in enumerate(ae_dist)
    for (idx_s2, s) in enumerate(ae_dist)
        pi_r=(idx_s>100)+1
        pi_c=(idx_s2>100)+1
        a_idx=mod(idx_s,100)
        if a_idx == 0
            a_idx=100
        end
        if g[a_idx,pi_r]==idx_s2
            Transition[idx_s, idx_s2]=pi[pi_r,pi_c]
        end
    end
end
 =#
 #using CSV, DataFrames


 ae_dist_init=ones(1,200)
ae_dist_init=ae_dist_init/200
ae_dist_new=deepcopy(ae_dist_init)

for i=1:99999
    ae_dist_new=ae_dist_init*Transition
    diff=maximum(abs.(ae_dist_init-ae_dist_new))
    ae_dist_init=ae_dist_new
    if diff<0.00001
        break
    end
end    

sum(ae_dist_init)
plot(grid_A,ae_dist_init[1,1:100],label="e_H")
plot!(grid_A,ae_dist_init[1,101:200],label="e_L")
plot!(grid_A,(ae_dist_init[1,1:100]+ae_dist_init[1,101:200]),label="sum")

avg_a=mean(ae_dist_init)
avg_a_idx=findfirst(index->index>=avg_a,grid_A)
grid_A[avg_a_idx]
avg_a
grid_A[avg_a_idx-1]

#Part D
ae_dist_sim=zeros(2,1000)
ae_dist_sim[1,1:500].=avg_a+1
ae_dist_sim[1,501:1000].=avg_a+0.1
ae_dist_sim[2,1:500].=1
ae_dist_sim[2,501:1000].=0
ae_dist_sim_init=deepcopy(ae_dist_sim)
ae_dist_hist=zeros(1001,1000)
ae_dist_hist[1,1:1000]=ae_dist_sim_init[1,1:1000]
ae_dist_sim_mean=zeros(1001)
ae_dist_sim_mean[1]=avg_a+(0.1+1)/2
function findclose(to_find,lin_range)
    c=lin_range.-to_find
    d=collect(c)
    e=broadcast(abs, d)
    findmin(e)[2]
end

for t = 1:1000
    print(t)
    random=rand(1,1000)
    for i = 1:1000
        if random[i]<pi[(floor(Int,ae_dist_sim_init[2,i]+1)),1] ae_dist_sim[2,i]=1 else (ae_dist_sim[2,i]=0) end
        if ae_dist_sim[2,i]==1 e=1 else e=0.1 end
        ae_dist_sim[1,i]=p[findclose(ae_dist_sim_init[1,i],grid_A),floor(Int,(ae_dist_sim_init[2,i]+1))]+e
    end
    ae_dist_hist[t+1,1:1000]=ae_dist_sim[1,1:1000]
    ae_dist_sim_mean[t+1]=mean(ae_dist_sim[1,1:1000])
    ae_dist_sim_init=deepcopy(ae_dist_sim)
end
plot(ae_dist_sim_mean)
savefig("plot_mean.png")
plot(ae_dist_hist[:,2])

