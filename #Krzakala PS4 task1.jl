#LabourInequality PS4 Krzakala
#Task 1
using Statistics
n_h=1
n_l=1

n=n_h+n_l
t_skip=150
t_interest=480
t_total=t_skip+t_interest
year=Int(t_interest/6)

pi=[0.925 0.075;0.5 0.5]

#states of all individuals in array:

s=zeros(n,t_total)
s[1:n_h,1]=ones(n_h,1)

sz = size(s)

random=rand(sz[1],sz[2])

for t = 2:t_total
    for i = 1:n
        if s[i,t-1]==1 && random[i,t-1] < pi[1,1]
            s[i,t]=1
        elseif s[i,t-1]==0 && random[i,t-1]< pi[2,1]
            s[i,t]=1
        end
    end
end

skipped=s[:,(t_skip+1):t_total]
#to endowments
e=deepcopy(skipped)
e[e.==0] .= 0.1

#calculate yearly endowments of people starting at 1 and 0.1

rich_yr=zeros(1,year)

for yr=1:year
    rich_yr[1,yr]=sum(e[1,(yr):(yr+6)])
end

avg_e_rich=mean(rich_yr)
std_e_rich=std(rich_yr)
stdm_e_rich=std_e_rich/avg_e_rich*100

poor_yr=zeros(1,year)

for yr=1:year
    poor_yr[1,yr]=sum(e[2,(yr):(yr+6)])
end

avg_e_poor=mean(poor_yr)
std_e_poor=std(poor_yr)
stdm_e_poor=std_e_poor/avg_e_poor*100

n_negative_shocks_rich=sum((diff(skipped[1,1:t_interest])).==(-1))
n_negative_shocks_poor=sum((diff(skipped[2,1:t_interest])).==(-1))


avg_dur_low_rich=(t_interest-sum(skipped[1,1:t_interest]))/n_negative_shocks_rich
avg_dur_low_poor=(t_interest-sum(skipped[2,1:t_interest]))/n_negative_shocks_poor