# Static Discrete Choice exercise
function sdc(N, numSch;
                     γ = 1.2,
                     σ = sqrt((pi^2)/6.0),
                     seed = 0)

#seed = 0
#N=500; #number of students
#numSch=5; #number of schools
#γ    = 1.2;
#σ    = sqrt((pi^2)/6.0) #I am fixing the variance of the normal to be equal of the T1EV(0,1)

# Generating the data
# Main parameters:
Random.seed!(seed+1);
δ    = randn(numSch,1); #generating the parameter deltas (mean utility of each school)
dist = abs.(randn(N,numSch)); #generating the observable distance

# Draw the utilities
Random.seed!(seed+2);
ϵ=σ*randn(N,numSch);
util = repeat(δ',N,1) .- γ*dist .+ ϵ;

#Let's introduce a variation in the menu
menuSize = zeros(Int64,N,1);
choice   = zeros(Int64,N,1);
menu     = zeros(Bool, N,numSch);

#Loop
Random.seed!(seed+3);
for i=1:N
    #Number of schools in the menu
    menuSize[i]=rand(DiscreteUniform(1,numSch-1)) + 1;
    # Generate a random integer uniformaly distributed between 1 and (numSch-1)
    # I add 1 to assure that they will have at least 2 options

    #Identity of schools in the menu
    temp=randperm(numSch);
    # return a random permutation of integers between 1 and NumSch

    for j=1:numSch
        if any(temp[1:menuSize[i]].==j)
            # any() returns equal 1 if any element of the array is non-zero
            # 1:menuSize[i] is the menu of individual i
            # Thus I will put one in the matrix menu when the i have the
            # option j available
            menu[i,j]=true;
        end
    end

    #Choose the school that provides highest utility
    choice[i]=findall(util[i,:].== findmax(util[i,menu[i,:]])[1])[1]; # I am assuming no ties!
    # I am returning the position of the option with the highest utility
    # that is available in the menu

end

# Let's compute the likelihood funtion
# We will assume ϵ∼T1EV(0,1)
function dcobj_t1ev(x,choice,dist,N,numSch,menu)
    δ = [x[1:numSch-1];0.0];
    γ = x[numSch];
    v = exp.(δ'.-γ*dist).*menu;
    PredProb = v ./ repeat(sum(v,dims=2),1,numSch);
    chocProb = [PredProb[i,choice[i]] for i in 1:N];
    objVal   =-sum(log, chocProb);

    return objVal
end

#Now lets assume that  ϵ∼N(0,π^2/6)
#In this case we have to simulate the distribution in order to get the likelihood function.
function dcobj_normal(x,choice,dist,N,numSch, menu,seed)
    δ = [x[1:numSch-1];0.0];
    γ = x[numSch];
    σ = sqrt((pi^2)/6);

    # Draw the epsilon to simulate the choice probabilities
    Random.seed!(seed+4);
    numDraws=numSch*20;
    ϵsim = σ*randn(numDraws,numSch, N);

    #Compute the simulated utilities and the choice for each draws
    # Note that here we need only the probability of the best option, the other
    # probabilities do not matter here
    simChoice = zeros(Int64,numDraws,1);
    simProb   = zeros(N,1);

    for i=1:N
        simUtil = repeat(δ',numDraws,1) - γ*repeat(dist[i,:]',numDraws,1) + ϵsim[:,:,i];
        for k=1:numDraws
            simChoice[k]=findall(simUtil[k,:].== findmax(simUtil[k,menu[i,:]])[1])[1];
        end
        simProb[i] = sum(simChoice.==choice[i])/numDraws;
    end

    objVal   =-sum(log, simProb);
    return objVal
end

# Let's call the optimizer
xtrue =     [δ[1:numSch-1].-δ[numSch];γ];
xinit = 0.9*[δ[1:numSch-1].-δ[numSch];γ];
#dcobj_t1ev(xinit,choice,dist,N,numSch,menu)
#dcobj_normal(xinit,choice,dist,N,numSch,menu,seed)

out_t1ev   = optimize(x -> dcobj_t1ev(x,choice,dist,N,numSch, menu) , xinit)
out_normal = optimize(x -> dcobj_normal(x,choice,dist,N,numSch,menu,seed), xinit)

println("ϵ-T1EV:    minimum = $(out_t1ev.minimum) with argmin = $(out_t1ev.minimizer) in $(out_t1ev.iterations) iterations")
println("ϵ-Normal:  minimum = $(out_normal.minimum) with argmin = $(out_normal.minimizer) in $(out_normal.iterations) iterations")
println("θtrue = $(xtrue), θhat_t1ev = $(out_t1ev.minimizer), θhat_normal = $(out_normal.minimizer)")

# Organizing the outputs
comparison = (θtrue = xtrue, θhat_t1ev = out_t1ev.minimizer, θhat_normal = out_normal.minimizer);
data = ( δ=δ, dist=dist, ϵ=ϵ, choice=choice, menu=menu);
params = (γ = γ, σ = σ);
output = (comparison = comparison, data = data, params=params);

return output;
end
