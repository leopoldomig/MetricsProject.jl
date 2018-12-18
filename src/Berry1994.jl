# Monte Carlo Simulation of Berry 1994

# Selecting the needed packages
function Berry1994(R, repNum;
                     βo = 5.00,
                     βx = 2.00,
                     α  = 1.00,
                     γo = 1.00,
                     γx = 0.50,
                     γw = 0.25,
                     σc = 0.25,
                     σw = 0.25,
                     σd = 3.00)



#Let's set up the remaining parameters
J  = 2

# Let's randomly draw observable and unobservable characteristics that affect
# demand: xₘⱼ and ξₘⱼ
xₘⱼ =  randn(R, J, repNum);
ξₘⱼ =  randn(R, J, repNum);

# Let's randomly draw observable and unobservable characteristics that affect
# marginal cost: wₘⱼ and ωₘⱼ
wₘⱼ =  randn(R, J, repNum);
ωₘⱼ =  randn(R, J, repNum);

# Let's compute the marginal cost (cₘⱼ) for each market m and product/firm j
# as Berry does:
cₘⱼ = exp.( γo.*ones(R,J,repNum) .+ γx.*xₘⱼ + σc.*ξₘⱼ .+ γw.*wₘⱼ .+ σw.*ωₘⱼ );

# Let's compute the equilibrium prices (pₘⱼ) and market share (sₘⱼ)
function eq_p_sh( xₘ, ξₘ, cₘ, sₘ)
    p = cₘ .+ ( α.*(1.0.-sₘ)).^(-1);
    v = βo .+ βx.*xₘ .+ σd.*ξₘ .- α.*p;
    F = (exp.(v))./(1.0 .+ sum(exp.(v)));
    return F
end

s₀ = [0.5, 0.5];
sₘⱼ = similar(xₘⱼ);
pₘⱼ = similar(xₘⱼ);

for k = 1:repNum
    for r = 1:R
        sol = nlsolve(y -> eq_p_sh( xₘⱼ[r,:,k], ξₘⱼ[r,:,k], cₘⱼ[r,:,k], y)-y, s₀)

        s_eq = sol.zero

        if ( sol.f_converged == false || sum(s_eq .< 0)>0 || sum(s_eq)>1 )
            s_eq = [NaN , NaN];
        end

        p_eq = cₘⱼ[r,:,k] .+ (α.*(1.0 .- s_eq) ).^(-1);

        sₘⱼ[r,:,k] = s_eq;
        pₘⱼ[r,:,k] = p_eq;

    end
end

@assert sum(isnan.(sₘⱼ)) == 0

# Let's calculate the mean utility of product j: (δₘⱼ)
δₘⱼ = similar(xₘⱼ);
for k = 1:repNum
for r=1:R
    δₘⱼ[r,1,k] = log(sₘⱼ[r,1,k])-log(1-sₘⱼ[r,1,k]-sₘⱼ[r,2,k]);
    δₘⱼ[r,2,k] = log(sₘⱼ[r,2,k])-log(1-sₘⱼ[r,1,k]-sₘⱼ[r,2,k]);
end
end

# Let's compute the OLS estimates
θhatᵒˡˢ = zeros(3,repNum);

for k = 1:repNum
    y = reshape(δₘⱼ[:,:,k],R*J,1);
    X = [ones(R*J,1)  reshape(xₘⱼ[:,:,k],R*J,1) reshape(pₘⱼ[:,:,k],R*J,1) ];
    θhatᵒˡˢ[:,k] = inv(X'*X)*(X'*y);
end

# Taking the average and the standard deviation of the simulations
m_θhatᵒˡˢ  = mean(θhatᵒˡˢ,dims=2);
sd_θhatᵒˡˢ = std(θhatᵒˡˢ,dims=2);

# Let's compute the IV estimates
θhatⁱᵛ = zeros(3,repNum);

aux_xₘⱼ = similar(xₘⱼ);
aux_xₘⱼ[:,1,:] = xₘⱼ[:,2,:];
aux_xₘⱼ[:,2,:] = xₘⱼ[:,1,:];

for k = 1:repNum
    y = reshape(δₘⱼ[:,:,k],R*J,1);
    X = [ones(R*J,1)  reshape(xₘⱼ[:,:,k],R*J,1)  reshape(pₘⱼ[:,:,k],R*J,1) ];
    Z = [ones(R*J,1)  reshape(xₘⱼ[:,:,k],R*J,1)  reshape(wₘⱼ[:,:,k],R*J,1)  reshape(aux_xₘⱼ[:,:,k],R*J,1) ];
    Xhat = Z*inv(Z'*Z)*(Z'*X);
    θhatⁱᵛ[:,k] = inv(Xhat'*Xhat)*(Xhat'*y);
end

# Taking the average and the standard deviation of the simulations
m_θhatⁱᵛ  = mean(θhatⁱᵛ, dims=2);
sd_θhatⁱᵛ = std(θhatⁱᵛ, dims=2);

# Organizing the outputs
θ = [βo ; βx ; -α];
comparison = (θ = θ, m_θhatᵒˡˢ= m_θhatᵒˡˢ,  sd_θhatᵒˡˢ=sd_θhatᵒˡˢ,  m_θhatⁱᵛ= m_θhatⁱᵛ,  sd_θhatⁱᵛ= sd_θhatⁱᵛ);
data = ( xₘⱼ=xₘⱼ, ξₘⱼ=ξₘⱼ, wₘⱼ=wₘⱼ, ωₘⱼ=ωₘⱼ, cₘⱼ=cₘⱼ, sₘⱼ=sₘⱼ, pₘⱼ=pₘⱼ, δₘⱼ=δₘⱼ );
params = (βo = βo, βx = βx, α  = α, γo = γo, γx = γx, γw = γw, σc = σc, σw = σw, σd = σd);
output = (comparison = comparison, data = data, params=params);

return output

end
