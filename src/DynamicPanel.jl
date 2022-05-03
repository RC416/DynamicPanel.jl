module DynamicPanel

export static_panel

# Packages.
using LinearAlgebra, Kronecker, DataFrames, CSV, Distributions, InvertedIndices


### Static Panel function
#=
Data: dataframe with observables
y_var: column name for dependent variable
x_var: vector of column names for independent variables
individual index: column name for variable that identifies individuals
period_index: column name for variable that identifies periods
method: one of ["pooling", "within", "Mundlak CRE", "Chamberlain CRE"]
balanced_panel: [0, 1], 1:= balanced panel, 0:= unbalanced
period_effects: [0, 1], 1:= include period dummies, 0:= no period dummies
=#

function static_panel(data; y_var, x_vars, individual_index, period_index, method, balanced_panel=0, period_effects=1)

    ### 1 - Prepare data.
    # Sort data by individual index (inner) and period index (outer)
    sort!(data, [period_index, individual_index])

    # Get vectors/matrices of observables
    Y = Vector(data[:,y_var])
    X = Matrix(data[:,x_vars])
    i = Vector(data[:,individual_index])
    t = Vector(data[:,period_index])

    # Get individual, period and other useful values.
    nᵢ = size(unique(i),1)
    nₜ = size(unique(t),1)
    t_min = minimum(t)
    t_max = maximum(t)
    n_covariates = size(x_vars,1)
    i = i .- (t_min - 1)             # Standardize time values to start at 1.

    # Optional: address unbalanced panel.
    if balanced_panel != 1
        # Create a "full" vector of individual and period indices as if panel was balanced / no missing observations.
        t_full = kron(t_min:t_max, ones(nᵢ))
        i_full = repeat(unique(i), t_max-t_min+1)
        n_full = size(t_full,1) # number of obserations if panel had been balanced.

        # Create "full" vector of other observables.
        Y_full = zeros(n_full)
        X_full = zeros(n_full, size(X,2))

        # Get index of rows (individual-period) where observation exists
        obs_ind = findall(in(i+im*t), (i_full+im*t_full))
        missing_ind = findall(!in(i+im*t), (i_full+im*t_full))
        # Note: i+im*t creates a unique imaginary number for each individual-period combination.

        # Fill in values where observations exist. 
        Y_full[obs_ind] = Y
        X_full[obs_ind,:] = X

        # Reassign to main variable names.
        Y = Y_full
        X = X_full 
        i = i_full
        t = t_full
    end

    # Optional: include period dummies.
    if period_effects == 1
        # Dummy variable matrix. 
        E = I(nₜ) ⊗ ones(nᵢ) # equivalent to: kron(I(nₜ),ones(nᵢ))
        # Include in X matrix.
        X = [X E]
        # Drop dummies from any rows with missing observations.
        if balanced_panel != 1
            X[missing_ind,:] .= 0
        end 
    end

    ### 2 - Perform estimation according to selected method. 
    # A. Pooling: ignore individual effects and perform OLS.
    if method == "pooling"
        # Estimate parameters and variance-covarnaice matrix.
        β = (X'X)\X'Y
        vcov = inv(X'X)*(X'Diagonal((Y-X*β).^2)*X)*inv(X'X)
        SE = sqrt.(diag(vcov))

        # Display results.

        # Save/export results.
        return (coefs = β, vcov = vcov, SE = SE, residuals = Y-X*β);
    end

    # B. Within transformation
    if method == "within"
        # Individual dummies.
        D = repeat(I(nᵢ), nₜ)

        # Estimate parameters using FWT theorem method.
        DᵀD⁻¹ = 1/nₜ * I(nᵢ)                    # precalculate large inverse term by exploiting structure of D
        êx = (I(size(X,1)) - D*(DᵀD⁻¹)*D')*X   # = xᵢₜ - x̄ᵢ , repeated for each period
        êy = (I(size(Y,1)) - D*(DᵀD⁻¹)*D')*Y   # = yᵢₜ - ȳᵢ , repeated for each period
        β = (êx'êx)\êx'êy
        vcov = inv(X'X)*(X'*Diagonal((êy - êx*β).^2)*X)*inv(X'X)
        SE = sqrt.(diag(vcov))

        return (coefs = β, vcov = vcov, SE = SE, residuals = Y-X*β);
    end

end

end
