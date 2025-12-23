# Symmetric Map Simulation in Julia with different parameters for A and B sides

using Plots
using LinearAlgebra
using LaTeXStrings
using DelimitedFiles

# f_A: Solve the implicit equation for R_A: R_A = 1 - exp(-c_A * p_A * R_A)
function f_A(p_A, c_A; tol=1e-8, maxiter=1000)
    R = 0.99  # initial guess (nonzero)
    for i in 1:maxiter
        R_new = 1 - exp(-c_A * p_A * R)
        if abs(R_new - R) < tol
            return R_new
        end
        R = R_new
    end
    return R  # if not converged, return last iterate
end

# f_B: Solve the implicit equation for R_B: R_B = 1 - exp(-c_B * p_B * R_B)
function f_B(p_B, c_B; tol=1e-8, maxiter=1000)
    R = 0.99  # initial guess (nonzero)
    for i in 1:maxiter
        R_new = 1 - exp(-c_B * p_B * R)
        if abs(R_new - R) < tol
            return R_new
        end
        R = R_new
    end
    return R
end

# Regulations: both interlayer and intralayer regulations

function g_A(R_A, R_B, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter)
    
    if c_A_plus_intra == Inf
        p_intra = exp(-c_A_minus_intra * R_A)
    else
        p_intra = (1 - exp(-c_A_plus_intra * R_A)) * exp(-c_A_minus_intra * R_A)
    end

    if c_B_plus_inter == Inf
        p_inter = exp(-c_B_minus_inter * R_B)
    else
        p_inter = (1 - exp(-c_B_plus_inter * R_B)) * exp(-c_B_minus_inter * R_B)
    end

    return p * p_intra * p_inter
end

function g_B(R_A, R_B, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra)

    if c_A_plus_inter == Inf
        p_inter = exp(-c_A_minus_inter * R_A)
    else
        p_inter = (1 - exp(-c_A_plus_inter * R_A)) * exp(-c_A_minus_inter * R_A)
    end

    if c_B_plus_intra == Inf
        p_intra = exp(-c_B_minus_intra * R_B)
    else
        p_intra = (1 - exp(-c_B_plus_intra * R_B)) * exp(-c_B_minus_intra * R_B)
    end

    return p * p_intra * p_inter
end

# -----------------
# Parameter settings
# -----------------
p = 0.72        # common parameter in g functions

# Parameters for f_A and f_B
c_A = 30       # parameter in f_A for R_A
c_B = 30       # parameter in f_B for R_B

# Parameters in g_A (affecting p_A, which is computed from R_B)
c_B_plus_intra = Inf
c_B_minus_intra = 0
c_B_plus_inter = 3
c_B_minus_inter = 0

c_A_plus_intra = 5
c_A_minus_intra = 1.5
c_A_plus_inter = Inf
c_A_minus_inter = 3



# -----------------
# Simulation settings
# -----------------
T = 20000  # number of iterations

# Preallocate arrays to store iterates for R_A and R_B
R_A = zeros(T)
R_B = zeros(T)

# Initial conditions:
# We choose R_A[1] and R_B[1] freely.

R_A[1] = 0.2  # initial guess for R_A
R_B[1] = 0.8

# -----------------
# Main iteration loop
# -----------------
for t in 2:T
    # Compute p_A^t and p_B^t based on previous values:
    p_A = g_A(R_A[t-1], R_B[t-1], p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter)
    p_B = g_B(R_A[t-1], R_B[t-1], p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra)

    # Solve for the next R_A and R_B using the implicit equations:
    R_A[t] = f_A(p_A, c_A)
    R_B[t] = f_B(p_B, c_B)
end

# h_A(R_B) = f_A(g_A(R_B))
function h_A(R_A, R_B, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter, c_A)
    p_A = g_A(R_A, R_B, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter)
    return f_A(p_A, c_A)
end

# h_B(R_A) = f_B(g_B(R_A))
function h_B(R_A, R_B, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra, c_B)
    p_B = g_B(R_A, R_B, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra)
    return f_B(p_B, c_B)
end


function orbit_diagram()
    R_A_list = []
    R_B_list = []
    p_list = []
    real_part = []
    imag_part = []
    collecting = true
    R_A_default = 0.2
    R_B_default = 0.6
    for p = 0.735:-0.0002:0.68
        collecting = true
        R_A = zeros(T)
        R_B = zeros(T)
        pp = zeros(T)
        R_A[1] = R_A_default  # initial guess for R_A
        R_B[1] = R_B_default
        for t in 2:T
            # Compute p_A^t and p_B^t based on previous values:
            p_A = g_A(R_A[t-1], R_B[t-1], p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter)
            p_B = g_B(R_A[t-1], R_B[t-1], p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra)

            # Solve for the next R_A and R_B using the implicit equations:
            R_A[t] = f_A(p_A, c_A)
            R_B[t] = f_B(p_B, c_B)
            pp[t] = p
           
        end
        R_A_default = R_A[T]
        R_B_default = R_B[T]
        # if R_A[T] != 0 && R_B[T] != 0
        #     R_A_default = R_A[T]
        #     R_B_default = R_B[T]
        # else
        #     R_A_default = 0.2
        #     R_B_default = 0.6
        # end

        # Store the last 100 points of the orbit
        R_A_list = [R_A_list; R_A[T-10000:T]]
        R_B_list = [R_B_list; R_B[T-10000:T]]
        p_list = [p_list; pp[T-10000:T]]

        # eigenvalues = compute_eigenvalues(R_A[T], R_B[T], p,
        #     c_A_plus_intra, c_A_minus_intra, c_A_plus_inter, c_A_minus_inter,
        #     c_B_plus_intra, c_B_minus_intra, c_B_plus_inter, c_B_minus_inter,
        #     c_A, c_B)
        # rea = real(eigenvalues)[1]
        # ima = imag(eigenvalues)[1]
        # if sqrt(rea^2 + ima^2) > 1 && collecting
        #     collecting = false
        # end
        # if collecting
        #     push!(real_part, real(eigenvalues)[1])
        #     push!(imag_part, imag(eigenvalues)[1])
        # end
    end
    # Plot the orbit diagram
    return R_A_list, R_B_list, p_list, real_part, imag_part
end

"""
    compute_eigenvalues(R, params; ε=1e-6)

Computes the eigenvalues of the Jacobian of the map defined by h_A and h_B at the point R = [R_A, R_B].

- `R` is a two-element vector, e.g. [R_A*, R_B*].
- `params` is a tuple containing the parameters required by h_A and h_B.
- `ε` is the finite-difference step size (default 1e-6).
"""
function compute_eigenvalues(R_A, R_B, p, 
    c_A_plus_intra, c_A_minus_intra, c_A_plus_inter, c_A_minus_inter, 
    c_B_plus_intra, c_B_minus_intra, c_B_plus_inter, c_B_minus_inter,
    c_A, c_B; ε=1e-6)

    # Define central difference approximations for the partial derivatives.
    d_hA_d_RA = (h_A(R_A + ε, R_B, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter, c_A) - h_A(R_A - ε, R_B, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter, c_A)) / (2ε)
    d_hA_d_RB = (h_A(R_A, R_B + ε, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter, c_A) - h_A(R_A, R_B - ε, p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter, c_A)) / (2ε)
    d_hB_d_RA = (h_B(R_A + ε, R_B, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra, c_B) - h_B(R_A - ε, R_B, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra, c_B)) / (2ε)
    d_hB_d_RB = (h_B(R_A, R_B + ε, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra, c_B) - h_B(R_A, R_B - ε, p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra, c_B)) / (2ε)

    # Form the Jacobian matrix.
    J = [d_hA_d_RA d_hA_d_RB;
        d_hB_d_RA d_hB_d_RB]

    # Compute and return the eigenvalues (which may be complex).
    return eigen(J).values
end


function time_series(p)
    R_A = zeros(T)
    R_B = zeros(T)
    R_A[1] = 0.95  # initial guess for R_A
    R_B[1] = 0.95
    for t in 2:T
        # Compute p_A^t and p_B^t based on previous values:
        p_A = g_A(R_A[t-1], R_B[t-1], p, c_A_plus_intra, c_A_minus_intra, c_B_plus_inter, c_B_minus_inter)
        p_B = g_B(R_A[t-1], R_B[t-1], p, c_A_plus_inter, c_A_minus_inter, c_B_plus_intra, c_B_minus_intra)
    
        # Solve for the next R_A and R_B using the implicit equations:
        R_A[t] = f_A(p_A, c_A)
        R_B[t] = f_B(p_B, c_B)
    end
    return R_A, R_B
end




# Plot the time series

# scatter(10:T, R_A[10:end], xlabel=L"t", ylabel=L"R", markersize=1, lw=2, label=L"R_A", markerstrokewidth=0, right_margin=10Plots.mm)
# scatter!(10:T, R_B[10:end], xlabel=L"t", ylabel=L"R", markersize=1, lw=2, label=L"R_B", markerstrokewidth=0)
# savefig("NS_timeseries.pdf")
# savefig("NS_timeseries.svg")









# For numerical calculation of Lyapunov exponent.

# R_A_list, R_B_list, p_list, real_part, imag_part = orbit_diagram()
# scatter(p_list, R_A_list, label=L"R_A", markersize=2, xlabel=L"p", ylabel=L"R", markerstrokewidth=0)
# scatter!(p_list, R_B_list, label=L"R_B", markersize=2, xlabel=L"p", ylabel=L"R", markerstrokewidth=0, right_margin=5Plots.mm, bottom_margin=2Plots.mm)
# plot!([p, p], [0, 1], lw=3, ls=:dash, color=:black, label="")

# writedlm("Lyapunov/NS_p_list.txt", p_list)
# writedlm("Lyapunov/NS_R_A_list.txt", R_A_list)
# writedlm("Lyapunov/NS_R_B_list.txt", R_B_list)
# savefig("orbit_diagram_NS.pdf")
# savefig("orbit_diagram_NS.svg")
# savefig("orbit_diagram_NS.png")







# Orbit diagram

# R_A_list, R_B_list, p_list, real_part, imag_part = orbit_diagram()
# scatter(p_list, R_A_list, label=L"R_A", markersize=2, xlabel=L"p", ylabel=L"R", markerstrokewidth=0)
# scatter!(p_list, R_B_list, label=L"R_B", markersize=2, xlabel=L"p", markerstrokewidth=0, right_margin=3Plots.mm)
# # savefig("NS_theory.pdf")
# # savefig("NS_theory.png")









