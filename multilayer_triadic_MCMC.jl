# Monte Carlo simulation of triadic percolation with intra- and interlayer regulations

using Plots
using Random
using DelimitedFiles
using Base.Threads
using LaTeXStrings

function poisson_multilayer_network(N, c_A,  c_B, cp_AA, cn_AA, cp_BB, cn_BB, cp_AB, cn_AB, cp_BA, cn_BA)

    # Structural Poisson network with size N and average degree c.
    # Positive and negative regulatory networks are Poisson networks with average degree cp and cn.

    dict_nodes_stru_A = Dict{Int64,Vector{Int64}}()
    dict_edges_stru_A = Dict{Int64,Vector{Int64}}()
    dict_nodes_stru_B = Dict{Int64,Vector{Int64}}()
    dict_edges_stru_B = Dict{Int64,Vector{Int64}}()

    dict_edges_pos_AA = Dict{Int64,Vector{Int64}}()
    dict_edges_neg_AA = Dict{Int64,Vector{Int64}}()
    dict_edges_pos_BB = Dict{Int64,Vector{Int64}}()
    dict_edges_neg_BB = Dict{Int64,Vector{Int64}}()
    dict_edges_pos_AB = Dict{Int64,Vector{Int64}}()   # Links of layer A that are regulated by nodes in layer B.
    dict_edges_neg_AB = Dict{Int64,Vector{Int64}}()   
    dict_edges_pos_BA = Dict{Int64,Vector{Int64}}()   # Links of layer B that are regulated by nodes in layer A.
    dict_edges_neg_BA = Dict{Int64,Vector{Int64}}()


    upper_bound_A = N * c_A * 2  # preallocation, max number of edges
    upper_bound_B = N * c_B * 2  # preallocation, max number of edges

    for i = 1:N
        dict_nodes_stru_A[i] = Int64[]
        dict_nodes_stru_B[i] = Int64[]
    end

    for j = N+1:N+upper_bound_A
        dict_edges_stru_A[j] = Int64[]
    end

    for j = N+1:N+upper_bound_B
        dict_edges_stru_B[j] = Int64[]
    end
    edgeID_A = N + 1
    edgeID_B = N + 1


    # Fill layer A

    for i = 1:N
        for j = i+1:N
            if rand() < c_A / N
                push!(dict_nodes_stru_A[i], edgeID_A)
                push!(dict_nodes_stru_A[j], edgeID_A)
                push!(dict_edges_stru_A[edgeID_A], i)
                push!(dict_edges_stru_A[edgeID_A], j)
                edgeID_A += 1
            end
        end
    end

    max_edgeID_A = edgeID_A - 1
    for i = max_edgeID_A+1:upper_bound_A+N
        delete!(dict_edges_stru_A, i)
    end

    M_A = max_edgeID_A - N  # Num of edges

    ###########################################################
    # Fill layer B
    for i = 1:N
        for j = i+1:N
            if rand() < c_B / N
                push!(dict_nodes_stru_B[i], edgeID_B)
                push!(dict_nodes_stru_B[j], edgeID_B)
                push!(dict_edges_stru_B[edgeID_B], i)
                push!(dict_edges_stru_B[edgeID_B], j)
                edgeID_B += 1
            end
        end
    end

    max_edgeID_B = edgeID_B - 1
    for i = max_edgeID_B+1:upper_bound_B+N
        delete!(dict_edges_stru_B, i)
    end

    M_B = max_edgeID_B - N  # Num of edges
    ###########################################################

    # Here `i` is the edgeID, hence

    # `i` ranges from N+1 to N+M_A for layer A and N+1 to N+M_B for layer B

    # Fill intra-layer positive and negative regulations for layer A

    for i = N+1:N+M_A
        dict_edges_pos_AA[i] = Int64[]
        dict_edges_neg_AA[i] = Int64[]
    end

    for edge = N+1:N+M_A
        for node = 1:N
            rand_p = rand()
            if !(node in dict_edges_stru_A[edge])  # This is only necessary for intra-layer regulations, even if it is not necessary on large network limit.
                if rand_p < cp_AA / N
                    push!(dict_edges_pos_AA[edge], node)
                elseif rand_p > cp_AA / N && rand_p < (cn_AA + cp_AA) / N
                    push!(dict_edges_neg_AA[edge], node)
                end
            end
        end
    end


    # Fill intra-layer positive and negative regulations for layer B
    for i = N+1:N+M_B
        dict_edges_pos_BB[i] = Int64[]
        dict_edges_neg_BB[i] = Int64[]
    end

    for edge = N+1:N+M_B
        for node = 1:N
            rand_p = rand()
            if !(node in dict_edges_stru_B[edge])
                if rand_p < cp_BB / N
                    push!(dict_edges_pos_BB[edge], node)
                elseif rand_p > cp_BB / N && rand_p < (cn_BB + cp_BB) / N
                    push!(dict_edges_neg_BB[edge], node)
                end
            end
        end
    end

    # Fill interlayer positive and negative regulations for layer A
    for i = N+1:N+M_A
        dict_edges_pos_AB[i] = Int64[]
        dict_edges_neg_AB[i] = Int64[]
    end
    for edge = N+1:N+M_A
        for node = 1:N
            rand_p = rand()
            if rand_p < cp_AB / N
                push!(dict_edges_pos_AB[edge], node)
            elseif rand_p > cp_AB / N && rand_p < (cn_AB + cp_AB) / N
                push!(dict_edges_neg_AB[edge], node)
            end
        end
    end

    # Fill interlayer positive and negative regulations for layer B
    for i = N+1:N+M_B
        dict_edges_pos_BA[i] = Int64[]
        dict_edges_neg_BA[i] = Int64[]
    end
    for edge = N+1:N+M_B
        for node = 1:N
            rand_p = rand()
            if rand_p < cp_BA / N
                push!(dict_edges_pos_BA[edge], node)
            elseif rand_p > cp_BA / N && rand_p < (cn_BA + cp_BA) / N
                push!(dict_edges_neg_BA[edge], node)
            end
        end
    end

    # Return the network structure
    return dict_nodes_stru_A, dict_edges_stru_A, dict_nodes_stru_B, dict_edges_stru_B,
           dict_edges_pos_AA, dict_edges_neg_AA,
           dict_edges_pos_BB, dict_edges_neg_BB,
           dict_edges_pos_AB, dict_edges_neg_AB,
           dict_edges_pos_BA, dict_edges_neg_BA
end


function bfs(node, dict_nodes, dict_edges, counted, retained, clusterID, ID_list)

    # Breadth first search algorithm
    # Given the network structure, return the size of component that `node` belongs to
    # and assign this component an ID called `clusterID`

    # counted: 0 not counted, 1 counted
    # retained: 0 not retained, 1 retained

    queue_edges = Int64[]
    queue_nodes = Int64[]

    push!(queue_nodes, node)

    counted[node] = 1
    ID_list[node] = clusterID
    cluster_size = 1

    while !isempty(queue_nodes)
        node = popfirst!(queue_nodes)

        for edge in dict_nodes[node]
            if (counted[edge] == 0) & (retained[edge] == 1)
                push!(queue_edges, edge)
                counted[edge] = 1
                ID_list[edge] = clusterID
            end
        end

        while !isempty(queue_edges)
            edge = popfirst!(queue_edges)
            for node in dict_edges[edge]
                if (counted[node] == 0) & (retained[node] == 1)
                    cluster_size += 1
                    push!(queue_nodes, node)
                    counted[node] = 1
                    ID_list[node] = clusterID
                end
            end
        end
    end
    return cluster_size, counted, ID_list, clusterID
end


function gc(dict_nodes, dict_edges, retained)

    # Using the BFS algorithm, given the network structure, return the size of the largest component and its ID

    num_nodes = length(dict_nodes)
    num_edges = length(dict_edges)
    counted = zeros(Int64, num_nodes + num_edges)
    ID_list = zeros(Int64, num_nodes + num_edges)
    max_cluster_size = 0
    max_clusterID = 0
    clusterID = 0
    for i = 1:num_nodes
        if (counted[i] == 0) & (retained[i] == 1)
            clusterID += 1
            cluster_size, counted, ID_list, clusterID = bfs(i, dict_nodes, dict_edges, counted, retained, clusterID, ID_list)
            if cluster_size > max_cluster_size
                max_cluster_size = cluster_size
                max_clusterID = clusterID
            end
        else
            continue
        end
    end
    return max_cluster_size, max_clusterID, ID_list
end


function simulation(dict_nodes_stru_A, dict_edges_stru_A, dict_nodes_stru_B, dict_edges_stru_B,
           dict_edges_pos_AA, dict_edges_neg_AA,
           dict_edges_pos_BB, dict_edges_neg_BB,
           dict_edges_pos_AB, dict_edges_neg_AB,
           dict_edges_pos_BA, dict_edges_neg_BA, p, retained_A, retained_B, Tmax; 
           flag_AA=false, flag_BB=false, flag_AB=false, flag_BA=false)


    # When flag_XX is true, cp_XX is infinity and the input cp_XX is not used. By default, all flags are false.
    # flags are used as keywards arguents.

    # TODO: - deal with infinity cp_pos cases: in both function generating the network, as well as simulation/


    #  Numerical simulation on the given network
    #  Retained: initial condition of edges. 1: retained, 0: damaged
    #  Tmax: duration of the time series

    N = length(dict_nodes_stru_A)

    active_node_A = ones(N)
    active_node_B = ones(N)

    # Finding the giant component and deactivate the nodes and links which are not in the giant component.

    # Collect the time series
    cluster_size_list_A = zeros(Tmax)
    cluster_size_list_B = zeros(Tmax)
    T_list = zeros(Tmax)

    for t = 1:Tmax

        # Find the giant component
        max_cluster_size_A, max_clusterID_A, ID_list_A = gc(dict_nodes_stru_A, dict_edges_stru_A, retained_A)
        max_cluster_size_B, max_clusterID_B, ID_list_B = gc(dict_nodes_stru_B, dict_edges_stru_B, retained_B)

        for i = 1:N
            if ID_list_A[i] != max_clusterID_A
                active_node_A[i] = 0      # Only nodes in the giant component are active.
            else
                active_node_A[i] = 1
            end
        end

        for i = 1:N
            if ID_list_B[i] != max_clusterID_B
                active_node_B[i] = 0      # Only nodes in the giant component are active.
            else
                active_node_B[i] = 1
            end
        end

        cluster_size_list_A[t] = (max_cluster_size_A / N)
        cluster_size_list_B[t] = (max_cluster_size_B / N)

        T_list[t] = t

        # Update the pL_A
        for edge in keys(dict_edges_stru_A)

            pos_regulation_intra = dict_edges_pos_AA[edge]
            neg_regulation_intra = dict_edges_neg_AA[edge]
            pos_regulation_inter = dict_edges_pos_AB[edge]
            neg_regulation_inter = dict_edges_neg_AB[edge]

            sum_pos_inter = 0
            sum_neg_inter = 0
            sum_pos_intra = 0
            sum_neg_intra = 0

            for pos_node_intra in pos_regulation_intra
                sum_pos_intra += active_node_A[pos_node_intra]
            end

            for neg_node_intra in neg_regulation_intra
                sum_neg_intra += active_node_A[neg_node_intra]
            end

            for pos_node_inter in pos_regulation_inter
                sum_pos_inter += active_node_B[pos_node_inter]
            end

            for neg_node_inter in neg_regulation_inter
                sum_neg_inter += active_node_B[neg_node_inter]
            end

            at_least_one_pos_each_layer = (sum_pos_intra > 0) && (sum_pos_inter > 0)
            no_negative_in_neither_layer = (sum_neg_intra == 0) && (sum_neg_inter == 0)

            if flag_AA || flag_AB   # If intra-layer or inter-layer positive regulation is infinite, then edges are retained as long as there is no negative regulation (and random damage)
                if (no_negative_in_neither_layer) && (rand() < p)
                    retained_A[edge] = 1
                else
                    retained_A[edge] = 0
                end
            else
                if (at_least_one_pos_each_layer) && (no_negative_in_neither_layer) && (rand() < p)
                    retained_A[edge] = 1
                else
                    retained_A[edge] = 0
                end
            end
        end


        # Update the pL_B
        for edge in keys(dict_edges_stru_B)
            pos_regulation_intra = dict_edges_pos_BB[edge]
            neg_regulation_intra = dict_edges_neg_BB[edge]
            pos_regulation_inter = dict_edges_pos_BA[edge]
            neg_regulation_inter = dict_edges_neg_BA[edge]

            sum_pos_intra = 0
            sum_neg_intra = 0
            sum_pos_inter = 0
            sum_neg_inter = 0

            for pos_node_intra in pos_regulation_intra
                sum_pos_intra += active_node_B[pos_node_intra]
            end

            for neg_node_intra in neg_regulation_intra
                sum_neg_intra += active_node_B[neg_node_intra]
            end

            for pos_node_inter in pos_regulation_inter
                sum_pos_inter += active_node_A[pos_node_inter]
            end

            for neg_node_inter in neg_regulation_inter
                sum_neg_inter += active_node_A[neg_node_inter]
            end

            at_least_one_pos_each_layer = (sum_pos_intra > 0) && (sum_pos_inter > 0)
            no_negative_in_neither_layer = (sum_neg_intra == 0) && (sum_neg_inter == 0)

            if flag_BB || flag_BA   # If intra-layer or inter-layer positive regulation is infinite, then edges are retained as long as there is no negative regulation (and random damage)
                if (no_negative_in_neither_layer) && (rand() < p)
                    retained_B[edge] = 1
                else
                    retained_B[edge] = 0
                end
            else
                if (at_least_one_pos_each_layer) && (no_negative_in_neither_layer) && (rand() < p)
                    retained_B[edge] = 1
                else
                    retained_B[edge] = 0
                end
            end
        end
        
    end
    return cluster_size_list_A, cluster_size_list_B, T_list, retained_A, retained_B
end

function simulation_orbit()
    N = 10000
    c_A = 30
    c_B = 30

    cp_BB = 0 #infty
    cn_BB = 0
    cp_AB = 3
    cn_AB = 0

    cp_AA = 5
    cn_AA = 1.5
    cp_BA = 0 #infty
    cn_BA = 3

    dict_nodes_stru_A, dict_edges_stru_A, 
    dict_nodes_stru_B, dict_edges_stru_B,
    dict_edges_pos_AA, dict_edges_neg_AA,
    dict_edges_pos_BB, dict_edges_neg_BB,
    dict_edges_pos_AB, dict_edges_neg_AB,
    dict_edges_pos_BA, dict_edges_neg_BA = poisson_multilayer_network(N, c_A, c_B, cp_AA, cn_AA, cp_BB, cn_BB, cp_AB, cn_AB, cp_BA, cn_BA)

    R_collection_A = Float64[]
    R_collection_B = Float64[]
    p_collection = Float64[]
    n = 50
    Tmax = 100

    # Create a range of p values
    p_values = collect(0.8:-0.005:0.65)
    num_p = length(p_values)

    # Preallocate result arrays
    R_collection_A = zeros(n * num_p)
    R_collection_B = zeros(n * num_p)
    p_collection = zeros(n * num_p)
    pA_0 = 0.04
    pB_0 = 0.04
   
    @threads for i in 1:num_p
        p = p_values[i]
        @show p
        
        # retained_A = [ones(length(dict_nodes_stru_A)); rand(length(dict_edges_stru_A)) .< pA_0]
        # retained_B = [ones(length(dict_nodes_stru_B)); rand(length(dict_edges_stru_B)) .< pB_0]


        retained_A = [ones(length(dict_nodes_stru_A)); rand(length(dict_edges_stru_A)) .< pA_0]
        retained_B = [ones(length(dict_nodes_stru_B)); rand(length(dict_edges_stru_B)) .< pB_0]


        R_A, R_B, T_list, retained_A, retained_B = simulation(
            dict_nodes_stru_A, dict_edges_stru_A,
            dict_nodes_stru_B, dict_edges_stru_B,
            dict_edges_pos_AA, dict_edges_neg_AA,
            dict_edges_pos_BB, dict_edges_neg_BB,
            dict_edges_pos_AB, dict_edges_neg_AB,
            dict_edges_pos_BA, dict_edges_neg_BA,
            p, retained_A, retained_B, Tmax, flag_BA=true, flag_BB=true
        )

        idx_start = (i - 1) * n + 1
        idx_end = i * n
        R_collection_A[idx_start:idx_end] .= R_A[end-n+1:end]
        R_collection_B[idx_start:idx_end] .= R_B[end-n+1:end]
        p_collection[idx_start:idx_end] .= p
    end

    # for p = 1:-0.01:0
    #     @show p
    #     pA_0 = 0.1
    #     pB_0 = 0.1
    #     retained_A = [ones(length(dict_nodes_stru_A)); rand(length(dict_edges_stru_A)) .< pA_0]
    #     retained_B = [ones(length(dict_nodes_stru_B)); rand(length(dict_edges_stru_B)) .< pB_0]
    #     R_A, R_B, T_list, retained_A, retained_B = simulation(dict_nodes_stru_A, dict_edges_stru_A, dict_nodes_stru_B, dict_edges_stru_B,
    #         dict_edges_pos_AA, dict_edges_neg_AA,
    #         dict_edges_pos_BB, dict_edges_neg_BB,
    #         dict_edges_pos_AB, dict_edges_neg_AB,
    #         dict_edges_pos_BA, dict_edges_neg_BA, p, retained_A, retained_B, Tmax,
    #         flag_BB=true, flag_BA=true)  # Intra-layer positive regulation is infinite, inter-layer positive regulation is infinite.
    #     R_collection_A = [R_collection_A; R_A[end-n+1:end]]
    #     R_collection_B = [R_collection_B; R_B[end-n+1:end]]
    #     p_collection = [p_collection; p * ones(n)]
    # end
    return p_collection, R_collection_A, R_collection_B
end

function generate_network()
    N = 10000
    c_A = 30
    c_B = 30

    cp_BB = 0 #infty
    cn_BB = 0
    cp_AB = 3
    cn_AB = 0

    cp_AA = 5
    cn_AA = 1.5
    cp_BA = 0 #inf
    cn_BA = 3

    dict_nodes_stru_A, dict_edges_stru_A,
    dict_nodes_stru_B, dict_edges_stru_B,
    dict_edges_pos_AA, dict_edges_neg_AA,
    dict_edges_pos_BB, dict_edges_neg_BB,
    dict_edges_pos_AB, dict_edges_neg_AB,
    dict_edges_pos_BA, dict_edges_neg_BA = poisson_multilayer_network(N, c_A, c_B, cp_AA, cn_AA, cp_BB, cn_BB, cp_AB, cn_AB, cp_BA, cn_BA)
    
    return dict_nodes_stru_A, dict_edges_stru_A,
           dict_nodes_stru_B, dict_edges_stru_B,
           dict_edges_pos_AA, dict_edges_neg_AA,
           dict_edges_pos_BB, dict_edges_neg_BB,
           dict_edges_pos_AB, dict_edges_neg_AB,
           dict_edges_pos_BA, dict_edges_neg_BA
end

function time_series(p, dict_nodes_stru_A, dict_edges_stru_A,
    dict_nodes_stru_B, dict_edges_stru_B,
    dict_edges_pos_AA, dict_edges_neg_AA,
    dict_edges_pos_BB, dict_edges_neg_BB,
    dict_edges_pos_AB, dict_edges_neg_AB,
    dict_edges_pos_BA, dict_edges_neg_BA)
   
    Tmax = 500
    pA_0 = 0.04
    pB_0 = 0.04

    retained_A = [ones(length(dict_nodes_stru_A)); rand(length(dict_edges_stru_A)) .< pA_0]
    retained_B = [ones(length(dict_nodes_stru_B)); rand(length(dict_edges_stru_B)) .< pB_0]

    R_A, R_B, T_list, retained_A, retained_B = simulation(
        dict_nodes_stru_A, dict_edges_stru_A,
        dict_nodes_stru_B, dict_edges_stru_B,
        dict_edges_pos_AA, dict_edges_neg_AA,
        dict_edges_pos_BB, dict_edges_neg_BB,
        dict_edges_pos_AB, dict_edges_neg_AB,
        dict_edges_pos_BA, dict_edges_neg_BA,
        p, retained_A, retained_B, Tmax, flag_BA=true, flag_BB=true
    )
    return R_A, R_B, T_list
end

dict_nodes_stru_A, dict_edges_stru_A,
dict_nodes_stru_B, dict_edges_stru_B,
dict_edges_pos_AA, dict_edges_neg_AA,
dict_edges_pos_BB, dict_edges_neg_BB,
dict_edges_pos_AB, dict_edges_neg_AB,
dict_edges_pos_BA, dict_edges_neg_BA = generate_network()

p = 0.72

@time R_A, R_B, T_list = time_series(p, dict_nodes_stru_A, dict_edges_stru_A,
    dict_nodes_stru_B, dict_edges_stru_B,
    dict_edges_pos_AA, dict_edges_neg_AA,
    dict_edges_pos_BB, dict_edges_neg_BB,
    dict_edges_pos_AB, dict_edges_neg_AB,
    dict_edges_pos_BA, dict_edges_neg_BA)

# plot(R_A[end-100:end], R_B[end-100:end], xlabel=L"R_A", ylabel=L"R_B",
#     markersize=1, lw=2, label="", markerstrokewidth=0, right_margin=10Plots.mm, framestyle=:box,
#     grid=:off,
#     tickfontsize=16, color=:red)
plot(T_list[end-100:end], R_A[end-100:end], xlabel=L"T", ylabel=L"R_A",
    markersize=5, lw=2, label="", markerstrokewidth=0, right_margin=10Plots.mm, framestyle=:box,
    grid=:off, marker=:circle,
    tickfontsize=16, color=:red)

plot!(T_list[end-100:end], R_B[end-100:end], xlabel=L"T", ylabel=L"R_A",
    markersize=5, lw=2, label="", markerstrokewidth=0, right_margin=10Plots.mm, framestyle=:box,
    grid=:off, marker=:circle,
    tickfontsize=16, color=:black)


# writedlm("R_A_simulation.txt", R_A)
# writedlm("R_B_simulation.txt", R_B)

# savefig("time_series_NS_MC_500k.png")


# @time p_collection, R_collection_A, R_collection_B = simulation_orbit()
# scatter(p_collection, R_collection_A, label=L"R_A", markersize=2, xlabel=L"p", ylabel=L"R", markerstrokewidth=0, color=:red)
# scatter!(p_collection, R_collection_B, label=L"R_B", markersize=2, xlabel=L"p", ylabel=L"R", markerstrokewidth=0, right_margin=3Plots.mm, color=:black, tickfontsize=16)

R_A = readdlm("R_A_simulation.txt")
R_B = readdlm("R_B_simulation.txt")
scale = 401:501
T_list = collect(1:length(R_A))
plot(T_list[scale], R_A[scale], marker=:circle, xlabel=L"t", ylabel=L"R",
     lw=2, markerstrokewidth=0, right_margin=10Plots.mm, framestyle=:box,
    grid=:off, label=L"R_A", color=:red,
    tickfontsize=18)
plot!(T_list[scale], R_B[scale], marker=:circle, xlabel=L"t", ylabel=L"R",
    lw=2, markerstrokewidth=0, right_margin=10Plots.mm, framestyle=:box,
    grid=:off, label=L"R_B",legendfont= 14,legend=:bottomright, color=:black,
    tickfontsize=18)
ylims!(-0.04, 1.02)
savefig("neimark_sacker_simulation_NEW.pdf")
savefig("neimark_sacker_simulation_NEW.svg")
savefig("neimark_sacker_simulation_NEW.png")