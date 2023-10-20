using Graphs
using ProgressMeter
using FastPercolation

# include(srcdir("importance_problems/importance_problems.jl"))
# include(srcdir("contagion_models/percolation_monte_carlo.jl"))
# include(srcdir("utils/read_graph.jl"))


function compute_micro(importance_problem, g, test_sets, observable_constructor; kwargs...)

    if importance_problem == "vaccination"
        micros = []
        for s in test_sets
            g_ = deepcopy(g)
            (rem_edge!(g_,u,v) for u in s for v in neighbors(g_,u))
            M = observable_constructor(s,g_)
            micro = percolation_MC(g_,M; verbose=true, kwargs...)
            push!(micros,micro)
        end
        return micros
    else
        M = observable_constructor(test_sets,g)
        micro = percolation_MC(g,M; verbose=true, kwargs...)
        if isa(M,SentinelObservable)
            micro = reduce_sentinel_Qs(micro)
        elseif  (isa(M,UtilityTransformationObservable) && isa(M.M,SentinelObservable))
            micro = map(reduce_sentinel_Qs,micro)
        end
        return micro
    end
end

function compute_qualities(importance_problem, micro,p,u_idx)
    if importance_problem == "vaccination"
        return [micro_to_canonical(p,micro_[u_idx]) for micro_ in micro]
    else
        return micro_to_canonical(p,micro[u_idx])
    end
end


function greed_percolation_all_p(
    g::Graph, k::Int64,
    observable_constructor::Function, 
    Us::Vector{Function},  ps::Vector{Float64}; 
    iteration_callback::Union{Function,Nothing}= nothing, 
    load_optimal_sets::Union{Function,Nothing}= nothing,
    importance_problem = nothing,
    kwargs...
)

    nodes = collect(1:nv(g))
    optimal_sets = [[Vector{Int64}() for u in Us] for p in ps]
    
    for k_ in 1:k

        x = load_optimal_sets(k_)
        if x !== nothing
            println("loading sets for $k_")
            optimal_sets = x
            continue
        end

        ## get all test sets for all ps
        test_sets = [[Vector{Vector{Int64}}() for u in Us] for p in ps]
        for p_idx in eachindex(ps)
            for u_idx in eachindex(Us)
                s = optimal_sets[p_idx][u_idx]
                test_sets[p_idx][u_idx] = [union(s,j) for j in filter(x -> x ∉ s, nodes)]
            end
        end
        
        ## get all unique sets
        unique_test_sets::Vector{Vector{Int64}} = []
        for p_idx in eachindex(ps)
            for u_idx in eachindex(Us)
                ts = test_sets[p_idx][u_idx]
                for s in ts
                    s_ = sort(s)
                    if s_ ∉ unique_test_sets
                        push!(unique_test_sets,s_)
                    end
                end
            end
        end
        println("\nnumber unique sets: $(length(unique_test_sets))")

        ## calculate quality for all unique sets and for all objectives
        micro = compute_micro(importance_problem,g,unique_test_sets,observable_constructor)

        for u_idx in eachindex(Us)
            for (p_idx,p) in enumerate(ps)

                # get qualities of all sets
                qualities = compute_qualities(importance_problem,micro, p, u_idx)
                uniqueset2quality = Dict(unique_test_sets .=> qualities)

                # get indices of all sets
                ts = test_sets[p_idx][u_idx]
                best_set = argmax( t -> uniqueset2quality[sort(t)], ts)
                optimal_sets[p_idx][u_idx] = best_set
            end
        end
        if iteration_callback !== nothing
            iteration_callback(optimal_sets,micro,k_)
        end
    end
    return optimal_sets

end

