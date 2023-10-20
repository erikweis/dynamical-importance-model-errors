using DrWatson
@quickactivate "node-importance-model-errors"

include(srcdir("greedy.jl"))
include(srcdir("utils.jl"))

using FastPercolation
using Graphs

function get_Us(u_exps,thresholds,importance_problem)

    Us::Vector{Function} = Vector{Function}()
    if importance_problem == "influence_maximization"
        Us = [ x->x^α for α in u_exps]
        append!(Us,[x->Int(x<t) for t in thresholds])
    elseif importance_problem == "vaccination"
        Us = [ x->-abs(x)^α for α in u_exps]
        append!(Us,[x->Int(x>t) for t in thresholds])
    else
        Us = [ x->x^α for α in u_exps]
        append!(Us,[x->Int(x<t) for t in thresholds])
    end
    return Us
end


function f(params)

    @unpack network_name, p_rels, u_exps, thresholds, importance_problem, num_samples = params

    p_c = read_threshold(network_name)
    ps = p_rels .* p_c
    ps = [p for p in ps if 0<=p<=1]

    g = read_graph(network_name)
    k = Int(round(sqrt(nv(g))))

    Us = get_Us(u_exps,thresholds,importance_problem)

    ip2observable = Dict(
        "influence_maximization"=>InfluenceMaxObservable,
        "vaccination"=>ImmunizationObservable,
        "sentinel_surveillance"=>SentinelObservable
    )
    observable_constructor = (ts,g) -> UtilityTransformationObservable(
        ip2observable[importance_problem](ts,nv(g),ne(g)),
        Us
    )

    function callback(optimal_sets,micro,k)

        result = copy(params)
        result["k"] = k
        path = datadir("optimal_sets",network_name,savename(result,"jld2"))

        result["optimal_sets"] = optimal_sets
        #result["micro"] = micro

        wsave(path,result)
    end

    function load_optimal_sets(k)
        d = deepcopy(params)
        d["k"] = k
        path = datadir("optimal_sets",network_name,savename(d,"jld2"))
        if isfile(path)
            return wload(path)["optimal_sets"]
        end
        nothing
    end

    greed_percolation_all_p(g, k, observable_constructor, Us, ps; iteration_callback=callback,num_samples = num_samples, importance_problem = importance_problem, load_optimal_sets = load_optimal_sets)
end

ip2num_samples = Dict(
    "influence_maximization"=>10^5,
    "sentinel_surveillance"=>10^3,
    "vaccination"=>10^3
)

function get_param_dicts()

    d = Dict{String,Any}(
        "network_name"=>NETWORKS,#["karate_78"],
        "p_rels"=>[collect(0.01:0.1:5)],
        "u_exps"=>[2 .^ collect(-2:0.5:2)],
        "thresholds"=>[[0.05,0.1,0.2,0.3]],
        "importance_problem"=>IPs
    )
    ds::Vector{Dict{String,Any}} = dict_list(d;)
    for d in ds
        d["num_samples"] = ip2num_samples[d["importance_problem"]]
        if d["importance_problem"] == "sentinel_surveillance"
            d["thresholds"] = [2,3,4]
        end
    end
    ds
end

function param_sweep(idx_min, idx_max)

    dicts = get_param_dicts()[idx_min:idx_max]
    for d in dicts
        f(d)
    end
end

println(length(get_param_dicts()))
param_sweep(1,3)

const idx_min = parse(Int,ARGS[1])
const idx_max = parse(Int,ARGS[2])

println("starting to run params from index $idx_min to $idx_max")
param_sweep( idx_min, idx_max)
    




