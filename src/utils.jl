using DrWatson
@quickactivate "uncertain-influence-max"

using DelimitedFiles
using Graphs
using JSON

function load_edgelist(network_name)

    file_path = datadir("networks/$network_name/edgelist.txt")
    edgelist = readdlm(file_path, ' ', Int)
    edgelist .+= 1
    return [row[:] for row in eachrow(edgelist)]
end


function read_graph(network_name)
    
    adj_list = read_adjlist(network_name)
    
    g = SimpleGraph(length(adj_list),0)
    for i in eachindex(adj_list)
        for j in adj_list[i]
            add_edge!(g,i,j)
        end
    end
    return g
end


function read_adjlist(network_name)
    
    path = datadir("networks/$network_name/adjlist.txt")
    open(path,"r") do f
        adj_list = JSON.parse(read(f, String))
    end
    adj_list = convert(Vector{Vector{Int64}},adj_list)
    adj_list = map(x -> x .+ 1,adj_list) # julia is one-indexed
    adj_list
end

function read_threshold(network_name)
    path = datadir("networks/$network_name/threshold.dat")
    s = read(open(path,"r"),String)
    return parse(Float64,s)
end


NETWORKS = [
    "facebook_friends",
    "product_space/SITC",
    "faa_routes",
    "celegans_2019/male_gap_junction",
    "dolphins",
    "train_terrorists",
    "urban_streets/vienna",
    "student_cooperation",
    "polbooks",
    "malaria_genes/HVR_1",
    "lesmis",
    "karate/78",
    "jazz_collab",
    "euroroad",
    "crime"
]

IPs = [
    "influence_maximization",
    "vaccination",
    "sentinel_surveillance"
]