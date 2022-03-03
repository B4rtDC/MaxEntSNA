module Analyse
    using CSV, DataFrames, JSON3 # Data import
    #using Distributions
    #using 
    using Graphs, MetaGraphs
    #using Plots
    #using Plots.PlotMeasures
    using PyCall # to load up python results
    pickle = pyimport("pickle")
    using SparseArrays
    

    import Dates: now, format
    #import StatsBase: countmap


    #IG = pyimport("igraph")
	#LA = pyimport("leidenalg")

    
    """
        DTG()

    return DTG string for logging
    """
    function DTG()
        return "$(format(now(), "YYYY-mm-dd HH:MM:SS"))"
    end

    """
        loaddata(userfile::Union{String,Array{String,1}}, msgfiles::Union{String,Array{String,1}})

    From a (list of) user and message CSV file(s), get the dataframes. Works on the raw files downloaded from the Twitter misinformation report.

    ### kwargs
    * :usrfields : `Array{Symbol,1}` - the columns of the user dataframe you want to use (default: all)
    * :msgfields : `Array{Symbol,1}` - the columns of the message dataframe you want to use (default: all)

    ### Examples
    ```Julia
    loaddata("/path/to/files/userfile.csv", 
                ["/path/to/files/messagefile_1.csv", 
                    "/path/to/files/messagefile_1.csv"])

    ```

    ### Notes
    * Makes use of [`DataFrame!`](https://juliadata.github.io/DataFrames.jl/stable/lib/functions/#DataFrames.DataFrame!) 
    combined with [`CSV.File`](https://juliadata.github.io/CSV.jl/stable/#CSV.File) with a selection of columns for increased performance.
    Column types can also be passed for an additional speedup (see also [`CSV.File`](https://juliadata.github.io/CSV.jl/stable/#CSV.File) kwargs)
    * If multiple treads are available, they will be used to read the CSV files in parallel.

    See also: ['dropdataframes!'](@ref), ['setgraph!'](@ref), ['describe!'](@ref)
    """
    function loaddata(userfiles::Array{String,1}, msgfiles::Array{String,1}; kwargs...)
        @info "$(DTG()) - loading Twitter information operations archive files"
        # default settings
        defusrfilter = [:userid, :user_screen_name]
        defmsgfilter = [:tweetid, :userid, :is_retweet, :retweet_tweetid, :retweet_userid, :in_reply_to_tweetid, :in_reply_to_userid]

        # Load up CSV data with only selected columns
        usrdf = vcat([DataFrame(CSV.File(f; select=get(kwargs,:usrfields, defusrfilter), 
                                            limit=get(kwargs,:limit, nothing))) 
                        for f in userfiles]...)
        @info "$(DTG()) - finished with users" 
        msgdf = vcat([DataFrame(CSV.File(f; select=get(kwargs,:msgfields, defmsgfilter), 
                                            limit=get(kwargs,:limit, nothing),
                                            truestrings=["true", "True"], falsestrings=["false", "False"], # added for stablity
                                            types=Dict(:is_retweet => Bool)))
                        for f in msgfiles ]...)
        @info "$(DTG()) - finished with messages" 
        @info "$(DTG()) - finished desinfo loading data files"
        return usrdf, msgdf
    end

    # additional methods for flexibility
    loaddata(userfile::String, msfgile::String; kwargs...) = loaddata([userfile], [msfgile]; kwargs...)
    loaddata(userfile::String, msfgiles::Array{String,1}; kwargs...) = loaddata([userfile], msfgiles; kwargs...)
    loaddata(userfiles::Array{String,1}, msfgile::String; kwargs...) = loaddata(userfiles, [msfgile]; kwargs...)

    """
        loaddata(path::String; kwargs...)

    Load data files from a path containing user and message CSV files.

    ### kwargs:
    * `:usrdata` : `String` - pattern used to recognize user data file(s) (default: "users_csv_hashed")
    * `:msgdata` : `String` - pattern used to recognize message data file(s) (default: "tweets_csv_hashed")
    * `:usrfields` : `Array{Symbol,1}` - the fields of the user dataframe you want to use (default: all)
    * `:msgfields` : `Array{Symbol,1}` - the fields of the message dataframe you want to use (default: all)
    * `:usrfilter` : `Array{String,1}` - an additional pattern that should found in the message data files (e.g. "_2020", default:[""])
    * `:msgfilter` : `Array{String,1}` - an additional pattern that should found in the message data files (e.g. "_2020", default:[""])

    See also: ['dropdataframes!'](@ref), ['setgraph!'](@ref), ['describe!'](@ref)

    """
    function loaddata(path::String; kwargs...)
        usrfiles = joinpath.(path, filter(x -> occursin(get(kwargs, :usrdata, "users_csv_hashed"), x) && 
                                                    !occursin(".jsonl",x) &&
                                                any([occursin(t,x) for t in get(kwargs, :usrfilter, [""]) ]), 
                                readdir(path)))
        msgfiles = joinpath.(path, filter(x -> occursin(get(kwargs, :usrdata, "tweets_csv_hashed"),x) && 
                                                    !occursin(".jsonl",x) &&
                                                any([occursin(t,x) for t in get(kwargs, :msgfilter, [""]) ]), 
                                readdir(path)))
        @info """$(DTG()) - files loaded from path <$(path)> :\nusers:\n\t$(join(usrfiles, "\n\t"))\nmessages:\n\t$(join(msgfiles, "\n\t"))""" 
        return loaddata(usrfiles, msgfiles; kwargs...)
    end

    """
        readjsonl(p::String; kind::Symbol=false)

    Read external tweets from jsonl file into dataframe holding user id, user description & post id

    *Note*: external tweets that do not contain a hashtag will not be included in this dataframe by default
    
    # Arguments: 
    - `p::String`: the file that contains the additional tweets that are required for the interaction graph
    - `keepall:::Bool`: keep all messages or only the ones with a hashtag

    """
    function readjsonlforhastags(p::String; keepall::Bool=false, kwargs...)
        @info "$(DTG()) - reading $(p)"
        userids = Vector{String}()         # user ids
        userdes = Vector{String}()         # user description
        tags    = Vector{Vector{String}}() # tags used
        for line in eachline(p)
            msg = JSON3.read(line)
            # work if extended tweet
            if haskey(msg, :extended_tweet)
                @warn "$(DTG()) - EXTENDED TWEET FOUND!"
                throw(ArgumentError("Parsing of extended_tweet not implemented yet..."))
                continue
            end
            # normal work
            if length(msg[:entities][:hashtags]) > 0
                push!(userids, msg[:user][:id_str])
                push!(userdes, msg[:user][:screen_name])
                push!(tags, [tag[:text] for tag in msg[:entities][:hashtags]])
            end
            # if we want to keep it all (users)
            if keepall
                push!(userids, msg[:user][:id_str])
                push!(userdes, msg[:user][:screen_name])
                push!(tags, String[])
            end
        end

        return DataFrame([  :userid => userids, 
                            :user_screen_name => userdes,
                            :tags => tags
                            ])
    end

    # regex constant for bipartitehashtaggraph
    const tagpat = r"""\'(\S+)\'"""

    """
        bipartitehashtaggraph(udf, mdf, externaldfs...;
                        dropisolatedusers::Bool=true,
                        kwargs...)

    Generate layer counts, node mapping and edge list for bipartite hashtag graph 

    ## Arguments:
    * `udf::DataFrame`: dataframe containing flagged user info
    * `mdf::DataFrame`: dataframe containing flagged user messages
    * `externaldfs::DataFrame`: additional dataframe(s) containing messages

    ## Kwargs:
    * `dropisolatedusers` : `Bool` - included all users in the final network, even the ones that are not connected (default:true)
    
    
    """
    function bipartitehashtaggraph(udf, mdf, externaldfs...;
                                    dropisolatedusers::Bool=true,
                                    kwargs...)

        users = Dict{String, Dict}()    # forward mapping userid::String => metadata
        tags  = Set{String}()			# set of tags
        edges = Set{NTuple{2,String}}()	# list of edges
        
        # 1. run over flagged users info
        for row in eachrow(unique(udf, :userid))
            users[row.userid] = Dict(:kind => :user,
                                    :userid => row.userid, 
                                    :Label => row.user_screen_name, 
                                    :flagged => true)
        end

        # 2. run over flagged users activity
        for row in eachrow(mdf)
            if !ismissing(row.hashtags) && length(row.hashtags) > 2
                for m in eachmatch(tagpat, row.hashtags)
                tag = m.captures[1]
                # add tag in list
                push!(tags, tag)
                # add edge
                push!(edges, (row.userid, tag))
                end
            end
        end

        # 3. run over external data
        for xdf in externaldfs
            for row in eachrow(xdf)
            # add users if required
            if !haskey(users, row.userid)
                users[row.userid] = Dict(:kind => :user,
                                    :userid => row.userid, 
                                    :Label => row.user_screen_name, 
                                    :flagged => false)
            end
            # add edges to set
            if length(row.tags) > 0
                for tag in row.tags
                    # add tag if required
                    push!(tags, tag)
                    # add edge
                    push!(edges, (row.userid, tag))
                    end
                end
            end
        end

        # 4. get node counts (by default only connected nodes)
        if dropisolatedusers
            actualusers = Set(x[1] for x in edges)
        else
            actualusers = keys(users)	
        end
        N_⊤ = length(actualusers)
        N_⊥ = length(tags)

        # 5. define actual nodes in the network 
        vertices = Dict{Int, Dict}()
        for user in actualusers
            vertices[length(vertices)+1] = users[user]
        end
        for tag in tags
            vertices[length(vertices)+1] = Dict(:kind => :hashtag,
                                                :Label => tag)
        end

        # 6. define reverse mapping
        invertices = Dict{String, Int}()
        for (key, val) in vertices # Int => Dict
            if isequal(val[:kind], :user)
                invertices[val[:userid]] = key
            end
            if isequal(val[:kind], :hashtag)
                invertices[val[:Label]] = key
            end
        end

        # 7. map edges to node numbers
        edges = [(invertices[x[1]], invertices[x[2]]) for x in edges]

        @info "$(DTG()) - built components of bipartite hashtag graph ($(N_⊤) users, $(N_⊥) hashtags)"
        return N_⊤, N_⊥, vertices, invertices, edges
    
    end

    function bipartiteadjacencymatrix(edges)
        I = Int[]
        J = Int[]
        V = Bool[]
        for e in edges
            push!(I, e[1])
            push!(J, e[2])
            push!(V, true)
            push!(I, e[2])
            push!(J, e[1])
            push!(V, true)
        end

        @info "$(DTG()) - built adjacency matrix"
        return I,J,V
    end

    # regex patterns for data extraction
    const r_m_id = r""""id_str": "(\d+)", "ful"""
	const r_user = r""""user": {"id": \d*, "id_str": "(\d*)", "name": ".*?", "screen_name": "(.*?)","""



    """
        write_edgelist(edges::Array{Tuple}, fname::String)

    write out edgeslist to a csv file for use with the Python NEMtropy scipt. 


    used for the Bipartite Graph!

    *Note:* automatically reduces nodes number by 1 to account for 0-based indexing in Python.

    """
    function write_edgelist(edges::Vector{Tuple{Int,Int}}, fname::String)
        file = occursin(".csv", fname) ? fname : "$(fname).csv"
        open(file, "w") do f
            for edge in edges
                println(f,"$(edge[1]-1),$(edge[2]-1)")
            end
        end
        @info "$(DTG()) - wrote edgelist $(file)"
    end

    """
        write_projected_edges(model::PyObject, fname::String)

    write out a monopartite projection from the NEMtropy computation.

    Uses 1-based index in writeout to match the user metadata and visualize in Gephi.

    # appears to be slow
    """
    function write_projected_edges(model::PyObject, fname::String)
		file = occursin(".csv", fname) ? fname : "$(fname).csv"
		
		edgelist = Set{NTuple{2,Int}}() # list of edges
		adjlist = PyDict(model.projected_rows_adj_list)
		# loop over all nodes in top layer
		for i = 0:model.n_rows-1
			if haskey(adjlist, i)
				for neigh in adjlist[i]
                    if neigh > i
					    push!(edgelist, (i+1, neigh+1))
                    end
				end
			end
		end

		# write to file
		open(file, "w") do f
            for edge in edgelist
                println(f,"$(edge[1]),$(edge[2])")
            end
        end
        @info "$(DTG()) - wrote edgelist $(file)"
	end

    """
        write_edges(A_star::Analyse.SparseArrays.AbstractSparseMatrixCSC, fname::String)

    write out a monopartite projection from the observed V-motifs.

    Uses 1-based index in writeout to match the user metadata and visualize in Gephi.
    """
    function write_edges(A_star::SparseArrays.AbstractSparseMatrixCSC, fname::String)
		file = occursin(".csv", fname) ? fname : "$(fname).csv"

		edges = Set{NTuple{2,Int}}() # list of edges
		# loop over all nodes in top layer
		for i = 1:size(A_star,1)
			v = A_star[:,i]
			for neigh in v.nzind
                # only keep smallest ones
                if neigh > i 
				    push!(edges, (i, neigh))
                end
			end
		end

		# write to file
		open(file, "w") do f
            for edge in edges
                println(f,"$(edge[1]),$(edge[2])")
            end
        end
		
        @info "$(DTG()) - wrote edgelist $(file)"
	end

    """
        write_users(udict::Dict, N_⊤::Int, fname::String)

    write out the user info for visualisation in Gephi.
    """
    function write_users(udict::Dict, Nu::Int, fname::String)
        file = occursin(".csv", fname) ? fname : "$(fname).csv"
        
        open(file, "w") do f
            println(f, "Id; Label; userid; flagged")

            for i = 1:Nu
                println(f, "$(i); $(udict[i][:Label]); $(udict[i][:userid]); $(udict[i][:flagged])")
            end
        end

       @info "$(DTG()) - wrote nodelist $(file)"
    end

    """
        write_flagged_users(udict::Dict, Nu::Int, fname::String)

    write out the user flagged status as a vector.
    """
    function write_flagged_users(udict::Dict, Nu::Int, fname::String)
        file = occursin(".csv", fname) ? fname : "$(fname).csv"
        flags = [Int(udict[i][:flagged]) for i = 1:Nu]

        open(file, "w") do f
            print(f, join(flags,","))
        end

        @info "$(DTG()) - wrote flagged nodelist vector to $(file)"
    end
    

    """
        unpickle(filename::String)

    Load up a finished computation from a pkl file.
    """
    function unpickle(filename::String)
        r = nothing
        @pywith pybuiltin("open")(filename,"rb") as f begin
            r = pickle.load(f)
        end
        return r
    end

    """
        global_loader(globalpath::String)

    function that loads up all data and prepares it for parsing by NEMtropy
    """
    function global_loader(globalpath::String; tasks::Vector{Symbol}=[:edgelist])
        @info "$(DTG()) - Sarting up in $(globalpath) for the following tasks: \n$(join(tasks,"\n\t   "))"
        @assert isdir(globalpath)
        failures = String[]
        successes = String[]
        # read subfolders
        for subfolder in filter(x -> isdir(joinpath(globalpath, x)), readdir(globalpath))
            currentpath = joinpath(globalpath, subfolder)
            @info "$(DTG()) - working on $(subfolder)"
            
            ## STEP 0 - get relevant files
            # Load up the information operations data
            msgfiles = filter(x->x[end-3:end] == ".csv" && occursin("tweets_csv_hashed",x) && !occursin(".log",x), readdir(currentpath))
            usrfiles = filter(x->x[end-3:end] == ".csv" && occursin("users_csv_hashed",x)  && !occursin(".log",x), readdir(currentpath))
            @info """$(DTG()) - files overview\n\tmsgfile(s):\n\t   $(join(msgfiles,"\n\t   "))\n\tuserfile(s):\n\t   $(join(usrfiles,"\n\t   "))"""
            try
                udf, mdf = loaddata(joinpath.(currentpath, usrfiles), 
                                    joinpath.(currentpath, msgfiles), 
                                    msgfields=[:tweetid, :userid, :hashtags])
        
                # Load up the external data
                # - external retweet files:
                rtfiles = filter(x->occursin("retweet", x) && occursin("jsonl", x), readdir(currentpath))
                @info """$(DTG()) - external RETWEET files:\n\t   $(join(rtfiles,"\n\t   "))"""
                rtfiles = joinpath.(currentpath, rtfiles)
                rtdf = vcat([readjsonlforhastags(file) for file in rtfiles]...)
        
                # - external reply files:
                rpfiles = filter(x->occursin("replies", x) && occursin("jsonl", x), readdir(currentpath))
                @info """$(DTG()) - external REPLY files:\n\t   $(join(rpfiles,"\n\t   "))"""
                rpfiles = joinpath.(currentpath, rpfiles)
                rpdf = vcat([readjsonlforhastags(file) for file in rpfiles]...)
        
                ## STEP 1 - build the bipartite network
                N_⊤, N_⊥, vertices, invertices, edges = bipartitehashtaggraph(udf, mdf, rpdf, rtdf, dropisolatedusers=true)


                ## STEP 2 - do required tasks
                if :edgelist ∈ tasks
                    #write out the edge list for the NEMtropy package
                    write_edgelist(edges, joinpath(currentpath, "bipartite_edgelist"))
                end
                if :nodelist ∈ tasks
                    write_users(vertices, N_⊤, joinpath(currentpath, "nodelist"))
                end
                if :flaggednodevector ∈ tasks
                    write_flagged_users(vertices, N_⊤, joinpath(currentpath, "flaggednodelist"))
                end
                push!(successes, subfolder)
            catch err
                @warn "Something went wrong in folder $(currentpath)\n$(err)"
                push!(failures, subfolder)
                continue
            end
        end
        
        if length(successes) > 0
            @info """$(DTG()) - following folders completed without errors:\n\t- $(join(successes,"\n\t- "))"""
        end
        if length(failures) > 0
            @warn """$(DTG()) - following folders generated errors:\n\t- $(join(failures,"\n\t- "))"""
        end

        @info "$(DTG()) - ALL DONE"
    end



    export loaddata
end
