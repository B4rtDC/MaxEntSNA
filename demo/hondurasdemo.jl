### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 0f2f1f32-7d46-11ec-2642-8b6757aaeda2
begin
	# use proper virtual environment 
	using Pkg
	cd(joinpath(@__DIR__,".."))
	Pkg.activate(pwd())

	# Revise.jl equivalent for this notebook
	function ingredients(path::String)
		name = Symbol(basename(path))
		m = Module(name)
		Core.eval(m,
			Expr(:toplevel,
				 :(eval(x) = $(Expr(:core, :eval))($name, x)),
				 :(include(x) = $(Expr(:top, :include))($name, x)),
				 :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
				 :(include($path))))
		m
	end

	using PlutoUI
	using Logging
	using Graphs
	import LinearAlgebra: diagind

	using Plots, StatsPlots, StatsBase, LaTeXStrings, Measures
	
	using PyCall
	# import NEMtropy Python module in the notebook
	nemtropy = pyimport("NEMtropy")

	Analyse = ingredients("./src/Analyse.jl").Analyse
	Logging.disable_logging(Logging.Info)
	
	paths = ingredients(joinpath(pwd(),"demo","secrets.jl"))
	nothing
end

# ╔═╡ 37891fb2-39e1-476c-9d59-2b7a62a5a1cd
using BenchmarkTools

# ╔═╡ 6d60cd24-a1b8-4a49-bd40-bcd1ecfb5724
md"""

## Honduras demo
"""

# ╔═╡ 8e499c36-5241-4a80-8fd8-609f6b7c5b95
md"""
## Bipartite hashtag network
We have two parts in the datasource:
1. the data provide by Twitter (Information Operations Report)
2. the data that we get externally:
   - "retweet\_tweetid" for external tweet
   - "in\_reply\_to\_tweetid" for external replies
   - "quoted\_tweet\_tweetid" for external quotes (no data found)

These have been downloaded for each external dataset in the Twitter information operations datasets.

We nee to extract all hashtags and users that are required in the network. This requires a subset of fields:
* Required fields in message dataframe: `:userid, :user_screen_name`
* Required fields in message dataframe: `:tweetid, :userid, :hashtags`
"""

# ╔═╡ d8c56438-b323-413f-b50e-b73773b747b5
md"""### Raw data processing"""

# ╔═╡ da3dfb9b-99d3-4b6f-8d5d-a0f1efcec5f1
md"""
```Julia
# paths to data
datapath = "path/to/desinfodataset"
p_rp ="path/to/externalreplies"
p_rt = "path/to/externalretweets"
```
"""

# ╔═╡ 61cae8ca-dbbc-4c16-98f0-1b9e21e8e6d8


# ╔═╡ 09c60307-81af-44c6-a2ad-af6d6a8bb8f8
md"""
We start by loading the provided data:
"""

# ╔═╡ 77bf7775-798d-4be9-92b0-19d87c6900a0
begin
	udf, mdf = Analyse.loaddata(paths.datapath, msgfields=[:tweetid, :userid, :hashtags]);
	nothing
end

# ╔═╡ a8a0903c-66e8-4dcc-be16-9ad1152a8bc3
md"""
Once this is done, we can load the external data:
"""

# ╔═╡ e16ea1fb-04c7-4b42-a06a-45d0f1130b4a
begin
	# replies
	rpdf = Analyse.readjsonlforhastags(paths.p_rp)
	# retweets
	rtdf = Analyse.readjsonlforhastags(paths.p_rt)
	nothing
end

# ╔═╡ 51d40cc5-1e1b-4171-8097-bcc151e80bd0
md"""
Now we have the raw data needed to build the bipartite network.
"""

# ╔═╡ 0dfda87f-582a-4d29-b1ba-0c96b8b80211
md"""
### Network building
The function `bipartitehashtaggraph` returns the following items:
- `N_⊤`: Number of users (top/rows layer)
- `N_⊥`: Number of hashtags (bottom/colums layer)
- `vertices`: Dict with a `Int` => `Dict` mapping. The inner dict contains the following information: node kind, node label, node flag status, node userid
- `invertices`: Dict with `String` => `Int` mapping. The string matches the userid (for vertices) or the label (for hashtags).
- `edges`: array of tuples with all unique connections in the bipartite model.
"""

# ╔═╡ d85212dd-edb1-4f2d-bf76-d073b8b9f2cd
begin
	N_⊤, N_⊥, vertices, invertices, edges = Analyse.bipartitehashtaggraph(udf, mdf, rpdf, rtdf, dropisolatedusers=true)
end

# ╔═╡ a24a1a7c-ee61-40fe-af0b-666caa524c17
md"""
You can see that the number of users in the graphs varies for different setting, but the amount of tags does not. 
"""

# ╔═╡ 9ecd9ed1-682c-43ba-bd02-de3e25800ba7
begin
	Analyse.bipartitehashtaggraph(udf, mdf, rpdf, rtdf, dropisolatedusers=false)
end

# ╔═╡ 7ab6ae51-5b03-4584-8b1e-3b0f9bf2cf79
md"""
We can write out the edgelist to a file for usage with the NEMtropy script:
"""

# ╔═╡ a39816dd-9c49-4d74-b0b1-9a6b10157ed5
Analyse.write_edgelist(edges, "honduras_bipartite_edgelist")

# ╔═╡ 69594912-9bca-43a2-8db3-8e6043ec3e6b
md"""
With these elements we can build a sparse adjacency matrix for the bipartite graph.
"""

# ╔═╡ 33bac962-446e-4cae-8512-fb1be4ce60d0
begin
	I,J,V = Analyse.bipartiteadjacencymatrix(edges)
	A = Analyse.sparse(I,J,V)
end

# ╔═╡ 95271430-ff31-45d0-b3eb-51fa4c950542
md"""
We can illustrate the degree distribution of the hashtags to identify the most used ones. In turn we will remove some of these to evaluate what is the impact on the remaining network.
"""

# ╔═╡ 50db518c-e60e-432d-a723-4f16a991803d
begin
	# get degree sequence of tags
	tagdegrees = sum(A[:,N_⊤+1:end], dims = 1)'[:,1]
	# determine the histogram
	bins = minimum(tagdegrees)-1:1:maximum(tagdegrees)
	myhist = fit(Histogram, tagdegrees, bins , closed=:right)
	# get non-zero weights:
	nzinds = .!iszero.(myhist.weights)
	x = collect(minimum(tagdegrees):maximum(tagdegrees))[nzinds]
	p=scatter(x, myhist.weights[nzinds]./length(tagdegrees), scale=:log10, marker=:cross, label="", size=(600,350), left_margin=2mm, bottom_margin=3mm)
	yticks!(10. .^collect(-5:0))
	ylims!(10. ^-5,1)
	xticks!(10. .^collect(0:4))
	xlims!(1,1e4)
	xlabel!(L"k_{tag}")
	ylabel!(L"P(k_{tag})")
	savefig("hondurastagdegreedistribution.pdf")
	title!("degree distribution of hashtags")
	p
end

# ╔═╡ f9e2bcdd-70c4-429d-ac58-487b8fe6c1be
md"""
The observed v-motifs between two users i & j can be computed as follows:

"""
# *note*: this will induce a self-edge for each node, we will not take this in account in the next part of the analysis.

# ╔═╡ 9f2d502d-de72-45bc-a7c2-15f6c23cbc28
begin
	V_star = A[1:N_⊤, :] * A[:, 1:N_⊤]
	V_star[diagind(V_star)] .= 0
	V_star
end

# ╔═╡ 0297bf63-996b-4be1-9dab-32b0536ea222
md"""
### Applying the BiCM
In the previous step we obtained the projection of the bipartite graph onto the user layer, but we included all interactions and not only the statistically significant ones.

From edgelist or the adjacency matrix, we can obtain the ML parameters of the BiCM with the NEMtropy package and subsequently find the projection. For most efficient use, we use a separate Pyhon script that runs the computations for each graph.

The illustration below just shows how we can obtain the ML parameters directly in this notebook.
"""

# ╔═╡ ae4e7f1c-ac29-4e96-944a-bcce112c75eb
begin
	model = nemtropy.BipartiteGraph(edgelist=edges)
	model.solve_tool(method="fixed-point", initial_guess="degrees")
end

# ╔═╡ a10d053c-0097-4e8a-b217-c030040c9ec5
model.adj_list[0]

# ╔═╡ 15a77847-d655-4972-a7ac-b97fd388da19
md"""
The projection is obtained by running
```
model.compute_projection()
```
This runs in a seperate script for each dataset. We can continue the analysis from the stored ("pickled") results.

*a small note on switching between languages:*
* Julia uses 1-based indexing and the entire adjacency matrix
  ```
  Analyse.findnz(A[1,:]) = $(Analyse.findnz(A[1,:]))
  ```
  $(Analyse.findnz(A[1,:]))
* The NEMtropy package uses 0 based indexing and stores each layer separatly:
  ```
  model.adj_list[0]
  ```
  $(model.adj_list[0])
* We find the same value when we add the number of nodes in the top layer to this result and account for the 0-based indexing: 
  ```
  model.adj_list[0] + model.n_rows + 1
  ```
   $(convert.(Int, model.adj_list[0])[1] + model.n_rows + 1)
"""

# ╔═╡ 03eb2b8f-1ba9-401a-aaef-67dbe0118f05
md"""
### Processing the results
We start by loading up the results:
"""

# ╔═╡ a19f7ca4-9ad1-4eb9-92db-5099d32ef9cf
processed_model = Analyse.unpickle(paths.pickledpath);

# ╔═╡ bebe928f-eb57-4330-9cf5-66d071d28486
md"""
We can see the effect of the filtering when we look at some nodes:
"""

# ╔═╡ 68d821ed-df1f-4808-8a52-dc5097e3f1c7
begin
	for i in 0:10
		old_adj = V_star[:,i+1].nzind
		if length(old_adj) > 10
			try 
				new_adj = processed_model.projected_rows_adj_list[i] 
				@warn "node $(i) - before: $(length(old_adj)); after: $(length(new_adj))"
			catch
				continue
			end
		end
	end
end

# ╔═╡ ee3a6914-5417-4a0d-b8e1-3f3d86224bab
md"""
We can also look at how this leads to differences in community detection using the Leiden algorithm.
"""

# ╔═╡ 6d05fe3a-8860-4c70-a650-f2e345700c47
md"""
### Visualising the results
We can write out the edge lists before and after the removal of the non-significant links. Combined with a list of user (node) information, this allows to illustrate the before and after situation in Gephi or another tool.
"""

# ╔═╡ f55eaf45-4e53-4ae3-88db-6eb56cd43ff5
# user information
Analyse.write_users(vertices, N_⊤, "honduras_monopartite_nodelist")

# ╔═╡ 2b9d6e0e-2693-48e9-abc9-659387a13e2d
# projection on user layer (all)
Analyse.write_edges(V_star, "honduras_monopartite_edgelist")

# ╔═╡ 2960ee04-aca0-4ded-9df9-7d56c75a6315
# projection on user layer (significant ones)
Analyse.write_projected_edges(processed_model, "honduras_monopartite_edgelist_projected")

# ╔═╡ e2c4ca7e-ce67-4882-9fd6-ac6f210bb990
md"""
### Impact of removing tags
It is likely you do not capture a specific tag. we remove some tags with the highest degrees and look at the impact on the remaining network.

We will define the network from its degree sequences. 

When removing edges, we need to make sure we keep all nodes in play
"""

# ╔═╡ fa968c18-4dfc-4377-8016-f2077881a492
function write_users_reduced(udict, Nu, fname, forbidden)
	file = occursin(".csv", fname) ? fname : "$(fname).csv"
    
	open(file, "w") do f
		# first line
		println(f, "Id; Label; userid; flagged")
		# set counter
		c = 0
		for i = 1:Nu
			if i ∉ forbidden
				c += 1
				println(f, "$(c); $(udict[i][:Label]); $(udict[i][:userid]); $(udict[i][:flagged])")
			end
		end
	end

   @info "wrote nodelist $(file)"
	
end

# ╔═╡ 452446ac-14ec-444a-a462-2afe654980de
begin
	# get actual bipartite adjacency matrix
	Abip = A[1:N_⊤,N_⊤+1:end]
	# degree sequence that will be used
	k_u = sum(Abip, dims=2)[:,1]
	k_t = sum(Abip, dims=1)'[:,1]
	# copy working vectors fro removal
	k_u_red = copy(k_u)
	k_t_red = copy(k_t)
	e_red = copy(edges)
	# number of nodes to remove
	N_remove = 600
	N_disconnect = Int64[]
	# vectors for plotting
	disconnected_flagged = [0]
	disconnected_nonflaggged = [0]
	

	# iterative removal
	for i = 1:N_remove
		# find max degree
		_, maxin = findmax(k_t_red)
		# find degree of nodes connected to the tag
		icon = Abip[:,maxin]
		# decrease connected user degrees by one
		k_u_red[icon] .-= 1
		# remove all relevant edges
		deleteat!(e_red, findall(x-> x[2]==maxin+N_⊤, e_red))
		# check for disconnected nodes and add to subset
		negatives = findall(x->x<=0, k_u_red)
		for u in negatives
			if u ∉ N_disconnect
				push!(N_disconnect,u)
			end
		end
		# set tag degree to zero
		k_t_red[maxin] = 0
		# plotting variables
		nf = count(x-> vertices[x][:flagged], negatives)
		push!(disconnected_flagged, nf)
		push!(disconnected_nonflaggged, length(negatives)-nf)

		rat = [0.0005; 0.001;0.01]
		if i in map(x->round(x*N_⊥), rat)
			# write out edgelist
			Analyse.write_edgelist(e_red, joinpath(paths.reducemodelpath,"removedtags", "edgelist$(round(i/N_⊥, digits=4))removed"))
			# write out nodelist
			write_users_reduced(vertices, N_⊤, joinpath(paths.reducemodelpath,"removedtags", "nodelist$(round(i/N_⊥, digits=4))removed"),N_disconnect )
		end
		
	end
end

# ╔═╡ 6a5694fb-1885-4179-8940-ce04f3802aab
begin
	nnf = count(x->vertices[x][:flagged], 1:N_⊤)
	xvals = collect(0:N_remove) ./ N_⊥*100
	ppp= plot(xvals, disconnected_nonflaggged./(N_⊤-nnf)*100, label="non-flagged users", legend=:topleft)#, marker=:circle)
	plot!(xvals, disconnected_flagged./nnf*100, label="flagged users")#, marker=:circle)
	plot!(xvals, (disconnected_flagged .+ disconnected_nonflaggged)./N_⊤*100, label="total", color=:black,line=:dash)#, marker=:circle
	xlabel!("% of high degree hashtags removed")
	ylabel!("% of disconnected nodes")
	xlims!(0,maximum(xvals))
	ylims!(0,80)
	savefig("hondurasremoval.pdf")
	ppp
end

# ╔═╡ 9afb55c9-f446-46ba-a1d0-7c35dd1054c3


# ╔═╡ Cell order:
# ╟─6d60cd24-a1b8-4a49-bd40-bcd1ecfb5724
# ╠═0f2f1f32-7d46-11ec-2642-8b6757aaeda2
# ╟─8e499c36-5241-4a80-8fd8-609f6b7c5b95
# ╟─d8c56438-b323-413f-b50e-b73773b747b5
# ╟─da3dfb9b-99d3-4b6f-8d5d-a0f1efcec5f1
# ╠═61cae8ca-dbbc-4c16-98f0-1b9e21e8e6d8
# ╟─09c60307-81af-44c6-a2ad-af6d6a8bb8f8
# ╠═77bf7775-798d-4be9-92b0-19d87c6900a0
# ╟─a8a0903c-66e8-4dcc-be16-9ad1152a8bc3
# ╠═e16ea1fb-04c7-4b42-a06a-45d0f1130b4a
# ╟─51d40cc5-1e1b-4171-8097-bcc151e80bd0
# ╟─0dfda87f-582a-4d29-b1ba-0c96b8b80211
# ╠═d85212dd-edb1-4f2d-bf76-d073b8b9f2cd
# ╟─a24a1a7c-ee61-40fe-af0b-666caa524c17
# ╠═9ecd9ed1-682c-43ba-bd02-de3e25800ba7
# ╟─7ab6ae51-5b03-4584-8b1e-3b0f9bf2cf79
# ╠═a39816dd-9c49-4d74-b0b1-9a6b10157ed5
# ╟─69594912-9bca-43a2-8db3-8e6043ec3e6b
# ╠═33bac962-446e-4cae-8512-fb1be4ce60d0
# ╟─95271430-ff31-45d0-b3eb-51fa4c950542
# ╟─50db518c-e60e-432d-a723-4f16a991803d
# ╟─f9e2bcdd-70c4-429d-ac58-487b8fe6c1be
# ╠═9f2d502d-de72-45bc-a7c2-15f6c23cbc28
# ╟─0297bf63-996b-4be1-9dab-32b0536ea222
# ╠═ae4e7f1c-ac29-4e96-944a-bcce112c75eb
# ╠═a10d053c-0097-4e8a-b217-c030040c9ec5
# ╟─15a77847-d655-4972-a7ac-b97fd388da19
# ╟─03eb2b8f-1ba9-401a-aaef-67dbe0118f05
# ╠═a19f7ca4-9ad1-4eb9-92db-5099d32ef9cf
# ╟─bebe928f-eb57-4330-9cf5-66d071d28486
# ╠═68d821ed-df1f-4808-8a52-dc5097e3f1c7
# ╟─ee3a6914-5417-4a0d-b8e1-3f3d86224bab
# ╟─6d05fe3a-8860-4c70-a650-f2e345700c47
# ╠═f55eaf45-4e53-4ae3-88db-6eb56cd43ff5
# ╠═2b9d6e0e-2693-48e9-abc9-659387a13e2d
# ╠═2960ee04-aca0-4ded-9df9-7d56c75a6315
# ╟─e2c4ca7e-ce67-4882-9fd6-ac6f210bb990
# ╠═452446ac-14ec-444a-a462-2afe654980de
# ╠═fa968c18-4dfc-4377-8016-f2077881a492
# ╟─6a5694fb-1885-4179-8940-ce04f3802aab
# ╟─9afb55c9-f446-46ba-a1d0-7c35dd1054c3
# ╟─37891fb2-39e1-476c-9d59-2b7a62a5a1cd
