# Test flow-based local graph clustering functions

include("LocalClusteringObjs.jl")

using MAT

mat = matread("Netscience_xy.mat")
A = mat["A"]
xy = mat["xy"]
s = 4


R = [s]
R = neighborhood(A,R,2)

S = mqi(A,R)

include("display_graph.jl")

edgebrightness = 0.8
f = display_graph(A,xy,edgebrightness)

# Show Seed set and output set
T = sink_nodes(F)
scatter!(f,[xy[R,1]],[xy[R,2]], color = :red ,markersize = 7)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :blue,markersize = 3)
