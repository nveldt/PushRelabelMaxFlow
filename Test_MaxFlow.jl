# Test the push-relabel max flow code

include("maxflow.jl")

using MAT

mat = matread("Netscience_xy.mat")
A = mat["A"]
xy = mat["xy"]
n = size(A,1)
s = 4
R,vals = findnz(A[:,s])
R,vals = findnz(A[:,R])

# Set up a weighted max-flow/min-cut problem via the standard auxiliary graph
# construction for SimpleLocal/FlowSeed/FlowImprove.
alpha = .7
rvec = zeros(n)
rvec[R] .= 1
svec = alpha*rvec
epsilon = .01
tvec = alpha*epsilon*(ones(n)-rvec)

C = [spzeros(1,1) sparse(svec') spzeros(1,1);
    spzeros(n,1) A sparse(tvec);
    spzeros(1,1) spzeros(1,n) spzeros(1,1)]

flowtolerance = 1e-10
F = maxflow(C,1, n+2, flowtolerance)
E = cut_edges(F)
S = source_nodes(F).-1
S = setdiff(S,0)

include("display_graph.jl")
f = display_graph(A,xy,.8)

# Color source set
scatter!(f,[xy[R,1]],[xy[R,2]], color = :red,markersize = 5)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :blue,markersize = 3)

sorig = source_nodes(F,flowtolerance)
torig = sink_nodes(F,flowtolerance)

# Check cut value and flow value, which should be equal
@show sum(C[sorig,torig])
@show F.flowvalue
