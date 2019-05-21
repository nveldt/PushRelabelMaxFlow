include("maxflow.jl")

using MAT

mat = matread("Netscience_xy.mat")
A = mat["A"]
xy = mat["xy"]

# Define a seed set
Rstart = [150, 272]
R = neighborhood(A,Rstart,2)

# Run both mqi and flowimprove and display the output
S = mqi(A,R)

include("display_graph.jl")
f = display_graph(A,xy,.8)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :green,markersize = 10)
scatter!(f,[xy[R,1]],[xy[R,2]], color = :red,markersize = 3)

S = flowimprove(A,R)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :blue,markersize = 5)
scatter!(f,[xy[R,1]],[xy[R,2]], color = :red,markersize = 3)


# Run a step of MQI
S = mqi_step(A,R,.15)

include("display_graph.jl")
f = display_graph(A,xy,.8)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :blue,markersize = 5)
scatter!(f,[xy[R,1]],[xy[R,2]], color = :red,markersize = 3)
