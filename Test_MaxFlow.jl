include("maxflow.jl")

using MAT

mat = matread("Netscience_xy.mat")
A = mat["A"]
xy = mat["xy"]
s = 4
t = 67

F = maxflow(A,s,t)

include("display_graph.jl")

f = display_graph(A,xy,.8)

# Highlight the source and sink node.
S = F.source_nodes
scatter!(f,[xy[F.s,1]],[xy[F.s,2]], color = :yellow,markersize = 10,markershape = :o)
scatter!(f,[xy[F.t,1]],[xy[F.t,2]], color = :yellow,markersize = 10,markershape = :diamond)

# Color source set and sink set
T = sink_nodes(F)
scatter!(f,[xy[S,1]],[xy[S,2]], color = :blue,markersize = 5)
scatter!(f,[xy[T,1]],[xy[T,2]], color = :red,markersize = 5)

## Highlight the cut edges in magenta
edgeset = cut_edges(F)
ei = edgeset[:,1]
ej = edgeset[:,2]
lx = [xy[ei,1]';xy[ej,1]';NaN*ones(1,length(ei))]
ly = [xy[ei,2]';xy[ej,2]';NaN*ones(1,length(ei))]
for i = 1:size(edgeset,1)-1
    plot!(lx[:,i],ly[:,i],color = :magenta, linewidth = 1)
end
i = size(edgeset,1)
plot!(lx[:,i],ly[:,i],color = :magenta, linewidth = 1)
