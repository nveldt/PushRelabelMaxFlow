using MatrixNetworks
using SparseArrays

mutable struct stFlow
    value::Float64 # gives you the max-flow value
    source_nodes::Vector{Int64} # give the indices of the nodes attached to the source
    C::SparseMatrixCSC # gives the original capacity matrix
    F::SparseMatrixCSC # gives the values of the flows on each edge
    s::Int64  # index of source node
    t::Int64 # index of sink node
end

"""
maxflow

Given a sparse matrix A representing a weighted and possibly directed graph,
a source node s, and a sink node t, return the maximum s-t flow.

Returns F, which is of type stFlow.
"""
function maxflow(B::Union{SparseMatrixCSC,MatrixNetwork},s::Int,t::Int)

    # The code actually assumes a SparseMatrixCSC input
    if typeof(B) <: SparseMatrixCSC
    else
        B = sparse(B)
    end

    N = size(B,1)

    # Extract weights from source s to non-terminal nodes,
    # and from non-terminal nodes to sink node t
    sWeights = Array(B[s,:])
    tWeights = Array(B[:,t])
    NonTerminal = setdiff(collect(1:N),[s t])

    sWeights = sWeights[NonTerminal]
    tWeights = tWeights[NonTerminal]

    # Extract the edges between non-terminal nodes
    A = B[NonTerminal,NonTerminal]

    # A = the matrix of capacities for all nodes EXCEPT the source and sink
    # sWeights = a vector of weights for edges from source to non-terminal nodes
    # tWeights = vector of weights from non-terminal nodes to the sink node t.

    # This is the map from the original node indices to the rearranged
    # version in which the source is the first node and the sink is the last
    Map = [s; NonTerminal; t]

    # Directly set up the flow matrix
    C = [spzeros(1,1) sparse(sWeights') spzeros(1,1);
         sparse(sWeights) A sparse(tWeights);
         spzeros(1,1) sparse(sWeights') spzeros(1,1)]

    # Allocate space for the flow we will calculate
    # In a flow problem, we will eventually need to send flow the reverse
    # direction, so it's important to allocate space for F[i,j] if C[j,i] is an
    # edge, even if C[i,j] is not directed
    Cundir = C+C'
    F = SparseMatrixCSC(N,N,Cundir.colptr,Cundir.rowval,zeros(length(Cundir.rowval)))
    ExcessNodes = vec(round.(Int64,findall(x->x!=0,sWeights).+1))

    # Initialize the Preflow and the excess vector
    for v = ExcessNodes
        F[1,v] = C[1,v]
        F[v,1] = -C[1,v]
    end
    excess = [0;sWeights;0]
    source_nodes, FlowMat, value = Main_Push_Relabel(C,F,ExcessNodes,excess)

    smap = sortperm(Map)
    F = stFlow(value, sort(Map[source_nodes]),C[smap,smap],FlowMat[smap,smap],s,t)
    return F
end

"""
This maxflow code assumes that A represents the adjacencies between
non-terminal nodes. Edges adjecent to source node s and sink node t
are given by vectors svec and tvec.

This code sets s as the first node, and t as the last node.
"""
function maxflow(A::Union{SparseMatrixCSC,MatrixNetwork},svec::Vector{Float64},tvec::Vector{Float64})

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end

    # Directly set up the flow matrix
    C = [spzeros(1,1) sparse(sWeights') spzeros(1,1);
         sparse(sWeights) A sparse(tWeights);
         spzeros(1,1) sparse(sWeights') spzeros(1,1)]

    # Allocate space for the flow we will calculate
    # In a flow problem, we will eventually need to send flow the reverse
    # direction, so it's important to allocate space for F[i,j] if C[j,i] is an
    # edge, even if C[i,j] is not directed.
    Cundir = C+C'
    F = SparseMatrixCSC(N,N,Cundir.colptr,Cundir.rowval,zeros(length(Cundir.rowval)))
    ExcessNodes = vec(round.(Int64,findall(x->x!=0,sWeights).+1))

    # Initialize the Preflow and the excess vector
    for v = ExcessNodes
        F[1,v] = C[1,v]
        F[v,1] = -C[1,v]
    end
    excess = [0;sWeights;0]
    source_nodes, FlowMat, value = Main_Push_Relabel(C,F,ExcessNodes,excess)

    F = stFlow(value,source_nodes,FlowMat,s,t)
end

maxflow(A::Union{SparseMatrixCSC,MatrixNetwork},svec::Vector{Int64},tvec::Vector{Int64}) =
    maxflow(A,float(svec),float(tvec))


flow(F::stFlow) =
    F.value

"""
Given a flow, stored in an stFlow object, return the set of nodes attached to
the source
"""
function source_nodes(F::stFlow)
    # Run a bfs from the sink node. Anything with distance
    # n is disconnected from the sink. Thus it's part of the minimium cut set
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,F.t)
    S = Vector{Int64}()
    for i = 1:n
        if finalHeight[i] == n
            push!(S,i)
        end
    end

    # Sanity checks
    @assert(~in(F.t,S))
    @assert(in(F.s,S))

    return S
end

"""
Given a flow, stored in an stFlow object, return the set of nodes attached to
the sink
"""
function sink_nodes(F::stFlow)
    # Run a bfs from the sink node. Anything with distance < n is sink-attached.
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,F.t)
    T = Vector{Int64}()
    for i = 2:n
        if finalHeight[i] < n
            push!(T,i)
        end
    end

    # Sanity checks
    @assert(in(F.t,T))
    @assert(~in(F.s,T))

    return T
end

"""
Gives the cut as a list of edges.
"""
function cut_edges(F::stFlow)
    # Run a bfs from the sink node to get source and sink sets
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,F.t)
    T = Vector{Int64}()
    S = Vector{Int64}()
    for i = 1:n
        if finalHeight[i] < n
            push!(T,i)
        else
            push!(S,i)
        end
    end

    I,J,V = findnz(F.C[S,T])
    return [S[I] T[J]]

end


## FlowSeed based methods
#
# S = mqi(A,R)
# S = flowimprove(A,R)
# S = simplelocal(A,R,delta)
#
#
# S,cond = mqi_step(A,R,phi) # one step of MQI

# CURRENTLY NOT IMPLEMENTED
# Have MQI and other methods return a flow (currently it just returns a cut set)
# Have not implemented general weighted version.
#
# e.g. implement:
#
# S,cond,F = mqi_maxflow(A,R,phi) # these give back intermediate information...
# # mqi(A,wvec,R) # give an arbitrary vertex weight set... = pi in FlowImprove ?


# Main_Push_Relabel returns a maximum flow F and the min s-t cut set S for the
# flow graph C.
#
# C = the capacity matrix for the flow problem.
#   This code assumes node 1 is the source, and node n is the sink.
#   the preflow immediately pushes all flow from the source to create an
#   excess on nodes in the graph.
#
# F = an initial flow. It can be initialize to zero.
#
# excess = the vector of excess values at the start of the algorithm. If F = 0,
#   this is the vector of edge capacities from the implicit source to the graph.
#   If F != 0, then it's the excess from a previous run of the algorithm
function Main_Push_Relabel(C::SparseMatrixCSC,
    F::SparseMatrixCSC,ExcessNodes::Array{Int64},excess::Array{Float64})

    # here, n includes only one terminal node, the sink
    n = size(C,1)

    height = zeros(Int64,n)      # label/height of each node
    inQ = zeros(Bool,n)          # list whether or not nodes are in the queue

    # Store adjacency list. Because flow can be sent either direction on an
    # arc during the course of the algorithm, it's important to list all neighbors
    # or each node, counting both incoming and outgoing edges
    Neighbs,d = ConstructAdj(C+C',n)

    # We will maintain a queue of active nodes.
    #   An actual queue implementation is available in the DataStructures.jl
    #   Julia package. The performane is nearly identical (and in some cases
    #   slightly slower), thus to minimize dependency on outside packages, we
    #   just use a Vector, rather than an actual implemented Queue
    Queue = Vector{Int64}()

    # Start by saturating edges from source to its neighbors
    # All nodes with nonzero excess are the first to be processed
    for v = ExcessNodes
        push!(Queue,v)
    end
    inQ[ExcessNodes] .= true

    # count the number of nodes that have been relabeled
    relabelings::Int64 = 0

    height = relabeling_bfs(C,F,n)     # compute initial distance from sink
    # In the code and comments, height = distance from sink = label of node

    # Continue until the queue no longer contains any active nodes.
    while length(Queue) > 0

        u = pop!(Queue)     # Select a new active node

        inQ[u] = false      # Take it out of the queue

        # discharge flow through node u
        relabelings += discharge!(C,F,Queue,u,Neighbs[u],height,excess,n,d[u],inQ)

        # if u is still active, put it back into the queue
        if excess[u] > 0 #&& height[u] < n
            prepend!(Queue,u)
            inQ[u] = true
        end

        # Global relabeling heuristic for push-relabel algorithm.
        # This periodically recomputes distances between nodes and the sink
        if relabelings == n
            relabelings = 0
            dist = relabeling_bfs(C,F)
            height = dist
        end

    end

    # Compute final distances from sink using BFS. Anything with distance
    # n is disconnected from the sink. Thus it's part of the minimium cut set
    finalHeight = relabeling_bfs(C,F,n)
    S = Vector{Int64}()
    push!(S,1)          # Include the source node
    for i = 2:n
        if finalHeight[i] == n
            push!(S,i)
        end
    end

    mflow = excess[n]     # the excess at the sink equals the maximum flow value

    return S, F, mflow

end

# Discharege operation: pushes flow away from node u across admissible edges.
# If excess[u] > 0 but no admissible edges exist, we relabel u.
function discharge!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    excess::Array{Float64},n::Int64,du::Int64,inQ::Array{Bool})

    vLocal::Int64 = 1           # Start at the first neighbor of node u
    hu = height[u]
    relabeled = 0

    # As long as there is excess at node u and there is another neighbor to explore...
    while excess[u] > 0 && vLocal <= du

            # ...grab the next neighbor of node u
            v = uNeighbs[vLocal]

            # ... if edge (u,v) is admissible, push more flow.
            # Otherwise, move to the next neighbor of u
            if hu > height[v] && C[u,v] - F[u,v] > 0
                pushflow!(C,F,Queue,u,v,excess,height,inQ,n)
                vLocal += 1
            else
                vLocal += 1
            end
    end

    # if we needed to visit every neighbor of u, we must relabel it,
    # so that at least one admissible edge is created
    if vLocal > du
        relabeled = 1
        relabel!(C,F,Queue,u,uNeighbs,height,du,n)
    end

    return relabeled
end

# Relabel sets the label/height of node u to be equal to the minimum label
# such that an admissible edge exists. An edge (u,v) is admissible if
# height[u] = height[v] + 1
function relabel!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    du::Int64,n::Int64)
   # find smallest new height making a push possible, if such a push is possible

   min_height = Inf
   # search through the neighbors of u
   # and relabel so that height[u] = height[v] + 1 for some v in the neighborhood
   for vLocal = 1:du
       v = uNeighbs[vLocal]
       if C[u,v] - F[u,v] > 0
           min_height = min(min_height, height[v])
           height[u] = min_height + 1
       end
   end

end

# Push flow from an active node u to a node v via an admissible edge (u,v)
function pushflow!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int},u::Int64,v::Int64,excess::Array{Float64},height::Array{Int64},
    inQ::Array{Bool},n::Int64)

    send = min(excess[u], C[u,v] - F[u,v])
    F[u,v] += send
    F[v,u] -= send
    excess[u] -= send
    excess[v] += send

    # If v isn't in the queue, isn't the sink, isn't the course,
    # and is active, then add it to the Queue
    if ~inQ[v] && v < n && v > 1 #&& height[v] < n
        prepend!(Queue,v)
        inQ[v] = true
    end
end

# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    d = zeros(Int64,n)
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
        d[i] = ci[i+1]-ci[i]
    end

    # d is the number of neighbors. This is the unweighted degree,
    # but note importantly that if the original graph is weighted this is
    # not the same as the degree vector d we will sometimes use
    return Neighbs, d

end

# Given initial capacity matrix C and flow matrix F, compute the distance
# from each node to the specified "start" node.
# Start defaults to node n, which is assumed to be the sink node
function relabeling_bfs(C::SparseMatrixCSC,F::SparseMatrixCSC,start::Int64=0)
    # To avoid subtraction cancellation errors that may have ocurred when pushing
    # flow, when computing a bfs we round edges to zero if they are under
    # a certain tolerance
    Cf = round.((C-F),digits = 6)
    n = size(Cf,1)

    if start == 0
        start = n
    end

    rp = Cf.colptr
    ci = Cf.rowval

    N=length(rp)-1

    d = n*ones(Int64,N)
    sq=zeros(Int64,N)
    sqt=0
    sqh=0 # search queue and search queue tail/head

    # start bfs at the node "start"
    u = start
    sqt=sqt+1
    sq[sqt]=u
    d[u]=0
    while sqt-sqh>0
        sqh=sqh+1
        v=sq[sqh] # pop v off the head of the queue
        for ri=rp[v]:rp[v+1]-1
            w=ci[ri]
            if d[w] > n-1
                sqt=sqt+1
                sq[sqt]=w
                d[w]= d[v]+1
            end
        end
    end

    return d
end

## FlowSeed-based code
function mqi(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64})
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    n = size(A,2)
    epsilon = Inf
    numR = length(R)
    emptyvec = zeros(numR)
    S, cond = FlowSeed(A,R,epsilon,emptyvec,emptyvec,false,true)
    return S
end

# One step of MQI
function mqi_step(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},alpha::Float64)
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    epsilon = Inf
    numR = length(R)
    emptyvec = zeros(numR)
    FlowSeed_Step(A,R,epsilon,emptyvec,emptyvec, alpha,false,true)
end

# FlowImprove
function flowimprove(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},localflag::Bool=true)
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    n = size(A,2)
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    epsilon = volR/(volA-volR)
    numR = length(R)
    emptyvec = zeros(numR)
    S, relcond = FlowSeed(A,R,epsilon,emptyvec,emptyvec,true,localflag)
    return S
end

# FlowImprove
function flowimprove_step(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},alpha::Float64,localflag::Bool=true)
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    n = size(A,2)
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    epsilon = volR/(volA-volR)
    numR = length(R)
    emptyvec = zeros(numR)
    FlowSeed_Step(A,R,epsilon,emptyvec,emptyvec,alpha,true,localflag)
end

# SimpleLocal
function simplelocal(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},delta::Float64)
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    n = size(A,2)
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    epsilon = volR/(volA-volR)+delta
    numR = length(R)
    emptyvec = zeros(numR)
    S, relcond = FlowSeed(A,R,epsilon,emptyvec,emptyvec,true,true)
    return S
end

# Simplelocal step
function simplelocal_step(A::Union{SparseMatrixCSC,MatrixNetwork},
    R::Vector{Int64},delta::Float64,alpha::Float64,localflag::Bool=true)
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    n = size(A,2)
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    epsilon = volR/(volA-volR)+delta
    numR = length(R)
    emptyvec = zeros(numR)
    FlowSeed_Step(A,R,epsilon,emptyvec,emptyvec,alpha,true,localflag)
    return S
end

# FlowSeed
# with simplified parameters
function FlowSeed(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},
    epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    relcondFlag::Bool= true,localFlag::Bool=true)

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    n = size(A,1)

    # Find one-hop neighbors of R, and get the complement set
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] .= 0
    Rc = findall(x->x!=0,inRc)             # complement of R

    S, relcond = FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,relcondFlag,localFlag)
    return S, relcond
end

# The most general function, FlowSeed, which minimizes a localized variant of
# conductance which penalizes the exclusion of seed nodes from the output set.
#
# Parameters:
#
#   A = adjacency matrix for a graph
#
#   R = node indices for a seed set,
#   Rn = immediate neighbors of R
#   Rc = complement set of R
#
#   epsilon = locality parameter
#
#   pR = a length(R) vector with penalties on exluding seed nodes in R from
#       the output set. pR[i] is the penalty or excluding R[i] from the output
#
#   RinS = a length(R) zero-one matrix indicating which nodes in R are stricly
#           required to be in the output set
#
#  relcondFlag = a boolean flag indicating whether to compute the relative
#                conductance score or the exact conductance score for each
#                intermediate improved set. Choosing false (i.e. updating with
#                exact conductance) will sometimes lead to fewer iterations and
#                lower conductance output, but will not actually minimize the
#                relative conductance or seed penalized conductance.
#
#  localFlag = a boolean flag indicating whether or not to use the local
#               computations. If volR is large and epsilon is small, in some
#               cases it may be better for the subroutine to perform one
#               global caluculations that multiple "local" computations.
#
# d = weighted degree vector of the graph
#
#
# volA, volR = volumes of the entire graph and seed set respectively
function FlowSeed(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},
    Rn::Vector{Int64},Rc::Vector{Int64},epsilon::Float64,pR::Array{Float64},
    RinS::Array{Float64},d::Array{Float64},volA::Float64=0.0,volR::Float64=0.0,
    relcondFlag::Bool= true,localFlag::Bool=true)

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    fR = volR/(volA - volR)
    if epsilon < fR
        println("Locality parameter epsilon was set to small. Setting it to
                    lower bound of $fR. Computations will not be local.")
        epsilon = fR
    end

    n = size(A,1)

    if localFlag
        if volA*epsilon/volR < 10
            println("Note that vol(R)/epsilon = O(vol(G)).
            For these parameters \nit may be faster to run the algorithm
            without the locality setting.")
        end
    end

    # Call nodes that must be S the "strong seed nodes"
    localStrong = findall(x->x!=0,RinS)

    StrongSeeds = R[localStrong]
    numstrong = length(StrongSeeds)

    # If something is marked as a strong seed, put an infinite penalty
    # on excluding it from the output set
    pR[localStrong] .= Inf

    # Conductance of R
    Stats = set_stats(A,R,volA)
    alphaCurrent = Stats[4]
    # Conductance of R is same as localized seed penalized conductance of R
    # alpha2 = cutval(A,R,R,d,1.0,epsilon,volA,pR,RinS)
    # println("$alpha2, $alphaCurrent")


    println("\nEpsilon = $epsilon");
    println("There are $numstrong strong seed nodes.")
    println("The full seed set has conductance $alphaCurrent ");
    println("-------------------------------------------------------")
    BestS = R
    alph0 = 2
    alphaBest = alphaCurrent

    source = zeros(n)
    sink = zeros(n)
    dr = d[R]
    drc = d[Rc]

    while alphaCurrent < alph0

        # Prepare source-side and sink-side edge weights for the augmented
        # local flow graph
        # Seed nodes have an edge to the source of the following weight
        source[R] = alphaCurrent*(pR .+ 1).*dr

        # Non-seed nodes have an edge to the sink
        sink[Rc] = alphaCurrent*epsilon*drc

        # Compute the new min s-t cut
        if localFlag
            # Do it by repeatedly solving smaller problems, starting
            # by looking at the immediate neighbors Rn
            S = LocalPushRelabel(A,R,source,sink,Rn)
        else
            # Run a single min-cut computation on the whole graph
            S = NonLocalPushRelabel(A,R,source,sink)
        end

        if length(S) > 0 && length(S) < n

            # Check stats for new set
            if relcondFlag
                alphaS = cutval(A,S,R,d,1.0,epsilon,volA,pR,RinS)
            else
                Stats = set_stats(A,S,volA)
                alphaS = Stats[4]
            end

            if alphaS < alphaCurrent
                numS = size(S,1)
                ra = round(alphaS,digits =4)
                println("Improvement found: R-Conductance = $ra, Size = $numS")
                BestS = S
                alphaBest = alphaS
            end

        else
            alphaS = alphaCurrent
        end

        alph0 = alphaCurrent
        alphaCurrent = alphaS

    end

    SL = BestS
    sizeSL = length(SL)
    cond = alphaBest
    println("------------------------------------------------------")
    println("Final Answer: Conductance = $cond, Size = $sizeSL ")

    return SL, cond
end


# Run a single step of FlowSeed, for a fixed alpha
function FlowSeed_Step(A::Union{SparseMatrixCSC,MatrixNetwork},R::Vector{Int64},
    epsilon::Float64,pR::Array{Float64}, RinS::Array{Float64}, alpha::Float64,
    relcondFlag::Bool= true,localFlag::Bool=true)

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    n = size(A,1)

    # Find one-hop neighbors of R, and get the complement set
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] .= 0
    Rc = findall(x->x!=0,inRc)             # complement of R

    fR = volR/(volA - volR)
    if epsilon < fR
        println("Locality parameter epsilon was set to small. Setting it to lower bound of $fR. Computations will not be local.")
        epsilon = fR
    end

    n = size(A,1)

    # Call nodes that must be S the "strong seed nodes"
    localStrong = findall(x->x!=0,RinS)

    StrongSeeds = R[localStrong]
    numstrong = length(StrongSeeds)

    # If something is marked as a strong seed, put an infinite penalty
    # on excluding it from the output set
    pR[localStrong] .= Inf

    BestS = R
    alphaBest = alpha
    alphaCurrent = alpha

    source = zeros(n)
    sink = zeros(n)
    dr = d[R]
    drc = d[Rc]

    # Prepare source-side and sink-side edge weights for the augmented
    # local flow graph
    # Seed nodes have an edge to the source of the following weight
    source[R] = alphaCurrent*(pR .+ 1).*dr

    # Non-seed nodes have an edge to the sink
    sink[Rc] = alphaCurrent*epsilon*drc

    # Compute the new min s-t cut
    if localFlag
        # Do it by repeatedly solving smaller problems, starting
        # by looking at the immediate neighbors Rn
        S = LocalPushRelabel(A,R,source,sink,Rn)
    else
        # Run a single min-cut computation on the whole graph
        S = NonLocalPushRelabel(A,R,source,sink)
    end

    # If things improved, return this set, otherwise we'll return R
    if length(S) > 0 && length(S) < n

        # Check stats for new set
        if relcondFlag
            alphaS = cutval(A,S,R,d,1.0,epsilon,volA,pR,RinS)
        else
            Stats = set_stats(A,S,volA)
            alphaS = Stats[4]
        end

        if alphaS < alphaCurrent
            BestS = S
            alphaBest = alphaS
        end
    end

    return BestS
end

# Starting from a set of seed nodes R, do a breadth first search to get
# a k-hop neighborhood of R
function neighborhood(A::Union{SparseMatrixCSC,MatrixNetwork},R::Array{Int64},k::Int64)

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    rp = A.rowval
    ci = A.colptr
    n = size(A,1)

    eS = zeros(n)
    eS[R] .= 1

    # For node i, the neighbors of i are rp[ci[i]:ci[i+1]-1]
    for i = R
        neighbs = rp[ci[i]:ci[i+1]-1]
        eS[neighbs] .= 1
    end

    # This could be more efficient, but recursively calling won't take too long
    # as long as k isn't too large
    if k == 1
        return findall(x->x!=0,eS)
    else
        return neighborhood(A,findall(x->x!=0,eS),k-1)
    end

end

# For a set S in a graph with adjacency matrix A, return some information about
# S including its conductance, number of interior edges, volume, and cut.
function set_stats(A::Union{SparseMatrixCSC,MatrixNetwork},
    S::Vector{Int64},volA::Float64)

    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end
    if volA == 0.0
        volA = sum(A.nzval)
    end

    if length(S) == size(A,1)
        # then we have an indicator vector
        S = findall(x->x!=0,eS)
        AS = A[:,S]
    else
        # then we have a subset
        @assert(minimum(S) >= 1)
        @assert(maximum(S) <= size(A,1))
        AS = A[:,S]
    end

    vol = sum(AS.nzval);
    SAS = AS[S,:]
    edges = sum(SAS.nzval);
    cut = vol-edges

    cond = cut/minimum([vol,volA-vol]);

    return cut, vol, edges, cond

end

# Compute the s-t cut score corresponding to a set S, in an augmented graph
# with source and sink node
function cutval(A::SparseMatrixCSC,S::Vector{Int64},
    R::Vector{Int64},d::Array{Float64,2},alpha::Float64,epsilon::Float64,
    volA::Float64,pR::Array{Float64},RinS::Array{Float64})

    n = size(A,1)
    if volA == 0.0
        volA = sum(A.nzval)
    end

    strongR = R[findall(x->x!=0,RinS)]
    @assert(length(setdiff(strongR,S)) == 0)    # S should contain strongR

    @assert(minimum(S) >= 1)
    @assert(maximum(S) <= size(A,1))
    AS = A[:,S];


    volS = sum(AS.nzval);
    SAS = AS[S,:]
    edges = sum(SAS.nzval);
    cutS = volS-edges

    volR = sum(d[R])

    # penalty vector, should only be nonzero for R nodes
    penalty = zeros(n)
    penalty[R] = pR.*d[R]

    RS = intersect(R,S)
    volRS = sum(d[RS])
    RnotinS = setdiff(R,RS)   # the set of nodes in R that aren't in S
    pRnotinS = sum(penalty[RnotinS])    # the penalty for excluding R nodes from A

    cutScore = cutS - alpha*volRS + alpha*volR + alpha*epsilon*(volS-volRS) + alpha*pRnotinS

    @assert(cutScore >= 0)

    relcond = cutS/(volRS - epsilon*(volS-volRS) - pRnotinS)

    return relcond
end

# LocalPushRelabel: computes the minimumn s-t cut for a flow graph in strongly-local
#               time. It repeatedly solves localized min-cut problems.
#
# Input Parameters:
#
# A = a symmetric matrix representing an undirected graph. It can be weighted.
#
# R = a list of nodes that share an edge with the source node
#
# sWeights and tWeight store the nonnegative weight of each node to the source
# and sink. For node i, exactly one of sWeights[i] and tWeights[i] is nonzero
#
# Rn = a list of nodes not in R that neighbor a node in R
function LocalPushRelabel(A::SparseMatrixCSC,R::Vector{Int64},
    sWeights::Array{Float64},tWeights::Array{Float64},Rn::Array{Int64})

    timer = 0.0

    n = size(A,1)
    rp = A.rowval
    ci = A.colptr

    # Now we want to locally compute maximum flows
    # C = indices of "complete" nodes in the local graph L, which are nodes
    #   whose degree in the local graph equals the degree in the global graph.
    # I = local indices of nodes that are in L, but not complete. These do
    #      share edges with one another, but only with complete nodes.

    # Initialize the complete set to be the set of nodes adjacent to the source
    C_global = R
    I_global = Rn           # everything else is incomplete
    Ac = A[C_global,:]      # set of edges from the complete set to the rest of the graph

    # We will maintain a map from indices in a local subgraph, to global indices in A.
    # These don't include the sink node in the flow graph, we are considering
    # just a growing local subgraph of A
    Local2Global = [C_global; I_global]
    # Node i in the local graph corresponds to the node with index
    # Local2Glocal[i] in the global graph A

    # Number of nodes in the local graph
    Lsize = length(Local2Global)

    # Indices, in the local graph, of complete and incomplete nodes
    C_local = collect(1:length(R))
    I_local = collect(length(R)+1:Lsize)
    numI = length(I_global)     # number of incomplete nodes

    # Build the initial local graph

    AcToI = Ac[:,I_global]     # edges between complete and incomplete nodes
    AcToc = Ac[:,C_global]     # edges between complete nodes
    L = [AcToc AcToI;
        AcToI' spzeros(numI,numI)]   # adjacency matrix for local graph

    # We distinguish between L the "local graph", and Lf, the "local flow graph"
    # which additionally contains the sink node t (as node 1).

    # In the local flow graph, each non-terminal node has either a source-side
    # or sink-side edge.
    tToL = reshape(tWeights[Local2Global],Lsize)
    sToL = reshape(sWeights[Local2Global],Lsize)

    # By adding the edges to the sink,
    # we transform the local graph L into the local flow graph Lf

    Lf = [spzeros(1,1) sparse(tToL');
         sparse(tToL) L]

    # Initialize the flow matrix; allocate space for non-zero flow values
    nLf = size(Lf,1)
    F = SparseMatrixCSC(nLf,nLf,Lf.colptr,Lf.rowval,zeros(length(Lf.rowval)))
    # Find the minimum cut for Lf.
    #
    # The first node in Lf is the sink, so offset indices of R by 1.
    start = time()
    S_local,F,excess = FlowSeed_Push_Relabel(Lf,F,collect(2:length(R)+1),[0; sToL])
    timer += time()-start

    # F is a preflow that is returned. It is NOT the maximum flow for Lf.
    # S is the set of nodes in the min s-t cut of Lf. S_local are the local
    # indices in L, (not the indices in A or Lf)

    # We "expand" L around nodes in S that were previously "incomplete"
    E_local = setdiff(S_local,C_local)         # Nodes to expand around
    E_global = Local2Global[E_local]           # their global indices

    # Keep track of which nodes are in the local graph L
    inL = zeros(Bool,n)
    inL[Local2Global] .= true

    # As long as we have new nodes to expand around, we haven't yet found
    # the global minimum s-t cut, so we continue.
    while length(E_local) > 0

        # Update which nodes are complete and which are incomplete
        C_local = [C_local; E_local]
        C_global = Local2Global[C_local]

        # Take these away from I_local
        I_local = setdiff(I_local,E_local)
        I_global = Local2Global[I_local]

        # To complete nodes in E, first add all the possible edges in the
        # current local graph, so that they match the global graph edges
        # (This is one of the most expensive parts of the expansion)
        L[E_local,E_local] = A[E_global,E_global]
        L[E_local,I_local] = A[E_global,I_global]
        L[I_local,E_local] = L[E_local,I_local]'

        # Now we must expand the local graph so that NEW neighbors of E
        # are added to L
        Lnew = Vector{Int64}()
        for v = E_global
            # This extracts the neighbor list of node v from the
            # rowval and colptr vectors of the adjacency matrix
            Neighbs_of_v = rp[ci[v]:ci[v+1]-1]
            for nv = Neighbs_of_v
                if ~inL[nv]
                    inL[nv] = true
                    push!(Lnew,nv)
                end
            end
        end
        numNew = length(Lnew)

        # We must add

        # Store local indices for new nodes added to L
        Lnew_local = collect((Lsize+1):(Lsize+numNew))

        # These are going to be "incomplete" nodes
        I_local = [I_local; Lnew_local]

        # Expand L by adding edges from the old local graph to Lnew.
        # Note that we don't include any edges between nodes in Lnew.
        P = A[Local2Global,Lnew]
        L = [L P;
            P' spzeros(numNew,numNew)]

        # Update the set of indices in L
        Local2Global = [Local2Global; Lnew]

        # excess stores the amount of "excess" flow after a flow computation.
        #
        # Extend the excess vector to accomodate the new size of L.
        # Since Lnew were not present in the last flow computation, they
        # have zero excess.
        excess = [excess; zeros(numNew)]

        # For the next local min-cut computation, we need to know which
        # nodes come with nonzero excess. These are "active" nodes.
        ExcessNodes = findall(x->x!=0,excess)

        # Update the capacity to the sink.
        tToL = [tToL; tWeights[Lnew]]
        # Now we construct a new local flow graph, and repeat

        Lf = [spzeros(1,1) sparse(tToL');
             sparse(tToL) L]

        Fold = F    # Old flow, saved as a warm start

        # Construct an initial flow F that includes the previous flow Fold
        # as a warm start. First, we allocate space for future
        # flow.
        # (This is one of the most expensive parts of the expansion)
        nLf = size(Lf,1)

        F = SparseMatrixCSC(nLf,nLf,Lf.colptr,Lf.rowval,zeros(length(Lf.rowval)))
        F[1:Lsize+1,1:Lsize+1] = Fold

        Lsize = size(L,1)

        # Compute min s-t cut for local flow graph and see if we need to expand
        S_local,F,excess = FlowSeed_Push_Relabel(Lf,F,ExcessNodes,excess)

        E_local = setdiff(S_local,C_local)     # the nodes that need completing
        E_global = Local2Global[E_local]       # their global indices

    end

    # return the global indices of the minimum cut set
    return Local2Global[S_local]

end

# A non-local version of the min-cut code that works by calling the same
# subroutine, but on the entire graph all at once
function NonLocalPushRelabel(A::SparseMatrixCSC,R::Vector{Int64},
    sWeights::Array{Float64},tWeights::Array{Float64})

        # Directly set up the flow matrix
        C = [spzeros(1,1) sparse(tWeights');
            sparse(tWeights) A]

        # Allocate space for the flow we will calculate
        F = SparseMatrixCSC(n+1,n+1,C.colptr,C.rowval,zeros(length(C.rowval)))

        # R is the set of nodes with excess, and the excess
        # will come from source-side edges that are immediately saturated
        S, F, excess = FlowSeed_Push_Relabel(C,F,R.+1,[0;sWeights])

        # The returned F is a preflow, not the maximum flow.
        # We are only interested in the cut.

        return S
end

# FlowSeed_Push_Relabel returns a preflow F and the min s-t cut set S for the
# flow graph C. It does not solve the maximum s-t flow problem.
#
# C = the capacity matrix for the flow problem.
#   Node 1 is the sink, and there is no explicit representation of a source,
#   the preflow immediately pushes all flow from the source to create an
#   excess on nodes in the graph.
#
# F = an initial flow. It can be initialize to zero.
#
# ExcessNodes = the set of nodes which at the start begin with some positive excess
#     These can be thought of as nodes that are adjacenct to the implicit source node
#     and the edges from the source are flooded. Or they may represent nodes that
#     have a nonzero excess from the initial flow F. The indices given here
#     should account for the fact that node 1 is already reserved for the sink.
#
# excess = the vector of excess values at the start of the algorithm. If F = 0,
#   this is the vector of edge capacities from the implicit source to the graph.
#   If F != 0, then it's the excess from a previous run of the algorithm
function FlowSeed_Push_Relabel(C::SparseMatrixCSC,
    F::SparseMatrixCSC,ExcessNodes::Array{Int64},excess::Array{Float64})

    # check excess node list
    # assert(countnz(excess) == length(ExcessNodes))

    # here, n includes only one terminal node, the sink
    n = size(C,1)

    height = zeros(Int64,n)      # label/height of each node
    inQ = zeros(Bool,n)          # list whether or not nodes are in the queue

    # Store adjacency list. There are ways to update this if calling
    # this function multiple times for growing local graphs, but it
    # does not appear to be a bottleneck to simply recompute frequently
    Neighbs,d = ConstructAdj(C,n)

    # We will maintain a queue of active nodes.
    Queue = Vector{Int64}()
    # An actual queue implementation is available in the DataStructures.jl
    # Julia package. The performane is nearly identical (and in some cases
    # slightly slower), thus to minimize dependency on outside packages, we
    # just use a Vector.

    # All nodes with nonzero excess are the first to be processed
    for v = ExcessNodes
        push!(Queue,v)
    end
    inQ[ExcessNodes] .= true

    # count the number of nodes that have been relabeled
    relabelings::Int64 = 0

    height = relabeling_bfs(C,F,1)     # compute initial distance from sink

    # In the code and comments, height = distance from sink = label of node

    # Continue until the queue no longer contains any active nodes.
    while length(Queue) > 0

        u = pop!(Queue)     # Select a new active node
        inQ[u] = false      # It's no longer in the queue

        if height[u] < n    # Check that the node is still active

            # discharge flow through node u
            relabelings += discharge_fs!(C,F,Queue,u,Neighbs[u],height,excess,n,d[u],inQ)

            # if u is still active, re-place it into the queue
            if excess[u] > 0 && height[u] < n
                prepend!(Queue,u)
                inQ[u] = true
            end

        end

        # Global relabeling heuristic for push-relabel algorithm.
        # This recomputes distances between nodes and the sink
        if relabelings == n
            relabelings = 0
            dist = relabeling_bfs(C,F,1)
            height = dist
        end

    end

    # Compute final distances from sink using BFS. Anything with distance
    # n will be the cut set.
    finalHeight = relabeling_bfs(C,F,1)
    S = Vector{Int64}()
    for i = 2:n
        if finalHeight[i] == n
            push!(S,i-1)
        end
    end

    excess[1] = 0.0     # ignore whatever excess there was at the sink.
    return S, F, excess

end

# FlowSeed version of discharge function.
#
# Discharege operation: pushes flow away from node u across admissible edges.
# If excess[u] > 0 but no admissible edges exist, we relabel u.
function discharge_fs!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    excess::Array{Float64},n::Int64,du::Int64,inQ::Array{Bool})

    vLocal::Int64 = 1
    hu = height[u]
    relabeled = 0
    while excess[u] > 0 && vLocal <= du
            v = uNeighbs[vLocal]
            if hu > height[v] && C[u,v] - F[u,v] > 0
                pushflow_fs!(C,F,Queue,u,v,excess,height,inQ,n)
                vLocal += 1
            else
                vLocal += 1
            end
    end

    if vLocal > du
        relabeled = 1
        relabel!(C,F,Queue,u,uNeighbs,height,du,n)
    end

    return relabeled
end

# Push flow from an active node u to a node v via an admissible edge (u,v)
function pushflow_fs!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int},u::Int64,v::Int64,excess::Array{Float64},height::Array{Int64},
    inQ::Array{Bool},n::Int64)

    send = min(excess[u], C[u,v] - F[u,v])
    F[u,v] += send
    F[v,u] -= send
    excess[u] -= send
    excess[v] += send

    # If v isn't in the queue, isn't the sink, is active, add it to the Queue
    if ~inQ[v] && v > 1 && height[v] < n
        prepend!(Queue,v)
        inQ[v] = true
    end
end
