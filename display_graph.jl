# Plot a graph. Code taken (and slightly edited) from Huda Nassar's Julia tutorial:
# https://github.com/nassarhuda/JuliaTutorials/blob/master/plotting.ipynb
using LinearAlgebra
using Plots

function display_graph(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},grayscale::Float64)
  f = plot(leg=false, axis = false,grid = false)
  ei,ej,w = findnz(triu(A))
  lx = [xy[ei,1]';xy[ej,1]';NaN*ones(1,length(ei))]
  ly = [xy[ei,2]';xy[ej,2]';NaN*ones(1,length(ei))]
  for i = 1:length(w)
      plot!(f,lx[:,i],ly[:,i],color = RGB(grayscale,grayscale,grayscale), linewidth = 1)
  end
  scatter!(f,xy[:,1],xy[:,2],color = RGB(grayscale,grayscale,grayscale),markerstrokecolor =  RGB(grayscale,grayscale,grayscale))
  return f
end
