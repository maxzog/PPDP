
abstract type Concentration end

function conc(ps::Vector{<:Particle}, ax::Mesh1D)
   cs=zeros(ax.ncells)
   for p in ps
      ind=ceil(Int64, (p.pos[2]-ax.xmin)/ax.dx)
      cs[ind]+=1
   end
   return cs
end

function conc(ps::Vector{<:Particle}, ax::Mesh3D)
   cs=zeros(ax.nCellsY)
   for p in ps
      ind=ceil(Int64, (p.pos[2]-ax.yMin)/ax.dy)
      cs[ind]+=1
   end
   return cs
end

#function conc(ps::Vector{Particle}, ax::Mesh3D, ::Type{Concentration})
#   cs=zeros(ax.ncells)
#   for p in ps
#      ind=ceil(Int64, (p.pos[2]-ax.xmin)/ax.dx)
#      cs[ind]+=1
#   end
#   return cs
#end
