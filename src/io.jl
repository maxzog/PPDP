
const ThreeVec{T} = SVector{3,T}
const ThreeMat{T} = SMatrix{3,3,T}

abstract type Particle end
abstract type Field end
abstract type Mesh end

mutable struct ChnlPartFull<:Particle
    # 64 bit integers
    id::Int64
    id0::Int64
    sub::Int64
    colWall::ThreeVec{Int64}
    # 64 bit floats
    d::Float64
    pos::ThreeVec{Float64}
    vel::ThreeVec{Float64}
    ang::ThreeVec{Float64}
    Acol::ThreeVec{Float64}
    Tcol::ThreeVec{Float64}
    dt::Float64
    delta::ThreeMat{Float64}
    omega::ThreeMat{Float64}
    coltime::Float64
    # 32 bit integers
    ind::ThreeVec{Int32}
    flag::Int32
    grp::Int32
end

mutable struct ChnlPart<:Particle
    # 64 bit integers
    id::Int64
    id0::Int64

    # 32 bit integers
    grp::Int32

    # 64 bit floats
    d::Float64
    coltime::Float64
    pos::ThreeVec{Float64}
    vel::ThreeVec{Float64}
    ang::ThreeVec{Float64}
end

mutable struct ChnlPartEns<:Particle
    # 64 bit integers
    id::Int64
    id0::Int64

    # 32 bit integers
    grp::Int32

    # 32 bit floats
    d::Float32
    coltime::Float32
    pos::ThreeVec{Float32}
    vel::ThreeVec{Float32}
    fld::ThreeVec{Float32}
end

function read_particles(fn::String, ::Type{ChnlPartFull})
    # Open the binary file
    io = open(fn, "r")

    # Read the number of particles and particle size
    header = Vector{Int32}(undef, 2)
    read!(io, header)
    np, size = header[1], header[2]

    # Initialize the particle array
    ps = Vector{ChnlPartFull}(undef, np)

    # Define reusable buffers
    int_buf = Vector{Int64}(undef, 1)
    intv_buf = Vector{Int64}(undef, 3)
    float_buf = Vector{Float64}(undef, 1)
    floatv_buf = Vector{Float64}(undef, 3)
    floata_buf = Array{Float64}(undef, (3, 3))
    int32_buf = Vector{Int32}(undef, 1)
    int32v_buf = Vector{Int32}(undef, 3)

    # Loop to read all particles
    for i = 1:np
        # Create an empty particle structure
        p = ChnlPartFull(
            0,
            0,
            0,
            zeros(Int64, 3),
            0.0,
            zeros(Float64, 3),
            zeros(Float64, 3),
            zeros(Float64, 3),
            zeros(Float64, 3),
            zeros(Float64, 3),
            0.0,
            zeros(Float64, (3, 3)),
            zeros(Float64, (3, 3)),
            0.0,
            zeros(Int32, 3),
            Int32(0),
            Int32(-1)
        )

        # Read particle fields directly into the structure
        read!(io, int_buf)
        p.id = int_buf[1]
        read!(io, int_buf)
        p.id0 = int_buf[1]
        read!(io, int_buf)
        p.sub = int_buf[1]
        read!(io, intv_buf)
        p.colWall = intv_buf

        read!(io, float_buf)
        p.d = float_buf[1]
        read!(io, floatv_buf)
        p.pos = floatv_buf
        read!(io, floatv_buf)
        p.vel = floatv_buf
        read!(io, floatv_buf)
        p.ang = floatv_buf
        read!(io, floatv_buf)
        p.Acol = floatv_buf
        read!(io, floatv_buf)
        p.Tcol = floatv_buf
        read!(io, float_buf)
        p.dt = float_buf[1]
        read!(io, floata_buf)
        p.delta = floata_buf
        read!(io, floata_buf)
        p.omega = floata_buf
        read!(io, float_buf)
        p.coltime = float_buf[1]

        read!(io, int32v_buf)
        p.ind = int32v_buf
        read!(io, int32_buf)
        p.flag = int32_buf[1]

        # Store the particle in the array
        ps[i] = p
    end

    # Close the file and return the particle array
    close(io)
    assign_groups!(ps)
    return ps
end

function read_particles(fn::String, ::Type{ChnlPart})
    # Open the binary file
    io = open(fn, "r")

    # Read the number of particles and particle size
    header = Vector{Int32}(undef, 2)
    read!(io, header)
    np, size = header[1], header[2]

    # Initialize the particle array
    ps = Vector{ChnlPart}(undef, np)

    # Define reusable buffers
    int_buf = Vector{Int64}(undef, 1)
    intv_buf = Vector{Int64}(undef, 3)
    float_buf = Vector{Float64}(undef, 1)
    floatv_buf = Vector{Float64}(undef, 3)
    floata_buf = Array{Float64}(undef, (3, 3))
    int32_buf = Vector{Int32}(undef, 1)
    int32v_buf = Vector{Int32}(undef, 3)

    # Loop to read all particles
    for i = 1:np
        # Create an empty particle structure
        p = ChnlPart(
            0,
            0,
            Int32(0),
            0.0,
            0.0,
            zeros(Float64, 3),
            zeros(Float64, 3),
            zeros(Float64, 3),
        )

        # Read particle fields directly into the structure
        read!(io, int_buf)
        p.id = int_buf[1]
        read!(io, int_buf)
        p.id0 = int_buf[1]
        read!(io, int_buf)
        read!(io, intv_buf)

        read!(io, float_buf)
        p.d = float_buf[1]
        read!(io, floatv_buf)
        p.pos = floatv_buf
        read!(io, floatv_buf)
        p.vel = floatv_buf
        read!(io, floatv_buf)
        p.ang = floatv_buf
        read!(io, floatv_buf)
        read!(io, floatv_buf)
        read!(io, float_buf)
        read!(io, floata_buf)
        read!(io, floata_buf)
        read!(io, float_buf)
        p.coltime = float_buf[1]

        read!(io, int32v_buf)
        read!(io, int32_buf)

        # Store the particle in the array
        ps[i] = p
    end

    # Close the file and return the particle array
    close(io)

    assign_groups!(ps)

    return ps
end

function read_particles(dir::String, step::Int64, ::Type{ChnlPartEns})
    suf = "."*"0"^(6-length(string(step)))*string(step)
    np = get_npart(dir*"/particle"*suf)
    X  = read_pos(dir*"/particle"*suf, np)
    U  = read_arr(dir*"/fld"*suf, np)
    V  = read_arr(dir*"/vel"*suf, np)
    ids= read_vec(dir*"/id"*suf, np)
    id0s= read_vec(dir*"/id0"*suf, np)
    ds = read_vec(dir*"/dp"*suf, np)
    cs = read_vec(dir*"/coltime"*suf, np)
    ps = Vector{ChnlPartEns}(undef, np)
    for i in 1:np
       ps[i] = ChnlPartEns(trunc(Int64, ids[i]), trunc(Int64, id0s[i]), -1, ds[i], cs[i], X[:,i], V[:,i], U[:,i])
    end
    perm = sortperm([p.id0 for p in ps])
    assign_groups!(ps)
    return ps[perm]
end

# Assign group number to particle based on diameter
function assign_groups!(ps::Vector{<:Particle}; tol = 1e-10)
    # Get diameter values
    dps = unique([p.d for p in ps])
    sort!(dps)
    # Get number of groups
    ngrp = lastindex(dps)
    # Assign
    for p in ps
        for i = 1:ngrp
            if abs(p.d - dps[i]) < tol
                p.grp = Int32(i)
            end
        end
    end
end

#######################################
###
### Ensight file reading (particles)
###
#######################################

function get_npart(fn::String)::Int32
    """
    Reads the ensight particle.****** files to get the number of particles
    There's probably a more elegant way to do this but hey this works 
    """
    if split(split(fn, '.')[end-1], '/')[end] != "particle"
        print(fn)
        error("[get_npart] Not the correct file for reading npart")
        return 
    end
    io = open(fn)
    skip(io, 240)
    n = Vector{Int32}(undef, 1)
    read!(io, n); close(io)
    return n[1]
end

function read_pos(fn::String, n::Int32)::Array{Float32}
    """
    Reads ensight particle geometry files
    outputs (3, npart) array of position data
    """
    arr = Array{Float32}(undef, (3, n))
    io = open(fn, "r")
    skip(io, 244+4*n)
    read!(io, arr)
    close(io); return arr
end

function read_vec(fn::String, n::Int32)::Vector{Float32}
    """
    Reads ensight particle data files and outputs a vector of particle scalar data
    The header size can vary (as far as I can tell) so read the particle.****** file first using get_npart()
    to determine how much to skip in the data file
    """
    vec = Vector{Float32}(undef, n)
    io = open(fn)
    skip(io, 80)
    read!(io, vec); close(io)
    return vec
end

function read_arr(fn::String, n::Int32)::Array{Float32}
    """
    Reads ensight particle data files and outputs a (3, npart) array of the vector data
    The header size can vary (as far as I can tell) so read the particle.****** file first using get_npart()
    to determine how much to skip in the data file
    """
    arr = Array{Float32}(undef, (3, n))
    io = open(fn)
    skip(io, 80)
    read!(io, arr); close(io)
    return arr
end

#######################################
###
### Mesh definitions and construction
###
#######################################

mutable struct Mesh1D{T<:AbstractFloat} <: Mesh
    nCells::Int64
    nEdges::Int64

    xMin::T
    xMax::T
    Lx  ::T
    dx  ::T

    centersX::Vector{T}
    edgesX  ::Vector{T}
end

mutable struct Mesh2D <: Mesh
    dim::Int64

    nCellsX::Int64
    nCellsY::Int64

    nEdgesX::Int64
    nEdgesY::Int64

    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64

    Lx::Float64
    Ly::Float64

    dx::Float64
    dy::Float64

    centersX::Vector{Float64}
    edgesX::Vector{Float64}
    centersY::Vector{Float64}
    edgesY::Vector{Float64}
end

mutable struct Mesh3D <: Mesh
    dim::Int64

    nCellsX::Int64
    nCellsY::Int64
    nCellsZ::Int64

    nEdgesX::Int64
    nEdgesY::Int64
    nEdgesZ::Int64

    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64
    zMin::Float64
    zMax::Float64

    Lx::Float64
    Ly::Float64
    Lz::Float64

    dx::Float64
    dy::Float64
    dz::Float64

    centersX::Vector{Float64}
    edgesX::Vector{Float64}
    centersY::Vector{Float64}
    edgesY::Vector{Float64}
    centersZ::Vector{Float64}
    edgesZ::Vector{Float64}
end

mutable struct NonUniformMesh3D <: Mesh
    dim::Int64

    nCellsX::Int64
    nCellsY::Int64
    nCellsZ::Int64

    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64
    zMin::Float64
    zMax::Float64

    Lx::Float64
    Ly::Float64
    Lz::Float64

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    centersX::Vector{Float64}
    centersY::Vector{Float64}
    centersZ::Vector{Float64}
end

function mesh(xmin::T, xmax::T, ncells::Int64)::Mesh1D{T} where {T <: AbstractFloat}
    @assert(xmin < xmax)
    L = xmax - xmin
    d = L / ncells
    xvcell = LinRange(xmin + d / 2, xmax - d / 2, ncells)
    xvedge = LinRange(xmin, xmax, ncells + 1)

    obj = Mesh1D{T}(ncells, ncells + 1, xmin, xmax, L, d, xvcell, xvedge)
    return obj
end

function mesh(minCorner::SVector{2,Float64}, maxCorner::SVector{2,Float64}, ncells::SVector{2,Int64})::Mesh2D
    for i in 1:2
       @assert(minCorner[i] < maxCorner[i])
    end
    Lv = maxCorner .- minCorner
    dv = Lv ./ ncells
    xcv = LinRange(minCorner[1] + dv[1] / 2, maxCorner[1] - dv[1] / 2, ncells[1])
    ycv = LinRange(minCorner[2] + dv[2] / 2, maxCorner[2] - dv[2] / 2, ncells[2])
    xev = LinRange(minCorner[1], maxCorner[1], ncells[1]+1)
    yev = LinRange(minCorner[2], maxCorner[2], ncells[2]+1)
    obj = Mesh2D(3, 
                 ncells[1], ncells[2],
                 ncells[1]+1, ncells[2]+1,
                 minCorner[1], maxCorner[1], 
                 minCorner[2], maxCorner[2], 
                 Lv[1], Lv[2],
                 dv[1], dv[2],
                 xcv, xev,
                 ycv, yev)
    return obj
end

#function mesh(minCorner::Point{2,Float64}, maxCorner::Point{2,Float64}, ncells::Vector{Int64})::Mesh2D
#    for i in 1:2
#       @assert(minCorner[i] < maxCorner[i])
#    end
#    Lv = maxCorner .- minCorner
#    dv = Lv ./ ncells
#    xcv = LinRange(minCorner[1] + dv[1] / 2, maxCorner[1] - dv[1] / 2, ncells[1])
#    ycv = LinRange(minCorner[2] + dv[2] / 2, maxCorner[2] - dv[2] / 2, ncells[2])
#    xev = LinRange(minCorner[1], maxCorner[1], ncells[1]+1)
#    yev = LinRange(minCorner[2], maxCorner[2], ncells[2]+1)
#    obj = Mesh2D(3, 
#                 ncells[1], ncells[2],
#                 ncells[1]+1, ncells[2]+1,
#                 minCorner[1], maxCorner[1], 
#                 minCorner[2], maxCorner[2], 
#                 Lv[1], Lv[2],
#                 dv[1], dv[2],
#                 xcv, xev,
#                 ycv, yev)
#    return obj
#end


function mesh(minCorner::SVector{3,Float64}, maxCorner::SVector{3,Float64}, ncells::SVector{3,Int64})::Mesh3D
    for i in 1:3
       @assert(minCorner[i] < maxCorner[i])
    end
    Lv = maxCorner .- minCorner
    dv = Lv ./ ncells
    xcv = LinRange(minCorner[1] + dv[1] / 2, maxCorner[1] - dv[1] / 2, ncells[1])
    ycv = LinRange(minCorner[2] + dv[2] / 2, maxCorner[2] - dv[2] / 2, ncells[2])
    zcv = LinRange(minCorner[3] + dv[3] / 2, maxCorner[3] - dv[3] / 2, ncells[3])
    xev = LinRange(minCorner[1], maxCorner[1], ncells[1]+1)
    yev = LinRange(minCorner[2], maxCorner[2], ncells[2]+1)
    zev = LinRange(minCorner[3], maxCorner[3], ncells[3]+1)
    obj = Mesh3D(3, 
                 ncells[1], ncells[2], ncells[3], 
                 ncells[1]+1, ncells[2]+1, ncells[3]+1,
                 minCorner[1], maxCorner[1], 
                 minCorner[2], maxCorner[2], 
                 minCorner[3], maxCorner[3],
                 Lv[1], Lv[2], Lv[3],
                 dv[1], dv[2], dv[3], 
                 xcv, xev,
                 ycv, yev,
                 zcv, zev)
    return obj
end

function mesh(minCorner::SVector{3,Float64}, maxCorner::SVector{3,Float64}, ncells::SVector{3,Int64},
              XV::Vector{Float64}, YV::Vector{Float64}, ZV::Vector{Float64},
              dXV::Vector{Float64}, dYV::Vector{Float64}, dZV::Vector{Float64})::NonUniformMesh3D
    for i in 1:3
       @assert(minCorner[i] < maxCorner[i])
    end
    Lv = maxCorner .- minCorner
    obj = NonUniformMesh3D(3, 
                           ncells[1], ncells[2], ncells[3], 
                           minCorner[1], maxCorner[1], 
                           minCorner[2], maxCorner[2], 
                           minCorner[3], maxCorner[3],
                           Lv[1], Lv[2], Lv[3],
                           dXV, dYV, dZV, 
                           XV,
                           YV,
                           ZV)
    return obj
end


#######################################
###
### Field/grid definitions and construction 
###
#######################################

mutable struct StretchedChannelGrid
   mesh :: NonUniformMesh3D

   U :: Array{Float32}
   V :: Array{Float32}
   W :: Array{Float32}
end

mutable struct ChannelGrid
   mesh :: Mesh3D

   U :: Array{Float32}
   V :: Array{Float32}
   W :: Array{Float32}
end

function read_field(fn::String, step::Int64, ::Type{StretchedChannelGrid})::StretchedChannelGrid
   # Read mesh info
   Ls,nx,ny,nz,xv,yv,zv=read_mesh(fn) 
   # Read velocity field
   fn_vel=fn*"/velocity/velocity."*string(step,pad=6)
   U,V,W = read_vel(fn_vel,nx,ny,nz) 
   
   # Convert to doubles
   nxDP=Int64(nx)
   nyDP=Int64(ny)
   nzDP=Int64(nz)
   LsDP=Float64.(Ls)

   minCorner=SVector{3,Float64}([LsDP[1], LsDP[3], LsDP[5]])
   maxCorner=SVector{3,Float64}([LsDP[2], LsDP[4], LsDP[6]])
   nCells=SVector{3,Int64}([nxDP,nyDP,nzDP])

   xvd=0.5*Float64.(xv[1:end-1] + xv[2:end])
   yvd=0.5*Float64.(yv[1:end-1] + yv[2:end])
   zvd=0.5*Float64.(zv[1:end-1] + zv[2:end])

   dx=Float64.(xv[2:end]-xv[1:end-1])
   dy=Float64.(yv[2:end]-yv[1:end-1])
   dz=Float64.(zv[2:end]-zv[1:end-1])

   obj = StretchedChannelGrid(mesh(minCorner,maxCorner,nCells,xvd,yvd,zvd,dx,dy,dz),
                     U, V, W)
   return obj
end

function read_field(fn::String, step::Int64, ::Type{ChannelGrid})::ChannelGrid
   # Read mesh info
   Ls,nx,ny,nz,xv,yv,zv=read_mesh(fn) 
   # Read velocity field
   fn_vel=fn*"/velocity/velocity."*string(step,pad=6)
   U,V,W = read_vel(fn_vel,nx,ny,nz) 
   
   # Convert to doubles
   nxDP=Int64(nx)
   nyDP=Int64(ny)
   nzDP=Int64(nz)
   LsDP=Float64.(Ls)

   minCorner=SVector{3,Float64}([LsDP[1], LsDP[3], LsDP[5]])
   maxCorner=SVector{3,Float64}([LsDP[2], LsDP[4], LsDP[6]])
   nCells=SVector{3,Int64}([nxDP,nyDP,nzDP])

   xvd=Float64.(xv)
   yvd=Float64.(yv)
   zvd=Float64.(zv)

   obj = StretchedChannelGrid(mesh(minCorner,maxCorner,nCells,xvd,yvd,zvd),
                     U, V, W)
   return obj
end

function read_vel(fn::String, nx::Int32, ny::Int32, nz::Int32)::Tuple{Array{Float32}, Array{Float32}, Array{Float32}}
    """
    Reads velocity field
    Returns three arrays (U, V, W)
    """
    U = Array{Float32}(undef, (nx,ny,nz))
    V = Array{Float32}(undef, (nx,ny,nz))
    W = Array{Float32}(undef, (nx,ny,nz))
    io = open(fn)
    skip(io, 244)
    read!(io, U)
    read!(io, V)
    read!(io, W)
    return U, V, W
end

function read_mesh(fn::String)
   ## Read the mesh
   io = open(fn*"/geometry", "r")
   # skip ensight crap
   cbuff = Vector{Char}(undef, 80)
   for _ in 1:6
      for i in eachindex(cbuff)
         cbuff[i] = read(io, Char)
      end
   end
   lengths=Vector{Float32}(undef, 6)
   # read extents
   read!(io, lengths)
   # skip more ensight crap
   for i in eachindex(cbuff)
      cbuff[i] = read(io, Char)
   end
   read!(io, Vector{Int32}(undef, 1))
   for _ in 1:2 
      for i in eachindex(cbuff)
         cbuff[i] = read(io, Char)
      end
   end

   # read mesh
   ncells = Vector{Int32}(undef, 3)
   read!(io, ncells)
   xv = Vector{Float32}(undef, ncells[1])
   yv = Vector{Float32}(undef, ncells[2])
   zv = Vector{Float32}(undef, ncells[3])
   read!(io, xv)
   read!(io, yv)
   read!(io, zv)
   close(io)
   nx = Int32(ncells[1]-1)
   ny = Int32(ncells[2]-1)
   nz = Int32(ncells[3]-1)
   return lengths,nx, ny, nz, xv, yv, zv 
end

function read_P(fn::String, nx::Int32, ny::Int32, nz::Int32)::Array{Float32}
    """
    Reads pressure... what else do you want from me?
    """
    arr = Array{Float32}(undef, (nx,ny,nz))
    io = open(fn)
    skip(io, 244)
    read!(io, arr); close(io)
    return arr
end
