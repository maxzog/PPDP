
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

function read_particles(fn::String, ::Type{ChnlPartFull})
    # Open the binary file
    io = open(fn, "r")

    # Read the number of particles and particle size
    header = Vector{Int32}(undef, 2)
    read!(io, header)
    np, size = header[1], header[2]
    println("Npart :: ", np)
    println("Size  :: ", size)

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
    println("Npart :: ", np)
    println("Size  :: ", size)

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

    assign_Tgroups!(ps)

    return ps
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
### Mesh definitions and construction
###
#######################################

mutable struct Mesh1D{T<:AbstractFloat} <: Mesh
    nCells::Int64
    nEdges::Int64

    xMin::T
    xMax::T
    Lx::T
    dx::T

    centersX::Vector{T}
    edgesX::Vector{T}
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

function mesh(xmin::Float64, xmax::Float64, ncells::Int64)::Mesh1D{Float64}
    @assert(xmin < xmax)
    L = xmax - xmin
    d = L / ncells
    xvcell = LinRange(xmin + d / 2, xmax - d / 2, ncells)
    xvedge = LinRange(xmin, xmax, ncells + 1)

    obj = Mesh1D{Float64}(ncells, ncells + 1, xmin, xmax, L, d, xvcell, xvedge)
    return obj
end

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


#######################################
###
### Field/grid definitions and construction 
###
#######################################

mutable struct ChannelGrid
   mesh :: Mesh3D

   U :: Array{Float32}
   V :: Array{Float32}
   W :: Array{Float32}
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

   obj = ChannelGrid(mesh(minCorner,maxCorner,nCells), 
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




