

function fftcorrgrid(U::Array{Float32}, comp::Int64)
    sumd = [1,2,3]
    deleteat!(sumd, findall(x->x==comp,sumd))
    uh = fft(U, [comp])
    uhs = conj(uh)
    u = sum(real.(ifft(uh.*uhs, [comp])), dims = sumd) / lastindex(U[:,1,1])^3
    return dropdims(u, dims = Tuple(sumd))[1:trunc(Int64, lastindex(u)/2)]
end

function R11(hit::HITGrid)
    corr = (fftcorrgrid(hit.U, 1) + fftcorrgrid(hit.V, 2) + fftcorrgrid(hit.W, 3))./3
    return corr
end

function R22(hit::HITGrid)
    corr = (fftcorrgrid(hit.U, 2) + fftcorrgrid(hit.V, 3) + fftcorrgrid(hit.W, 1))./3
    return corr
end

function getSpectrum(hit::HITGrid)
    Nx=hit.mesh.nCellsX
    n = trunc(Int, Nx/2)
    kv = fftshift(-n:n-1)
    dk = 2*π/hit.mesh.Lx
    Uk, Vk, Wk = spectral(hit)
    spec = zeros(Float32, ceil(Int64, norm([maximum(abs.(kv))+.5 for _ in 1:3])))
    counts = zeros(Int64, length(spec))
    @inbounds for k∈1:hit.mesh.nCellsZ, j∈1:hit.mesh.nCellsY, i∈1:hit.mesh.nCellsX
        kk = norm([kv[i], kv[j], kv[k]])
        ik = floor(Int64, kk/dk+0.5)+1
        if ik > 1e-6
            spec[ik] += 0.5*(real(Uk[i,j,k]*conj(Uk[i,j,k])) + real(Vk[i,j,k]*conj(Vk[i,j,k])) + real(Wk[i,j,k]*conj(Wk[i,j,k])))
        end
    end
    return spec/Nx^6 
end

function spectral(hit::HITGrid)
   p = plan_fft(hit.U)
   Uk = p * hit.U
   Vk = p * hit.V
   Wk = p * hit.W
   return Uk, Vk, Wk
end
