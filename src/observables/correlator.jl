"""
    propagator(H::MPOHamiltonian, bra::FiniteMPS, ket::FiniteMPS; rev::Bool=false, dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50))
    propagator(H::MPOHamiltonian, bras::Vector{<:FiniteMPS}, ket::FiniteMPS; rev::Bool=false, dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50))
"""
function propagator(H::MPOHamiltonian, bra::FiniteMPS, ket::FiniteMPS; rev::Bool=false, dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50))
    times = collect(0:dt:ft)
    propagators = zeros(ComplexF64, length(times))
    propagators[1] = dot(bra, ket)
    envs = environments(ket, H)
    for (i, t) in enumerate(times[2:end])
        alg = t > n * dt ? TDVP() : TDVP2(; trscheme=trscheme)
        ket, envs = timestep(ket, H, 0, dt, alg, envs)
        propagators[i+1] = dot(bra, ket)
    end
    rev ? propagators = conj.(propagators) : propagators = propagators
    return propagators
end

function propagator(H::MPOHamiltonian, bras::Vector{<:FiniteMPS}, ket::FiniteMPS; rev::Bool=false, dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50))
    times = collect(0:dt:ft)
    propagators = zeros(ComplexF64, length(bras), length(times))
    propagators[:,1] = [dot(bras[i], ket) for i in 1:length(bras)]
    envs = environments(ket, H)
    for (i, t) in enumerate(times[2:end])
        alg = t > n * dt ? TDVP() : TDVP2(; trscheme=trscheme)
        ket, envs = timestep(ket, H, 0, dt, alg, envs)
        for j in eachindex(bras)
            propagators[j,i+1] = dot(bras[j], ket)
        end
    end
    rev ? propagators = conj.(propagators) : propagators = propagators
    return propagators
end

struct RetardedGF{K} end
RetardedGF(::Type{RetardedGF{:f}}) = 1
RetardedGF(::Type{RetardedGF{:b}}) = -1


function dcorrelator(::Type{R}, H::MPOHamiltonian, gsenergy::Number, mps::Vector{<:FiniteMPS}; dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50)) where R<:RetardedGF
    t, half = collect(0:dt:ft), length(mps)÷2
    gf = SharedArray{ComplexF64, 3}(length(mps), half, length(0:dt:ft))
    @sync @distributed for i in 1:length(mps)
        if i <= half
            gf[i,:,:] = propagator(H, mps[1:half], mps[i]; rev=false, dt=dt, ft=ft, n=n, trscheme=trscheme) 
        else
            gf[i,:,:] = propagator(H, mps[(half+1):end], mps[i]; rev=true, dt=dt, ft=ft, n=n, trscheme=trscheme)
        end
    end
    for i in eachindex(t)
        factor₁, factor₂ = -im*exp(im*gsenergy*t[i]), -im*exp(-im*gsenergy*t[i])
        gf[1:half,:,i] = factor₁*gf[1:half,:,i]
        gf[(half+1):end,:,i] = factor₂*gf[(half+1):end,:,i]
    end
    return gf[1:half,:,:] + RetardedGF(R)*gf[(half+1):end,:,:]
end

struct GreaterLessGF end

function dcorrelator(::Type{GreaterLessGF}, H::MPOHamiltonian, gsenergy::Number, mps::Vector{<:FiniteMPS}; whichs=:both, dt::Number=0.05, ft::Number=5.0, n::Integer=3, trscheme=truncdim(50))
    t, half = collect(0:dt:ft), length(mps)÷2
    if whichs == :both
        gf = SharedArray{ComplexF64, 3}(length(mps), half, length(0:dt:ft))
        @sync @distributed for i in 1:length(mps)
            if i <= half
                gf[i,:,:] = propagator(H, mps[1:half], mps[i]; rev=false, dt=dt, ft=ft, n=n, trscheme=trscheme) 
            else
                gf[i,:,:] =  propagator(H, mps[(half+1):end], mps[i]; rev=true, dt=dt, ft=ft, n=n, trscheme=trscheme)
            end
        end
        for i in eachindex(t)
            factor₁, factor₂ = -im*exp(im*gsenergy*t[i]), -im*exp(-im*gsenergy*t[i])
            gf[1:half,:,i] = factor₁*gf[1:half,:,i]
            gf[(half+1):end,:,i] = factor₂*gf[(half+1):end,:,i]
        end
        return [gf[1:half,:,:], gf[(half+1):end,:,:]]
    elseif whichs == :greater
        gf = SharedArray{ComplexF64, 3}(length(mps), length(mps), length(0:dt:ft))
        @sync @distributed for i in 1:length(mps)
            gf[i,:,:] = propagator(H, mps[1:end], mps[i]; rev=false, dt=dt, ft=ft, n=n, trscheme=trscheme) 
        end
        for i in eachindex(t)
            gf[:,:,i] = -im*exp(im*gsenergy*t[i])*gf[:,:,i]
        end
        return gf
    elseif whichs == :less
        gf = SharedArray{ComplexF64, 3}(length(mps), length(mps), length(0:dt:ft))
        @sync @distributed for i in 1:length(mps)
            gf[i,:,:] =  propagator(H, mps[1:end], mps[i]; rev=true, dt=dt, ft=ft, n=n, trscheme=trscheme)
        end
        for i in eachindex(t)
            gf[:,:,i] = -im*exp(-im*gsenergy*t[i])*gf[:,:,i]
        end
        return gf
    end
end



struct GorkovGF end
struct MatsubaraGF end