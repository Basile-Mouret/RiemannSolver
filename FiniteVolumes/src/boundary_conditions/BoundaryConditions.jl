abstract type AbstractBC end

struct Outflow <: AbstractBC end

function apply_ghost(::Outflow, u_interior::AbstractVector{Float64}, x, ::Float64)
    return u_interior
end

# 1D boundary conditions

abstract type AbstractBC1D <: AbstractBC end

function apply_ghost(bc::AbstractBC1D, ::AbstractVector{Float64}, x, ::Float64)
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet <: AbstractBC1D
    value::Function
end

function apply_ghost(bc::Dirichlet, ::AbstractVector{Float64}, x, t::Float64)
    return bc.value(t)
end

struct Reflecting <: AbstractBC1D end


# 2D boundary conditions

abstract type AbstractBC2D <: AbstractBC end

function apply_ghost(bc::AbstractBC2D, u_interior::AbstractVector{Float64}, x, t::Float64)
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet2D{F<:Function} <: AbstractBC2D
    value::F
end

function apply_ghost(bc::Dirichlet2D{F}, uR::SVector{N,T}, x, t::Float64) where {F<:Function, N, T}
    return bc.value(x, t)::SVector{N,T}
end

struct Reflecting2D <: AbstractBC2D end
