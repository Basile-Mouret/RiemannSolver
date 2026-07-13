abstract type AbstractBC end

struct Outflow <: AbstractBC end

function apply_ghost(::Outflow, u_interior::AbstractVector{T}, x, t::T, normal::AbstractVector{T}) where {T<:Real}
    return u_interior
end

# 1D boundary conditions

abstract type AbstractBC1D <: AbstractBC end

function apply_ghost(bc::AbstractBC1D, ::AbstractVector{T}, x, ::T, normal::AbstractVector) where {T<:Real}
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet <: AbstractBC1D
    value::Function
end

function apply_ghost(bc::Dirichlet, ::AbstractVector{T}, x, t::T, normal::AbstractVector) where {T<:Real}
    return bc.value(t)
end

struct Reflecting <: AbstractBC1D end


# 2D boundary conditions

abstract type AbstractBC2D <: AbstractBC end

function apply_ghost(bc::AbstractBC2D, u_interior::AbstractVector{T}, x, t::T, normal::AbstractVector) where {T<:Real}
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet2D{F<:Function} <: AbstractBC2D
    value::F
end

function apply_ghost(bc::Dirichlet2D{F}, uR::SVector{N,T}, x, t::T, normal::AbstractVector) where {F<:Function, N, T<:Real}
    return bc.value(x, t)::SVector{N,T}
end

struct Reflecting2D <: AbstractBC2D end

