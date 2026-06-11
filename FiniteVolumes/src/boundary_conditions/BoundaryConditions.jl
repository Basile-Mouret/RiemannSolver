abstract type AbstractBC end

struct Outflow <: AbstractBC end

function apply_ghost(::Outflow, u_interior::Vector{Float64}, x, ::Float64)
    return copy(u_interior)
end


# 1D boundary conditions

abstract type AbstractBC1D <: AbstractBC end

function apply_ghost(bc::AbstractBC1D, u_interior::Vector{Float64}, x, t::Float64)::Vector{Float64}
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet <: AbstractBC1D
    value::Function
end

function apply_ghost(bc::Dirichlet, ::Vector{Float64}, x, t::Float64)
    return bc.value(t)
end

struct Reflecting <: AbstractBC1D end


# 2D boundary conditions

abstract type AbstractBC2D <: AbstractBC end

function apply_ghost(bc::AbstractBC2D, u_interior::Vector{Float64}, x, t::Float64)::Vector{Float64}
    error("apply_ghost not implemented for $(typeof(bc))")
end

struct Dirichlet2D <: AbstractBC2D
    value::Function
end

function apply_ghost(bc::Dirichlet2D, ::Vector{Float64}, x, t::Float64)
    return bc.value(x, t)
end

struct Reflecting2D <: AbstractBC2D end
