include("grid.jl")

Base.@kwdef mutable struct Domain
    grid::Grid
    particles::Vector{Particle}
    materials::Vector{Material}
    b::Float64
    Domain(grid) = new(grid, [], [], 0.0)
end

function addMaterial(domain::Domain, material::Material)
    push!(domain.materials, material)
end

function addParticle(domain::Domain, particle::Particle)
    push!(domain.particles, particle)
end

function genGrid(L, xₘᵢₙ, xₘₐₓ)

    nₑ = convert(Int64, ceil(xₘₐₓ-xₘᵢₙ)/L)
    nₙ = nₑ + 1

    nodes::Vector{Node} = []
    elements::Vector{Element} = []

    for i=1:nₙ
        push!(nodes, Node(xₘᵢₙ+(i-1)*L))
    end

    for i=1:nₑ
        push!(elements, Element([nodes[i], nodes[i+1]]))
    end

    for i=1:nₑ

        v₁ = i-1
        v₂ = i+1

        if v₁ >= 1
            elements[i].left = elements[v₁]
        end

        if v₂ <= nₑ
            elements[i].right = elements[v₂]
        end
    end

    return Grid(nodes, elements, L, xₘᵢₙ, xₘₐₓ)
end

function genParticles(domain::Domain, mat::Material, ppe::Int64, xₘᵢₙ::Float64, xₘₐₓ::Float64)

    dₓₚ::Float64 = domain.grid.L/ppe
    m::Float64 = dₓₚ*mat.ρ
    nₚ::Int64 = convert(Int64, ceil(xₘₐₓ-xₘᵢₙ)/dₓₚ)

    for i=1:nₚ
        addParticle(domain, Particle(xₘᵢₙ + (i-0.5)*dₓₚ, m, mat, 0.5*dₓₚ, 0))
    end

    addMaterial(domain, mat)
end