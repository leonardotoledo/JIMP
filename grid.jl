include("element.jl")
include("particle.jl")

Base.@kwdef struct Grid
    nodes::Vector{Node}
    elements::Vector{Element}
    lₑ::Float64
    xₘᵢₙ::Float64
    xₘₐₓ ::Float64
end

function reset(grid::Grid)
    for e in grid.elements
        reset(e)
    end
end

function lockNodeByPosition(grid::Grid, x::Float64)
    for e in grid.elements
        for n in e.nodes
            if n.x == x
                n.isLocked = true
                return
            end
        end
    end
end

function mapParticleToElement(grid::Grid, particle::Particle)
    id = convert(Int64, ceil((particle.x - grid.xₘᵢₙ)/grid.lₑ))

    @assert (id >= 1 && id <= length(grid.elements)) "A particle has gone outside the grid domain limits."
    return grid.elements[id]
end