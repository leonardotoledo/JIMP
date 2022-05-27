Base.@kwdef mutable struct Node
    x::Float64
    m::Float64
    f::Float64
    p::Float64
    isLocked::Bool
    Node(x) = new(x, 0, 0, 0, false)
end

function reset(node::Node)
    node.m = 0
    node.f = 0
    node.p = 0
end

function N(node::Node, x::Float64, lₚ::Float64, lₑ::Float64)
    
    Δx = x - node.x

    if (-lₑ - lₚ < Δx) && (Δx <= -lₑ + lₚ)
        return (lₑ + lₚ + Δx) * (lₑ + lₚ + Δx) / (4 * lₑ * lₚ)	
    end

    if (-lₑ + lₚ < Δx) && (Δx <= -lₚ)
        return 1 + Δx / lₑ	
    end

    if (-lₚ < Δx) && (Δx <= lₚ)
        return 1 - (Δx * Δx + lₚ * lₚ) / (2 * lₑ * lₚ)		
    end

    if (lₚ < Δx) && (Δx <= lₑ - lₚ)
        return 1 - Δx / lₑ			
    end

    if (lₑ - lₚ < Δx) && (Δx <= lₑ + lₚ)
        return (lₑ + lₚ - Δx) * (lₑ + lₚ - Δx) / (4 * lₑ * lₚ)
    end
    
    return 0
end

function ∇N(node::Node, x::Float64, lₚ::Float64, lₑ::Float64)

    Δx = x - node.x

    if (-lₑ - lₚ < Δx) && (Δx <= -lₑ + lₚ)
        return (lₑ + lₚ + Δx) / (2 * lₑ * lₚ)
    end

    if (-lₑ + lₚ < Δx) && (Δx <= -lₚ)
        return 1 / lₑ	
    end

    if (-lₚ < Δx) && (Δx <= lₚ)
        return -Δx / (lₑ * lₚ)	
    end

    if (lₚ < Δx) && (Δx <= lₑ - lₚ)
        return -1 / lₑ			
    end

    if (lₑ - lₚ < Δx) && (Δx <= lₑ + lₚ)
        return -(lₑ + lₚ - Δx) / (2 * lₑ * lₚ)
    end

    return 0

end