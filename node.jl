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

function N(node::Node, x::Float64, lₚ::Float64, L::Float64)
    
    Δx = x - node.x

    if (-L - lₚ < Δx) && (Δx <= -L + lₚ)
        return (L + lₚ + Δx) * (L + lₚ + Δx) / (4 * L * lₚ)	
    end

    if (-L + lₚ < Δx) && (Δx <= -lₚ)
        return 1 + Δx / L	
    end

    if (-lₚ < Δx) && (Δx <= lₚ)
        return 1 - (Δx * Δx + lₚ * lₚ) / (2 * L * lₚ)		
    end

    if (lₚ < Δx) && (Δx <= L - lₚ)
        return 1 - Δx / L			
    end

    if (L - lₚ < Δx) && (Δx <= L + lₚ)
        return (L + lₚ - Δx) * (L + lₚ - Δx) / (4 * L * lₚ)
    end
    
    return 0
end

function ∇N(node::Node, x::Float64, lₚ::Float64, L::Float64)

    Δx = x - node.x

    if (-L - lₚ < Δx) && (Δx <= -L + lₚ)
        return (L + lₚ + Δx) / (2 * L * lₚ)
    end

    if (-L + lₚ < Δx) && (Δx <= -lₚ)
        return 1 / L	
    end

    if (-lₚ < Δx) && (Δx <= lₚ)
        return -Δx / (L * lₚ)	
    end

    if (lₚ < Δx) && (Δx <= L - lₚ)
        return -1 / L			
    end

    if (L - lₚ < Δx) && (Δx <= L + lₚ)
        return -(L + lₚ - Δx) / (2 * L * lₚ)
    end

    return 0

end