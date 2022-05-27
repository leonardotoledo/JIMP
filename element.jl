include("node.jl")

Base.@kwdef mutable struct Element
    nodes::Vector{Node} # Element - nodes
    left::Element       # Element - left neighbor
    right::Element      # Element - right neighbor
    Element(nodes) = new(nodes)
end

function reset(element::Element)
    for n in element.nodes
        reset(n)
    end
end