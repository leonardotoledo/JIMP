include("node.jl")

Base.@kwdef mutable struct Element
    nodes::Vector{Node}
    left::Element
    right::Element
    Element(nodes) = new(nodes)
end

function reset(element::Element)
    for n in element.nodes
        reset(n)
    end
end