include("material.jl")

Base.@kwdef mutable struct Particle
    x::Float64
    v::Float64
    m::Float64
    mat::Material
    lₚ::Float64
    Vₒ::Float64
    V::Float64
    σ::Float64
    ε::Float64
    F::Float64

    function Particle(x, m, mat, lₚ, v)
        Vₒ = 2*lₚ
        V = Vₒ
        σ = 0 
        ε = 0
        F = 1
        return new(x, v, m, mat, lₚ, Vₒ, V, σ, ε, F)
    end
end