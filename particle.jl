include("material.jl")

Base.@kwdef mutable struct Particle
    x::Float64      # Particle - position
    v::Float64      # Particle - velocity
    m::Float64      # Particle - mass
    mat::Material   # Particle - material
    lₚ::Float64     # Particle - half-size
    Vₒ::Float64     # Particle - initial volume
    V::Float64      # Particle - current volume
    σ::Float64      # Particle - stress
    ε::Float64      # Particle - strain
    F::Float64      # Particle - deformation gradient

    function Particle(x, m, mat, lₚ, v)
        Vₒ = 2*lₚ
        V = Vₒ
        σ = 0 
        ε = 0
        F = 1
        return new(x, v, m, mat, lₚ, Vₒ, V, σ, ε, F)
    end
end