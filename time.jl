include("domain.jl")

Base.@kwdef mutable struct Time
    domain::Domain
    tₙ::Float64
    tₒ::Float64
    t::Float64
    Δt::Float64
    pct::Float64
    Time(domain, tₙ, pct) = new(domain, tₙ, 0, 0, 0, pct)
end

function computeΔt(time::Time)
    c = -1
    for mat in time.domain.materials
        cc = sqrt(mat.E/mat.ρ)
        if (cc > c)
            c = cc
        end
    end

    time.Δt = time.pct*time.domain.grid.L/c
end

function advance(time::Time)

    # 1. Particle-To-Grid Mapping (P2G)
    reset(time.domain.grid)

    b = time.domain.b
    L = time.domain.grid.L

    for p in time.domain.particles

        e = mapParticleToElement(time.domain.grid, p)

        for n in e.nodes
            if !n.isLocked
                S = N(n, p.x, p.lₚ, L)
                ∇S = ∇N(n, p.x, p.lₚ, L)

                n.m += p.m * S
                n.p += (p.m * p.v) * S

                # Volume force
                n.f += p.m * b * S

                # Internal force
                n.f -= p.σ * p.V * ∇S
            end
        end
    end

    # 2. Update of nodal momentum
    for n in time.domain.grid.nodes
        n.p += time.Δt*n.f
    end

    # 3. Grid-To-Particle (G2P)
    for p in time.domain.particles
        e = mapParticleToElement(time.domain.grid, p)

        for n in e.nodes
            if n.m > 0
                if !n.isLocked
                    S = N(n, p.x, p.lₚ, L)

                    p.v += time.Δt * n.f/n.m * S
                    p.x += time.Δt * n.p/n.m * S
                end
            end
        end
    end

    # 4. Computation of particles stresses
    for p in time.domain.particles
        e = mapParticleToElement(time.domain.grid, p)

        Lₚ = 0
        for n in e.nodes
            if n.m > 0
                if !n.isLocked
                    ∇S = ∇N(n, p.x, p.lₚ, L)
                    Lₚ += n.p/n.m*∇S
                end
            end
        end

        p.F *= 1 + time.Δt * Lₚ

        p.V = p.F * p.Vₒ

        Δε = time.Δt * Lₚ

        p.ε += Δε

        p.σ += p.mat.E * Δε
    end

    time.t += time.Δt
end