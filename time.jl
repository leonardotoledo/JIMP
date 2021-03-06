include("domain.jl")

Base.@kwdef mutable struct Time
    domain::Domain  # Time - domain
    tₙ::Float64     # Time - final time
    tₒ::Float64     # Time - initial time
    t::Float64      # Time - current time
    Δt::Float64     # Time - time step
    pct::Float64    # Time - percentage of critical time step
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

    time.Δt = time.pct*time.domain.grid.lₑ/c
end

function advance(time::Time)

    # 1. Particle-To-Grid Mapping (P2G)
    reset(time.domain.grid)

    b = time.domain.b
    lₑ = time.domain.grid.lₑ

    for p in time.domain.particles

        e = mapParticleToElement(time.domain.grid, p)

        for n in e.nodes
            if !n.isLocked
                S = N(n, p.x, p.lₚ, lₑ)
                ∇S = ∇N(n, p.x, p.lₚ, lₑ)

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
                    S = N(n, p.x, p.lₚ, lₑ)

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
                    ∇S = ∇N(n, p.x, p.lₚ, lₑ)
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