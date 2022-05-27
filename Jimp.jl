include("time.jl")

export Jimp

# MATERIAL DATA
E = 200.0e9         # Elastic modulus
ρ = 7800.0          # Mass density

# PHYSICAL DOMAIN DATA
l = 25.0            # Length
vₒ = .1             # Initial velocity multiplier
n = 1               # Mode of vibration
βₙ = (2*n-1)*π/2/l  # Eigenvalue (mode of vibration n)
ωₙ = βₙ * sqrt(E/ρ) # Frequency of oscilation (mode of vibration n)
Tₙ = 2*π/ωₙ         # Period (mode of vibration n)

# GRID DATA
lₑ = 1.0            # Grid element size

# PARTICLE DATA
ppe = 2             # Number of particles per element

# TIME DATA
tₒ = 0              # Initial time
tₙ = 5*Tₙ           # Final time
pct = 0.1           # Percentage of critical time step


# MODEL CREATION
material = Material(E, ρ)

domain = Domain(genGrid(lₑ, -lₑ, l+2*lₑ))

lockNodeByPosition(domain.grid, 0.0)
lockNodeByPosition(domain.grid, -lₑ)

genParticles(domain, material, ppe, 0.0, l)
for p in domain.particles
    p.v = vₒ*sin(βₙ*p.x) 
end

# TIME INTEGRATION
time = Time(domain, tₙ, pct)
computeΔt(time)

num_steps = convert(Int64, ceil((tₙ-tₒ)/time.Δt))
num_print_steps = 100
step = convert(Int64, ceil(num_steps / num_print_steps))
step_acc = step
progress = 100/num_print_steps

t = zeros(num_steps)
analytical = zeros(num_steps)
vCM = zeros(num_steps)

M = 0
for p in time.domain.particles
    global M += p.m
end

for i=1:num_steps

    # Validation metrics
    t[i] = time.t
    analytical[i] = vₒ/βₙ/l*cos(ωₙ*t[i])
    vCM[i] = 0
    for p in time.domain.particles
        vCM[i] += p.v*p.m/M
    end

    advance(time)

    if i > step_acc
        println("PROGRESS: ",  round(progress, digits=2), "\t|\tTIME: ", round(time.t, digits=4), " s")
        global progress += 100/num_print_steps
        global step_acc += step
    end
end

using Plots

plt = plot(t, analytical, label="Analytical", lw=10)
plot!(plt, t, vCM, label="JIMP", lw=1, seriestype = :scatter)

Plots.display(plt)

print("FINISHED")