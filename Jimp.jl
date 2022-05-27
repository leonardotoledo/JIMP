include("time.jl")

export Jimp

# MATERIAL DATA
E = 200.0e9
ρ = 7800.0

# GRID DATA
Δx = 1.0
lₑ = 25.0

# PARTICLE DATA
ppe = 2
vₒ = .1
n = 1
β = (2*n-1)*π/2/lₑ
c = sqrt(E/ρ)
ω = β * c
T = 2*π/ω

# TIME DATA
t0 = 0
tf = 5*T
pct = 0.1

# MODEL CREATION
material = Material(E, ρ)

domain = Domain(genGrid(Δx, -Δx, lₑ+2*Δx))

lockNodeByPosition(domain.grid, 0.0)
lockNodeByPosition(domain.grid, -Δx)

genParticles(domain, material, ppe, 0.0, lₑ)
for p in domain.particles
    p.v = vₒ*sin(β*p.x) 
end

# TIME INTEGRATION
time = Time(domain, tf, pct)
computeΔt(time)

num_steps = convert(Int64, ceil((tf-t0)/time.Δt))
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
    analytical[i] = vₒ/β/lₑ*cos(ω*t[i])
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