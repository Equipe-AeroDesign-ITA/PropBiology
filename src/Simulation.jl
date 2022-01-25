_interpolate(
    ȳ::AbstractVector, 
    v::Real, 
    stations::AbstractVector;
    interp = LinearInterpolation
) = v .* ones(Float64, length(stations))
_interpolate(
    ȳ::AbstractVector, 
    v::AbstractVector, 
    stations::AbstractVector;
    interp = LinearInterpolation
) = let A = interp(ȳ, v)
    A.(stations)
end
_interpolate(
    ȳ::AbstractVector, 
    v::AbstractMatrix, 
    stations::AbstractVector;
    interp = LinearInterpolation
) = permutedims(
    hcat(
        [
            _interpolate(ȳ, v[i, :], stations; interp = interp) for i = 1:size(v, 1)
        ]...
    )
)

_trapz(
    ys::AbstractVector, xs::AbstractVector
) = let dx = xs[2:end] .- xs[1:(end - 1)]
    sum((ys[1:(end - 1)] .+ ys[2:end]) .* dx) / 2
end
_trapz(
    ys::AbstractMatrix, xs::AbstractVector
) = [
    _trapz(ys[i, :], xs) for i = 1:size(ys, 1)
]

_rowcross(
    u::AbstractMatrix,
    v::AbstractMatrix
) = permutedims(
    hcat(
        (u[2, :] .* v[3, :] .- u[3, :] .* v[2, :]),
        (u[3, :] .* v[1, :] .- u[1, :] .* v[3, :]),
        (u[1, :] .* v[2, :] .- u[2, :] .* v[1, :])
    )
)

const default_stations = collect(0.0:0.01:1.0)

"""
    ```
    function simulation(
        ys::AbstractVector,
        cs,
        Ω::Real;
        n_blades::Int64 = 2,
        xs = 0.0,
        CL0 = 0.0,
        CLα = 2 * π,
        Cm0 = 0.0,
        Cmα = 0.0,
        polar = zeros(3),
        e = 0.0,
        γ = 0.0,
        V∞::Real = 0.0,
        a::Real = 340.0,
        ρ::Real = 1.225,
        n_iterations::Int64 = 20,
        h::Real = 1e-3,
        tip_loss::Bool = true,
        stations::AbstractVector = default_stations,
        interp = LinearInterpolation
    )
    ```

Simulate propeller flow

* `ys`: positions of stations in the radial coordinate
    in which properties are defined
* `cs`: chords for each control station in `ys` (or a scalar,
    for a chord applicable to all stations)
* `Ω`: rotational speed (rad/s)
* `n_blades`: number of blades. Defaults to 2
* `xs`: displacements of control section quarter chords from 
    radial line. Can be used to apply sweep to the blade, and leads
    to more accurate load predictions
* `CL0, CLα, Cm0, Cmα`: airfoil properties at each control section
* `polar`: polars for airfoils at each control section 
    (as a `(3, N)` matrix), or a single polar for 
    all sections (as a vector).

    Should follow format: `CD = polar[1] * CL ^ 2 + polar[2] * CL + polar[3]`
* `e`: displacement of elastic axis from quarter-chord, as a fraction 
    of the chord. Scalar for the blade, or a vector for each control section
* `γ`: local incidences, in radians. Scalar for the whole blade, or
    vector for each control section
* `V∞`: freestream velocity
* `a`: sound speed
* `ρ`: air density
* `n_iterations`: number of Newton iterations for BEMM
* `h`: step for linearization in Newton method for BEMM
* `tip_loss`: boolean indicating whether or not to apply
    Glauert's tip loss factor
* `stations`: vector ranging from 0 to 1 with positions for
    integration points
* `interp`: interpolation mechanism for obtaining sectional 
    properties
"""
function simulation(
    ys::AbstractVector,
    cs,
    Ω::Real;
    n_blades::Int64 = 2,
    xs = 0.0,
    CL0 = 0.0,
    CLα = 2 * π,
    Cm0 = 0.0,
    Cmα = 0.0,
    polar = zeros(3),
    e = 0.0,
    γ = 0.0,
    V∞::Real = 0.0,
    a::Real = 340.0,
    ρ::Real = 1.225,
    n_iterations::Int64 = 20,
    h::Real = 1e-3,
    tip_loss::Bool = true,
    stations::AbstractVector = default_stations,
    interp = LinearInterpolation
)

    polar = (
        isa(polar, AbstractVector) ?
        hcat(
            [polar for _ = 1:length(ys)]...
        ) :
        polar
    )

    R = maximum(ys) - minimum(ys)
    D = 2 * R

    A = R ^ 2 * π

    ȳ = ys ./ R
    ys = _interpolate(ȳ, ys, stations)

    cs = _interpolate(ȳ, cs, stations; interp = interp)
    xs = _interpolate(ȳ, xs, stations; interp = interp)
    CL0 = _interpolate(ȳ, CL0, stations; interp = interp)
    CLα = _interpolate(ȳ, CLα, stations; interp = interp)
    Cm0 = _interpolate(ȳ, Cm0, stations; interp = interp)
    Cmα = _interpolate(ȳ, Cmα, stations; interp = interp)
    γ = _interpolate(ȳ, γ, stations; interp = interp)
    e = _interpolate(ȳ, e, stations; interp = interp)

    polar = _interpolate(ȳ, polar, stations; interp = interp)

    dy = ys[2:end] .- ys[1:(end - 1)]

    Vr = @. ys * Ω

    f = Vi -> begin
        ṁ = A * (Vi + V∞) * ρ

        ϕ = atan.(V∞ + Vi, Vr)

        α = γ .- ϕ

        # Glauert's tip loss factor
        F = @. 2 * acos(exp(- n_blades * ((R - ys + 1e-10) / (ys + 1e-10)) / (2 * sin(ϕ)))) / π

        Vmod = @. sqrt((V∞ + Vi) ^ 2 + Vr .^ 2)
        M = Vmod ./ a

        β = sqrt.(1.0 .- M .^ 2)

        CL = (α .* CLα .+ CL0) .* (tip_loss ? F : 1.0)
        CD = polar[1, :] .* CL .^ 2 .+ polar[2, :] .* CL .+ polar[3, :]

        CL = CL ./ β
        CD = CD ./ β

        qc = @. ((V∞ + Vi) ^ 2 + Vr ^ 2) * ρ * cs / 2

        CX = @. CL * cos(ϕ) - CD * sin(ϕ)
        CY = @. CL * sin(ϕ) + CD * cos(ϕ)

        T = let qcCX = CX .* qc
            sum((qcCX[1:(end - 1)] .+ qcCX[2:end]) .* dy) * n_blades / 2
        end

        T - ṁ * Vi * 2
    end

    h = h * Ω * R

    fp = Vi -> (f(Vi + h) - f(Vi)) / h

    Vi = 0.0

    for _ = 1:n_iterations
        Vi = Vi - f(Vi) / fp(Vi)
    end

    ṁ = A * (Vi + V∞) * ρ

    ϕ = atan.(V∞ + Vi, Vr)

    α = γ .- ϕ

    # Glauert's tip loss factor
    F = @. 2 * acos(exp(- n_blades * ((R - ys + 1e-10) / (ys + 1e-10)) / (2 * sin(ϕ)))) / π

    Vmod = @. sqrt((V∞ + Vi) ^ 2 + Vr .^ 2)
    M = Vmod ./ a

    β = sqrt.(1.0 .- M .^ 2)

    CL = (α .* CLα .+ CL0) .* (tip_loss ? F : 1.0)
    CD = polar[1, :] .* CL .^ 2 .+ polar[2, :] .* CL .+ polar[3, :]

    CL = CL ./ β
    CD = CD ./ β

    qc = @. ((V∞ + Vi) ^ 2 + Vr ^ 2) * ρ * cs / 2

    CX = @. CL * cos(ϕ) - CD * sin(ϕ)
    CY = @. CL * sin(ϕ) + CD * cos(ϕ)

    Cm = (α .* Cmα .+ Cm0) .* (tip_loss ? F : 1.0) ./ β

    T = let qcCX = qc .* CX
        sum(dy .* (qcCX[1:(end - 1)] .+ qcCX[2:end])) * n_blades / 2
    end
    Q = let qcCYy = qc .* CY .* ys
        sum(dy .* (qcCYy[1:(end - 1)] .+ qcCYy[2:end])) * n_blades / 2
    end

    M = Cm .* cs .* qc

    quarter_chords = permutedims(
        [
            xs ys zero(xs)
        ]
    )
    elastic_centers = quarter_chords .+ permutedims(
        [
            (cs .* cos.(γ) .* e) zero(ys) (- cs .* sin.(γ) .* e)
        ]
    )

    dF = permutedims([(CY .* qc) zero(ys) (CX .* qc)])

    forces = hcat(
        [
            _trapz(
                dF[:, i:end], ys[i:end]
            ) for i = 1:length(ys)
        ]...
    )
    moments = let dM = permutedims([zero(ys) M zero(ys)])
        hcat(
            [
                _trapz(
                    dM[:, i:end] .+ _rowcross(
                        quarter_chords[:, i:end] .- elastic_centers[:, i],
                        dF[:, i:end]
                    ),
                    ys[i:end]
                ) for i = 1:length(ys)
            ]...
        )
    end

    My = moments[2, :]
    Mx = moments[1, :] .* cos.(γ) .- moments[3, :] .* sin.(γ)
    Mz = moments[3, :] .* cos.(γ) .+ moments[1, :] .* sin.(γ)

    Vx = - forces[3, :] .* sin.(γ) .+ forces[1, :] .* cos.(γ)
    Vz = forces[3, :] .* cos.(γ) .+ forces[1, :] .* sin.(γ)

    (
        T = T,
        Q = Q,
        CL = CL,
        Cm = Cm,
        CD = CD,
        forces = forces,
        moments = moments,
        Mx = Mx, 
        My = My,
        Mz = Mz,
        Vx = Vx,
        Vz = Vz,
        ϕ = ϕ,
        r = ys,
        c = cs,
        Vi = Vi
    )

end
