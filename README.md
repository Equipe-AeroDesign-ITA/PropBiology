# PropBiology.jl

A simple, blade element momentum theory code for propeller performance estimation.

To run a simulation, use:

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
