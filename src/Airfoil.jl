"""
Data type defining an airfoil

* `polar` should be an array `[a, b, c]` such that `CD = CL ^ 2 * a + CL * b + c`
"""
struct Airfoil{Fg <: Real}
	CL0::Fg
	CLα::Fg
	Cm0::Fg
	Cmα::Fg
	polar::Vector{Fg}
end
Airfoil(CL0::Fg, args...) where Fg = Airfoil{Fg}(CL0, args...)

Base.:+(afl1::Airfoil, afl2::Airfoil) = Airfoil(
    (afl1.CL0 + afl2.CL0) / 2,
    (afl1.CLα + afl2.CLα) / 2,
    (afl1.Cm0 + afl2.Cm0) / 2,
    (afl1.Cmα + afl2.Cmα) / 2,
    (afl1.polar .+ afl2.polar) ./ 2
)
Base.:*(η::Number, afl::Airfoil) = Airfoil(
    η * afl.CL0,
    η * afl.CLα,
    η * afl.Cm0,
    η * afl.Cmα,
    η .* afl.polar
)

"""
Default, inviscid flat plate airfoil
"""
flat_plate = Airfoil{Float64}(0.0, 2 * π, 0.0, 0.0, [0.0, 0.0, 0.0])
