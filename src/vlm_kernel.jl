function hshoe_kernel(
	cx,
	cy,
	cz,
	x1,
	x2,
	x3,
	p1x,
	p1y,
	p1z,
	p2x,
	p2y,
	p2z;
	M∞::Real = 0.0
)

	β = sqrt(1.0 - M∞ ^ 2)

	ax = cx .- p1x'
	ay = (cy .- p1y') .* β
	az = (cz .- p1z') .* β
	
	bx = cx .- p2x'
	by = (cy .- p2y') .* β
	bz = (cz .- p2z') .* β

	na = @. sqrt(ax ^ 2 + ay ^ 2 + az ^ 2) + 1e-10
	nb = @. sqrt(bx ^ 2 + by ^ 2 + bz ^ 2) + 1e-10

	ab = @. ax * bx + ay * by + az * bz

	mult = @. (1.0 / na + 1.0 / nb) / (ab + na * nb)

	d2 = @. na - ax * x1 - ay * x2 - az * x3
	d3 = @. nb - bx * x1 - by * x2 - bz * x3

	mult2 = @. 1.0 / (na * d2)
	mult3 = @. 1.0 / (nb * d3)
	
	vx = @. (
		(ay * bz - az * by) * mult +
		(ay * x3 - az * x2) * mult2 -
		(by * x3 - bz * x2) * mult3
	) / (4 * π * β ^ 2)
	vy = @. (
		(az * bx - ax * bz) * mult +
		(az * x1 - ax * x3) * mult2 -
		(bz * x1 - bx * x3) * mult3
	) / (4 * π * β)
	vz = @. (
		(ax * by - ay * bx) * mult +
		(ax * x2 - ay * x1) * mult2 -
		(bx * x2 - by * x1) * mult3
	) / (4 * π * β)

	[vx, vy, vz]
	
end

function vlm_kernel(
	cx,
	cy,
	cz,
	p1x,
	p1y,
	p1z,
	p2x,
	p2y,
	p2z;
	M∞::Real = 0.0
)

	β = sqrt(1.0 - M∞ ^ 2)

	ax = cx .- p1x'
	ay = (cy .- p1y') .* β
	az = (cz .- p1z') .* β
	
	bx = cx .- p2x'
	by = (cy .- p2y') .* β
	bz = (cz .- p2z') .* β

	na = @. sqrt(ax ^ 2 + ay ^ 2 + az ^ 2) + 1e-10
	nb = @. sqrt(bx ^ 2 + by ^ 2 + bz ^ 2) + 1e-10

	ab = @. ax * bx + ay * by + az * bz

	mult = @. (1.0 / na + 1.0 / nb) / (ab + na * nb)

	vx = @. (
		ay * bz - az * by
	) * mult / (4 * π * β ^ 2)
	vy = @. (
		az * bx - ax * bz
	) * mult / (4 * π * β)
	vz = @. (
		ax * by - ay * bx
	) * mult / (4 * π * β)

	[vx, vy, vz]
	
end
