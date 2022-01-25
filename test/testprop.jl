@assert begin
    res = PropBiology.simulation(
        [0.0, 1.0] ./ 4,
        [0.2, 0.1] ./ 4,
        9000.0 / 60.0;
        γ = [deg2rad(45.0), deg2rad(15.0)] .+ deg2rad(15.0),
        V∞ = 10.0
    )

    true
end
