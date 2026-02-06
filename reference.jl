


dot(x1, y1, x2, y2) = x1*x2 + y1*y2
mix(p1, p2, z) = (1-z)*p1 + z*p2
distance(x1, y1, x2, y2) = sqrt((x1 - x2)^2 + (y1 - y2)^2)

parabola(a, b, c, x) = a*x^2 + b*x + c

solve_quadratic(a, b, c; ϵ=0) = (-b - sqrt(max(0, b^2 - 4*a*c))) / (2a+ϵ)


function sd_ellipse(px, py, x_radius, y_radius=1, elps_center_x=0, elps_center_y=0)
    @assert x_radius > 0
    @assert y_radius > 0

    px = px - elps_center_x
    py = py - elps_center_y

    px = abs(px)
    py = abs(py)

    if y_radius > x_radius
        px, py = py, px
        x_radius, y_radius = y_radius, x_radius
    end

    @assert x_radius >= y_radius

    px = px / y_radius
    py = py / y_radius

    r = x_radius / y_radius

    parabola_peak_x = (r^2 - 1) / r

    A = -r / (parabola_peak_x + 1e-6)
    B = 2r
    C = 1 - r^2

    a = -A
    b = 2*A*px
    c = B*px + C - py

    tx = solve_quadratic(a, b, c)
    tx = clamp(tx, 0, parabola_peak_x - 1e-4)
    ty = parabola(A, B, C, tx)

    p_shift_x = px / r
    p_shift_y = py

    t_shift_x = tx / r
    t_shift_y = ty

    dif_x = t_shift_x - p_shift_x
    dif_y = t_shift_y - p_shift_y

    a2 = dif_x^2 + dif_y^2
    b2 = 2*dot(dif_x, dif_y, p_shift_x, p_shift_y)
    c2 = p_shift_x^2 + p_shift_y^2 - 1

    z = solve_quadratic(a2, b2, c2, ϵ=1e-20)

    # the point between t_shift and p_shift on the unit circle
    cx = mix(p_shift_x, t_shift_x, z)
    cy = mix(p_shift_y, t_shift_y, z) 

    # the point between t and p on the ellipse
    ex = r * cx
    ey = cy 

    return sign(c2) * y_radius * distance(ex, ey, px, py)

end


using Test
@testset "Test CheapEvolute Ellipse SDF" begin
    @test isapprox(sd_ellipse(0.0, 0.0, 4.0) , -1.0, atol=1e-4) 
    @test isapprox(sd_ellipse(4.0, 0.0, 4.0) , 0.0, atol=1e-4) 
    @test isapprox(sd_ellipse(0.0, 1.0, 4.0) , 0.0, atol=1e-4) 
    @test isapprox(sd_ellipse(3.75, 0.0, 4.0), -0.25, atol=1e-4) 

    @test isapprox(sd_ellipse(0.0, 0.0, 1.0) , -1.0, atol=1e-4) 
    @test isapprox(sd_ellipse(4.0, 0.0, 1.0) , 3.0, atol=1e-4) 
    @test isapprox(sd_ellipse(0.0, 1.0, 1.0) , 0.0, atol=1e-4) 
    @test isapprox(sd_ellipse(1.0, 1.0, 1.0) , sqrt(2) - 1, atol=1e-4) 
end