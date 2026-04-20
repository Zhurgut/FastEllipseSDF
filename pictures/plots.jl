using Plots
using LinearAlgebra

include("../reference.jl")

ellipse(x; x_radius=4) = sqrt(max(0, 1 - (x/x_radius)^2))
ddx_ellipse(x; x_radius=4) = -x / (x_radius^2 * ellipse(x, x_radius=x_radius) + 1e-9)
ddx2_ellipse(x; x_radius=4) = ( -ellipse(x, x_radius=x_radius) + x*ddx_ellipse(x; x_radius=x_radius) ) / ( x_radius^2 * ellipse(x, x_radius=x_radius)^2 + 1e-9)
ellipse_normal_slope(x; x_radius=4) = -1/ddx_ellipse(x; x_radius=x_radius)

dot(p1, p2) = p1.x*p2.x + p1.y*p2.y

distance(p1, p2) = sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)

evolute(x; x_radius=4) = if 0 <= x <= (x_radius^2-1)/x_radius
     -((x_radius^2 - 1)^(2/3) - x_radius^(2/3)*x^(2/3))^(3/2)
else
    0.0
end

function parabola1(x; x_radius=4)
    s = x_radius
    A = s^2 / (1 - s^2)
    B = 2s
    C = 1 - s^2
    return parabola(A, B, C, x)
end

function parabola_coeffs(x1, y1, x2, y2, x3, y3)

    A = [x1^2  x1  1.0;
         x2^2  x2  1.0;
         x3^2  x3  1.0]

    y = [y1; y2; y3]
    coeffs = A \ y   # solves the linear system
    return coeffs[1], coeffs[2], coeffs[3]   # a, b, c
end

function parabola2(x; x_radius=4)
    s = x_radius
    y1 = 1 - s^2
    x3 = (s^2-1)/s
    x2 = 0.5x3
    a, b, c = parabola_coeffs(0, y1, x2, evolute(x2, x_radius=x_radius), x3, 0)
    return parabola(a, b, c, x)
end

begin
R = 4
X = 0:0.01:R+1
plot(X, [ellipse.(X, x_radius=R) evolute.(X, x_radius=R) parabola1.(X, x_radius=R) parabola2.(X, x_radius=R)], xlims=(0, R), ylims=(1 - R^2, 1))
end

function approximate_normal_slope(x; x_radius=4)

    y = ellipse(x, x_radius=x_radius)

    @assert x_radius >= 1
    @assert x >= 0
    @assert y >= 0

    y_radius = 1

    s = x_radius
    A = s^2 / (1 - s^2)
    B = 2s
    C = 1 - s^2
    parabola_peak_x = (s^2-1)/s

    a = -A
    b = 2*A*x
    c = B*x + C - y

    tangent_x = solve_quadratic(a, b, c)
    tangent_x = clamp(tangent_x, 0, parabola_peak_x)

    p1 = (x=tangent_x, y=parabola(A, B, C, tangent_x))
    p2 = (x=x, y=y)

    rise = p2.y - p1.y
    run = p2.x - p1.x

    return rise/run

end

X = 0:0.01:4.5
plot(X, ellipse.(X), ratio=1)
plot(X, [ellipse_normal_slope.(X) approximate_normal_slope.(X)], 
    labels=["target" "approximation"], 
    ylims=(0,25),
    xlabel="x",
    ylabel="slope",
    title="Slope of Normal"

)




function true_ellipse_sdf(x_radius, x, y)
    x = abs(x)
    y = abs(y)


    @assert x_radius >= 1
    @assert x >= 0
    @assert y >= 0

    y_radius = 1

    f(xk) = -2(x - xk + y * ddx_ellipse(xk, x_radius=x_radius) + xk / x_radius^2)
    df(xk) = -2(-1 + y * ddx2_ellipse(xk, x_radius=x_radius) + 1 / x_radius^2)

    xₖ = sqrt(0.5) * x_radius

    prev_L = f(xₖ) 

    for i=1:30

        L = f(xₖ) 
        
        if i >= 5 && abs(L) >= abs(prev_L)
            # println("x= ", xₖ)
            break
        end

        prev_L = L

        xₖ = clamp(xₖ - L / df(xₖ), 1e-14, x_radius-1e-14)

        
        

    end

    p = (x=x, y=y)

    shifted_p = (x=p.x / x_radius, y=p.y)

    c2 = dot(shifted_p, shifted_p) - 1

    signum = sign(c2)

    return signum * distance((x=xₖ, y=ellipse(xₖ, x_radius=x_radius)), p)

end



begin
    R = 4
    X = -0.1:0.01:6
    Y = -0.1:0.01:3

    out = Matrix{Float64}(undef, length(Y), length(X))
    for r in 1:length(Y), c in 1:length(X)
        x = X[c]
        y = Y[r]

        T = true_ellipse_sdf(R, x, y)
        F = sd_ellipse(x, y, R)
        out[r, c] = 2 * abs(T - F) / (abs(T + F) + 1e-12)

    end

    

    p = heatmap(X, Y, out, title="Relative Error", size=(1200, 800), xlims=(X[1], X[end]), ylims=(Y[1], Y[end]))
    plot!(p, X, ellipse.(X, x_radius=R), legend=false, color=:green)
    savefig(p, "pictures/rel_error.png")
    p

end


begin
    R = 4
    X = -0.1:0.01:6
    Y = -0.1:0.01:3

    out = Matrix{Float64}(undef, length(Y), length(X))
    for r in 1:length(Y), c in 1:length(X)
        x = X[c]
        y = Y[r]

        # println("$r, $c -> $y, $x")
        T = true_ellipse_sdf(R, x, y)
        F = sd_ellipse(x, y, R)
        # out[r, c] = 2 * abs((T - F) / (T + F))
        # out[r, c] = log(2 * abs((T - F) / (T + F)) + 1e-12)
        out[r, c] = abs(T - F)
        # out[r, c] = T
    end


    p = heatmap(X, Y, out, title="Absolute Error", size=(1200, 800), xlims=(X[1], X[end]), ylims=(Y[1], Y[end]))
    plot!(p, X, ellipse.(X, x_radius=R), legend=false, color=:green)
    savefig(p, "pictures/abs_error.png")
    p

end





data = [
    (32.5, 34.7, "Newton Trig, 5 it."), # https://www.shadertoy.com/view/4lsXDN
    (33.4, 35.2, "Newton No-Trig, 5 it."), # https://www.shadertoy.com/view/tttfzr
    (45.4, 55.4, "Newton No-Trig optim. 3 it."), # https://www.shadertoy.com/view/tt3yz7
    (59,   56.4, "Analytical"), # https://www.shadertoy.com/view/4sS3zz
    (64.5, 89, "CheapEvolute")
]

fps_general = [x[1] for x in data]
fps_fixed = [x[2] for x in data]
algorithm = [x[3] for x in data]

begin
    B = bar(title="Performance Comparison General", ylabel="FPS", xrotation=30, xlabel=" ", ylims=(0,100), legend=:topleft)
    bar!(B, algorithm[1:3], fps_general[1:3], label="Newton")
    bar!(B, [algorithm[4]], [fps_general[4]], label="Analytical")
    bar!(B, [algorithm[5]], [fps_general[5]], label="CheapEvolute")
    savefig(B, "pictures/fps_general.png")
    B
end

begin
    B = bar(title="Performance Comparison Fixed", ylabel="FPS", xrotation=30, xlabel=" ", ylims=(0,100))
    bar!(B, algorithm[1:3], fps_fixed[1:3], label="Newton")
    bar!(B, [algorithm[4]], [fps_fixed[4]], label="Analytical")
    bar!(B, [algorithm[5]], [fps_fixed[5]], label="CheapEvolute")
    savefig(B, "pictures/fps_fixed.png")
    B
end