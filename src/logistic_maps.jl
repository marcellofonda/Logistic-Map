"Module containing functions to treat discrete nonlinear dynamical systems. By Marcello Fonda"
module LogisticMaps

using Plots

export logistic_map, draw_cobweb, find_period, find_fixed_points, fixed_points_as_tuples, Lyapunov_average

"""
Dictionary containing all the available maps.

Implemented keys: "logistic", "tent", "sine", "bimodal".
"""
maps = Dict(
    "logistic" => (x,r) -> 4r * (1-x)x,
    "tent" => (x,r) -> x ≤ .5 ? 2r * x : 2r - 2r * x,
    "sine" => (x,r) -> r * sin( π * x ),
    "bimodal" => (x,r) -> x < .5 ? 4r * (1-2x)*2x : 4r * (1-(2x-1))*(2x-1),
    "m_like" => (x,r) -> 3.8 * r * (.1 + abs(x - 0.5)) * (1 - (2x - 1)^2)#(x,r) -> 2.99r * (1.1 - sin(π * x)) * (1 - (2x - 1)^2)
)

"""
Dictionary containing available map derivatives.

Implemented keys: "logistic", "tent", "sine".
"""
derivatives = Dict(
    "logistic" => (x,r) -> 4r * (1 - 2x),
    "tent" => (x,r) -> x < .5 ? 2r : -2r,
    "sine" => (x,r) -> r * π * cos(π * x),
    "bimodal" => (x,r) -> x < 0.5 ? -16r*x + 8(1 - 2x)r : 8(2 - 2x)r - 8(2x - 1)r,
    "m_like" => (x,r) -> x > 0.5 ? 3.8 * r * (1 * (1 - (2x - 1)^2) - 8 * (2x - 1) * (.1 + x - 0.5)) :
    3.8 * r * (-1 * (1 - (2x - 1)^2) - 8 * (2x - 1) * (.1 + 0.5 - x))#x < 0.5 ? 3.8 * r * (1 * (1 - (2x - 1)^2) - 8 * (2x - 1) * (.1 + 0.5 - x)) : 3.8 * r * (1 * (1 - (2x - 1)^2) - 8 * (2x - 1) * (.1 + x - 0.5))#2.99r * ((-cos(π * x) * π) * (1 - (2x - 1)^2) + (1.1 - sin(π * x)) * -8 * (2x - 1))
)

"""
    logistic_map(x₀, r, N_steps; map="logistic")

Run the dynamical system associated to a `map` from [`maps`](@ref) (defaults to `"logistic"`), starting from point `x₀`, with parameter `r`, for `N_steps` discrete steps. 
Return a vector containing all visited positions.
"""
function logistic_map(x₀, r, N_steps; map="logistic")
    haskey(maps, map) || throw(ArgumentError("""Selected map "$map" not available"""))

    positions = [x₀]
    for _ in 1:N_steps
        push!(positions,maps[map](positions[end], r))
    end
    positions
end

"""
    cobweb_path(points, r; map="logistic")

Create a path for the cobweb plot. Given a set of `points`, the parameter `r`, 
and optionally the `map` type, this function generates a sequence of points 
that represent the cobweb plot path for the given dynamical system.
"""
function cobweb_path(points, r; map="logistic")
    haskey(maps, map) || throw(ArgumentError("""Selected map "$map" not available"""))

    point_and_following(point,r) = [(point,point), (point, maps[map](point,r))]
    vcat(point_and_following.(points,r)...)
end

"""
    draw_cobweb(points, r; start=1, map="logistic")

Generate and display a cobweb plot for a given dynamical system. 

The function takes a vector of `points` that represent the trajectory of the system, 
and the parameter `r`. Optional arguments include `start`, which specifies the starting 
point in the trajectory for the plot, and `map`, which defaults to "logistic" and determines 
the type of map used for the dynamical system.

The cobweb plot visually represents the iterative process of the map, showing how each point 
in the trajectory is mapped to the next. It's particularly useful for understanding the 
behavior of iterative systems, such as convergence to fixed points or periodic orbits.

# Examples
```julia
points = logistic_map(0.1, 3.5, 100)
draw_cobweb(points, 3.5)
```
"""
function draw_cobweb(points, r; start=1, map="logistic")
    haskey(maps, map) || throw(ArgumentError("""Selected map "$map" not available"""))

    x = 0:.01:1
    plot(x,maps[map].(x,r), aspect_ratio=1, legend=false, title="x₀=$(round(points[1],digits=5)), r=$(round(r,digits=5)), n=$(length(points)-1)")
    plot!(x,x)
    path = cobweb_path(points[start:end], r, map=map)
    plot!(path)
end


const MAX_N_POINTS=10_000

"""
    find_period(r; map="logistic")

Find the period of the equilibrium of the dynamical system associated to the given map (defaults to the logistic map), with the given parameter.
Return the period and the points visited.
If the period is greater than [`MAX_N_POINTS`](@ref), return `Inf`.
"""
function find_period(r; map="logistic")
    points = logistic_map(.2, r, MAX_N_POINTS, map=map)
    for i in 1:MAX_N_POINTS+1
        k = findfirst(==(points[i]), points)
        if k != i
            return i-k, points
        end
    end
    #println("infinito!")
    return Inf, points
end

"""
    find_fixed_points(r; map="logistic", n=1000)

Find the fixed points for the given map and parameter `r`. The `map` defaults to "logistic", 
and, if `include_non_periodic` is `true`, `n` specifies the number of points to consider if the period is infinite.
Returns a vector of fixed points.

See also [`find_period`](@ref)
"""
function find_fixed_points(r; map="logistic", n=1000, include_non_periodic=true)
    period, points = find_period(r, map=map) 
    fixed_points = period < Inf ? round.(points[end-period+1:end], digits=5) : include_non_periodic ? points[end-n:end] : []
end

"""
    fixed_points_as_tuples(r; map="logistic", n=1000)

Convert the fixed points obtained from `find_fixed_points` into tuples of 
the form (r, x), where `r` is the parameter and `x` is the fixed point.
"""
fixed_points_as_tuples(r; map="logistic", n=1000, include_non_periodic=true) = [(r, x) for x in find_fixed_points(r, map=map, n=n, include_non_periodic=include_non_periodic)]



# This section contains functions relative to the calculation of Lyapunov exponents

"""
    numerical_Lyapunov_exponent(r, N_steps, x₀; start=1, map="logistic", δx = .0001)

Calculate the numerical Lyapunov exponent for the given map and parameter `r`, 
starting from `x₀` for `N_steps`. `start` specifies the starting step for calculation, 
`map` is the type of map, and `δx` is the initial perturbation size.
Return the Lyapunov exponent for the map.
"""
function numerical_Lyapunov_exponent(r, N_steps, x₀; start=1, map="logistic", δx = .0001)
    # Initialize the positions as a (N_steps+1)×2 matrix
    positions = hcat(
        logistic_map(x₀,r, N_steps, map=map)[start:end], 
        logistic_map(x₀+δx,r, N_steps, map=map)[start:end]
    )'

    # Compute Deltaf map differences, obtain a vector and fetch nonzero values
    df = diff(positions, dims=1) |> filter(x -> (x≠0))
    # Stagger and divide difference vector
    df_ratio = df[2:end] ./ df[1:end-1]
    # Obtain logarithm and Lyapunov exponent
    sum(log.(abs.(df_ratio))) / length(df_ratio)
end

"""
    analytical_Lyapunov_exponent(points, r; start=1, stop=lastindex(points), map = "logistic")

Calculate the analytical Lyapunov exponent for the given map and parameter `r`.
`points` is the vector of points in the trajectory, `start` and `stop` define 
the range of points to consider, and `map` specifies the map type.
Returns the Lyapunov exponent based on the analytical derivative.
"""
function analytical_Lyapunov_exponent(points, r; start=1, stop=lastindex(points), map = "logistic")
    haskey(derivatives, map) || throw(ArgumentError("""Selected map "$map" does not have an analytically computed derivative"""))
    # Compute derivatives and fetch nonzero values
    derivs = derivatives[map].(points[start:stop], r) |> filter(x -> (x≠0))

    λ = sum(log.(abs.(derivs))) / length(derivs)
end

"""
    Lyapunov_average(x_range, r; start=1, stop=100, map = "logistic", mode="auto", δx=.0001)

Calculate the average Lyapunov exponent over a range of initial conditions `x_range` 
for the parameter `r`. `start` and `stop` define the step range, `map` is the map type, 
`mode` can be "auto", "numerical", or "analytical", and `δx` is the initial perturbation size.
Return the average Lyapunov exponent for the specified range.
"""
function Lyapunov_average(x_range, r; start=1, stop=100, map = "logistic", mode="auto", δx=.0001)
    haskey(maps, map) || throw(ArgumentError("""Selected map "$map" not available"""))
    mode ∈ ["auto", "numerical", "analytical"] || throw(ArgumentError("""Selected mode "$mode" not available"""))

    # Generate Lyapunov exponents for each starting point
    if (mode == "numerical") || ((mode == "auto") && (!haskey(derivatives,map))) 
        generated_exponents = numerical_Lyapunov_exponent.(r, stop, x_range, start=start, map=map, δx=δx)
    else
        generated_exponents = analytical_Lyapunov_exponent.(logistic_map.(x_range, r, stop, map=map), r, start=start, map=map)
    end
    
    # Return the averaged value
    return sum(generated_exponents)/length(generated_exponents)
end



end