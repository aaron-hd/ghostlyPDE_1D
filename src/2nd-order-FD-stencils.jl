# First derivative (4th order)
function first_derivative(f, n, dx)
    return (
        + f[n + 1]
        - f[n - 1]
    ) / (2 * dx)
end

# Second derivative (4th order)
function second_derivative(f, n, dx)
    return (
        + f[n + 1]
        - 2 * f[n]
        + f[n - 1]
    ) / (dx^2)
end
