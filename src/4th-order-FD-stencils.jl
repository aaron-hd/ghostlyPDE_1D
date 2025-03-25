# First derivative (4th order)
function first_derivative(f, n, dx)
    return (
        - f[n + 2]
        + 8 * f[n + 1]
        - 8 * f[n - 1]
        + f[n - 2]
    ) / (12 * dx)
end

# Second derivative (4th order)
function second_derivative(f, n, dx)
    return (
        - f[n + 2]
        + 16 * f[n + 1]
        - 30 * f[n]
        + 16 * f[n - 1]
        - f[n - 2]
    ) / (12 * dx^2)
end
