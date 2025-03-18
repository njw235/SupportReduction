using Random
import Clarabel
using SumOfSquares
using DynamicPolynomials
using Plots

times = zeros(length(10:25:250;))
for i in [10:25:250;]

    @polyvar x
    f = rand(i) .* x.^[1:1:i;]
    g = rand(5) .* x.^[0:1:4;]
    S = @set x>=0 && x <= 1
    model = SOSModel(Clarabel.Optimizer)
    @variable(model, s)
    @objective(model, Max, s)
    @constraint(model, c, f >= s*g,S)
    optimize!(model)
    times[Int((i-10)/25) + 1] = solve_time(model)
end

plot([1:1:length([10:25:250;]);], times)
savefig("timeplt.png")