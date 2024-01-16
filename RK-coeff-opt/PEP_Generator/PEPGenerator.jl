using BSeries, SymPy

#set order
q = 7
# generate series
series = bseries(q) do t, series
    return symbols("a_$(butcher_representation(t))", real=true)
end

# Our interest relies on the modified equation
modi = modified_equation(series)
renormalize!(modi)

# show the result for every condition (this expression involves from 1 to 4 trees).
# The conditions can be obtained from David's Jupyter Notebook 'Energy-preserving
# conjugates'. 
println(modi[rootedtree([1,2,3,4,4,4,3])] + modi[rootedtree([1,2,3,4,3,2,2])] - modi[rootedtree([1,2,3,3,3,2,3])])