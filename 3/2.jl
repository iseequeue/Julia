### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 89eb2182-cf2f-432f-bb83-65569e1a7d15
begin
	using LinearAlgebra
	using LaTeXStrings
end

# ╔═╡ de5b06e2-78b3-11ed-33d8-7377f554fd63
md"""
##### 1.Прямая и обратная подстановка
"""

# ╔═╡ d84a0ed1-07bc-4d85-9cb4-0765c5985537
begin
function LorU_solver(A, b)
    n = size(A, 1)
    x = zeros(n)

    if (size(A,1)!= size(A,2)) || (size(A,1)!=size(b, 1))
		error("wrong dimensionality")
        return x
    end
    
    if ( isapprox(prod(diag(A)), 0.0) )
        error("degenerate matrix")
        return x        
    end

	isL = istril(A)
	isU = istriu(A)

	println(isL, isU)
	if (isL)
		print(1)
		x[1] = b[1]/A[1,1]
		for i in 2:n
			s = 0.0
			for j in 1:i-1
				s+=A[i,j]*x[j]
			x[i] = (b[i] - s)/A[i,i]
			end
		end
		return x
	end
		
	if(isU)
		print(2)
		x[n] = b[n]/A[n,n]
		for i in 1:n-1
            s = 0.0            
            for j in 1:i                
                s+=A[n-i,n-j+1]*x[n-j+1]   
            end
            x[n-i] = (b[n-i] - s)/A[n-i,n-i]
        end
		return x
	end
	if(!isL && !isU)
		print(3)
		error("not lower- or upperdiagonal matrix")
		return x
	end		

end
end

# ╔═╡ e7e3f90a-9901-40a8-8d0f-8b02ade6abeb
let
	try
	# test forward(L,b)

	a = [
	    1 0 0; 
	    4 5 0; 
	    1 4 1]
	
	b = [1, 2, 3]
	
	x = LorU_solver(a, b)
		println(x)
	
	println("forward")
	println("solution = ", x)
	
	c =b - a*x	
	println("diff = ",  c)
	println()
	
	#test backward(U,b)
	a = [
	    1 9 3; 
	    0 1 1; 
	    0 0 1]
	
	b = [1, 2, 3]
	
	x = LorU_solver(a, b)
	
	println("backward")
	println("solution = ", x)
	
	c =b - a*x	
	println("diff = ",  c)
	println()
		#test errors
		a = [
		1 9 3; 
		0 1 1; 
		0 0 1]
		
		b = [1, 2, 3, 4]
		
		x = LorU_solver(a, b)
		
		println("backward")
		println("solution = ", x)
		
		c =b - a*x	
		println("diff = ",  c)

	catch e
		println(e)

	end

		
end

# ╔═╡ 552bd38c-3637-4cca-8ee7-467117ac27f6
md"""
##### 2.Метод прогонки
"""

# ╔═╡ 3833b8dc-2469-4195-9674-401bc3b545c2
begin

function tridiagsolve(a::AbstractVector, b::AbstractVector, c::AbstractVector, f::AbstractVector)
    
    n = size(b,1)
    x = similar(b, Float64, n)
    
    b_new = similar(b, Float64, n)
    f_new = similar(b, Float64, n)
    
    b_new[1] = b[1]
    f_new[1] = f[1]

    
    for i in 2:n
        if (isapprox(b_new[i-1], 0.0))
			error("degenerate matrix")
			return x
		end
        b_new[i] = b[i] - a[i-1]*c[i-1]     / b_new[i-1]
        f_new[i] = f[i] - a[i-1]*f_new[i-1] / b_new[i-1]
  
    end
    
    x[n] = f_new[n]/b_new[n]
    
    for i in 1:n-1
        
        x[n-i] = (f_new[n-i] - c[n-i]*x[n-i+1])/b_new[n-i]
    
    end
    
    return x

end



function tridiagsolve(A::Tridiagonal, f)
    
	return tridiagsolve(diag(A,-1), diag(A,0), diag(A, 1), f)
end

end

# ╔═╡ 822c7633-5c80-41ec-84f1-16f2de4d6760
let
	try
	#test2 for 1st method
	a = [8, 1, 0, 6, 0]
	b = [1, 1, 1, 1, 1, 1]
	c = [8, 0, 0, 5, 0]
	f = [1, 2, 3, 4, 5, 6]
	A = Tridiagonal(a, b, c)
	
	x = tridiagsolve(a, b, c, f)
	
	println("tridiagsol: ", x)
	println("tridiagdiff: ", A*x - f)
	println()

	
	#test2 for 2nd method
	a1 = [8, 9, 0]
	b1 = [5, 6, 7, 8]
	c1 = [1, 0, 3]
	A1 = Tridiagonal(a1, b1, c1)
	f = [1, 2, 3, 4]
	println(A1)
	x = tridiagsolve(A1, f)
		
	println("tridiagsol: ", x)
	println("tridiagdiff: ", A1*x - f)

	catch e
		println(e)
	end

end

# ╔═╡ 8c637239-f86d-4dd5-96b6-83812f330ad8
md"""
##### 3.Системы уравнений
"""

# ╔═╡ 17f94817-0459-4559-9973-d3ab90c1036b
begin
function solution_a()
    A = [
        8 9 4 1;
        0 4 1 0;
        0 0 -1 6;
        0 0 0 11 ]
    
    b = [9, 3, -1, 2]
    x = LorU_solver(A,b)
    return x, b - A * x
end

println("a: solution", solution_a())


function solution_b()
    A = [
        -2 1 0 0 0;
        1 -2 1 0 0;
        0 1 -2 1 0;
        0 0 1 -2 1;
        0 0 0 1 -2]
    
    b = [1, 1, 1, 1, 1]
    x = tridiagsolve(A,b)
    return x, b - A * x
end

println("b: solution", solution_b())


function solution_c()
    A = [
        1 8 -3 9;
        0 4 10 -2;
        8 2 -5 1;
        3 1 6 12]
    
    b = [3, 6, 1, 4]
    x = A \ b
    return x, b - A * x
end

println("c: solution", solution_c())

end

# ╔═╡ 733d0b49-ccb5-4ea9-b34b-b2cdb07fe182
let
	a = [2, 1, 1, 1, 1]
	b = [1, 5, 1, 6, 1, 1]
	c = [1, 1, 5, 1, 5]
	f = [1, 2, 3, 4, 5, 6]
	A = Tridiagonal(diagm(-1 => a, 0 => b, 1 => c))
	
	println(diag(A,-1))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
LaTeXStrings = "~1.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.1"
manifest_format = "2.0"
project_hash = "d1890f26107c4055623ee2649391bcee691e3c42"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╠═89eb2182-cf2f-432f-bb83-65569e1a7d15
# ╠═de5b06e2-78b3-11ed-33d8-7377f554fd63
# ╠═d84a0ed1-07bc-4d85-9cb4-0765c5985537
# ╠═e7e3f90a-9901-40a8-8d0f-8b02ade6abeb
# ╠═552bd38c-3637-4cca-8ee7-467117ac27f6
# ╠═3833b8dc-2469-4195-9674-401bc3b545c2
# ╠═822c7633-5c80-41ec-84f1-16f2de4d6760
# ╠═8c637239-f86d-4dd5-96b6-83812f330ad8
# ╠═17f94817-0459-4559-9973-d3ab90c1036b
# ╠═733d0b49-ccb5-4ea9-b34b-b2cdb07fe182
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
