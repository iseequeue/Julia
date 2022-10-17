### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 0b3448c0-4e66-11ed-285d-0bf6234f2318
begin
	using LinearAlgebra

function romberg(f, a, b; atol=1e-6, maxstep::Integer=100)
    maxstep = max(1, maxstep)  # хотя бы одно разбиение
    I = Matrix{Float64}(undef, maxstep+1, maxstep+1)
    I[1, 1] = (b - a) * (f(a) + f(b)) / 2
    for i in 2:maxstep+1
        let hc = (b - a) / 2^(i-1), np = 2^(i-2)
            I[i, 1] = I[i-1, 1] / 2 + hc * sum(f, (a + hc * (2i-1) for i in 1:np))
        end
        for k in i-1:-1:1
            I[k, i-k+1] = (2^i*I[k+1, i-k] - I[k, i-k]) / (2^i - 1)
        end
        abs(I[1, i] - I[2, i-1]) < atol && return I[1, i]
    end
    error("Точность не удовлетворена.")
end


	foo(x) = sin(100*x*exp(-x^2))
	a, b, atol = -1/3, 3, 1e-10
	println("5.5.1 I = ", romberg(foo, a, b; atol=atol))




	function intadapt(f, a, b, tol, xtol=eps(), fa=f(a), fb=f(b), m=(b-a)/2, fm=f(m))
    if a > b; a, b = b, a; end
    
    xl = (a + m)/2; fl = f(xl)  # расположение:
    xr = (m + b)/2; fr = f(xr)  # a -- xl -- m -- xr -- b
    
    T = Vector{Float64}(undef, 3)
    h = b - a
    T[1] = h * (fa + fb)/2
    T[2] = T[1]/2 + h/2 * fm
    T[3] = T[2]/2 + h/4 * (fl + fr)
    S = (4*T[2:end] - T[1:2]) / 3

    err = (S[2] - S[1]) / 15
    
    if abs(err) < tol * (1 + tol * abs(S[2]))
        Q = S[2]
        nodes = [a, xl, m, xr, b]
    else
        b - a ≤ xtol && error("Достигнут предел точности отрезка интегрирования `xtol`.")
        Ql, nodesl = intadapt(f, a, m, tol, xtol, fa, fm, xl, fl)
        Qr, nodesr = intadapt(f, m, b, tol, xtol, fm, fb, xr, fr)
        Q = Ql + Qr
        nodes = [nodesl; nodesr[2:end]]
    end
    return (Q, nodes)
end
	
	acc = (exp(pi/2)-1)/2
	boo(x) = exp(x)*cos(x)
	a, b = 0, pi/2

	for i in 2:12

		Q, nodes = intadapt(boo, a, b, 1/10^(i))
		diff = abs(Q-acc)

		println(diff, ' ', length(nodes))
		
	end

	
	

end

# ╔═╡ 6919ed24-00c3-412f-99b7-2f85764abaab


# ╔═╡ d155a5aa-fd97-4291-ba2a-8624ecc9ae77


# ╔═╡ 563c35a4-71b8-4da2-86aa-b677fc2cdfb0


# ╔═╡ a4752dd0-4cdf-4c44-b3b2-94a4349b16dd


# ╔═╡ 2de55366-878d-4325-9642-a2e20200eb03


# ╔═╡ d34a832c-0df1-404e-9c6b-721d517dc689


# ╔═╡ 050c2356-98a6-4bde-b65f-0e193245bd25


# ╔═╡ f9cc5c47-f5d2-4483-9bfe-f5ef1f32b5d9


# ╔═╡ cde05094-7aeb-4d90-894a-301509faa61c


# ╔═╡ 443b7cfd-c231-47e5-8236-2eecd77e4c76


# ╔═╡ 64869539-3903-4505-9181-bad170c2a8e1


# ╔═╡ b9aecebc-16a5-4eb9-b757-1d9c20f0e119


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.1"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

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
# ╠═0b3448c0-4e66-11ed-285d-0bf6234f2318
# ╠═6919ed24-00c3-412f-99b7-2f85764abaab
# ╠═d155a5aa-fd97-4291-ba2a-8624ecc9ae77
# ╠═563c35a4-71b8-4da2-86aa-b677fc2cdfb0
# ╠═a4752dd0-4cdf-4c44-b3b2-94a4349b16dd
# ╠═2de55366-878d-4325-9642-a2e20200eb03
# ╠═d34a832c-0df1-404e-9c6b-721d517dc689
# ╠═050c2356-98a6-4bde-b65f-0e193245bd25
# ╠═f9cc5c47-f5d2-4483-9bfe-f5ef1f32b5d9
# ╠═cde05094-7aeb-4d90-894a-301509faa61c
# ╠═443b7cfd-c231-47e5-8236-2eecd77e4c76
# ╠═64869539-3903-4505-9181-bad170c2a8e1
# ╠═b9aecebc-16a5-4eb9-b757-1d9c20f0e119
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
