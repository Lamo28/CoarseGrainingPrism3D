using LinearAlgebra
using Plots

function numchop(x::Real)
    if abs(x) <= 1000*eps()
        x = 0
    end
    return x
end
function numchop(x::Complex)
    ansreal = real(x)
    ansim = imag(x)
    if abs(ansreal) <= 1000*eps()
        ansreal = 0
    end
    if abs(ansim) <= 1000*eps()
        ansim = 0
    end
    return x= ansreal + ansim*im
end
function numchop(x::Array)
    for i in 1:length(x)
        if x[i] != 0
            x[i] = numchop(x[i])
        end
    end
    return x
end

function delta(i::Real,j::Real,k::Real,K::Int)
	sol = 0
	if i <= (j+k) && j<= (i+k) && i+j+k <= K &&2*(i+j+k)%2 == 0
		sol =1
	end
	return sol
end

function qn(n::Real,K::Int)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2)))
    return real(sol)
end

function qnfact(n::Real,K::Int)
    sol = 1
    for i in 1:n
        sol *= qn(i,K)
    end
    return sol
end

function visqrt(i::Real,K::Int)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1,K))
    return sol
end

function trian(i::Real,j::Real,k::Real,K::Int)
    sol = 0
    if delta(i,j,k,K) == 1
        sol = delta(i,j,k,K)*sqrt(qnfact(i+j-k,K)*qnfact(i-j+k,K)*qnfact(-i+j+k,K)/qnfact(i+j+k+1,K))
    end
    return sol
end

function RacWig6j(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real,K::Int)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sumz = 0
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1,K)/
                ((qnfact(e-z,K)*qnfact(f-z,K)*qnfact(g-z,K))* (qnfact(z-a,K)*qnfact(z-b,K)*qnfact(z-c,K)*qnfact(z-d,K)))
        end
        sol = trian(i,j,m,K)*trian(i,l,n,K)*trian(k,j,n,K)*trian(k,l,m,K)*sumz
    end
    return sol
end

function Fsymb(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real,K::Int)
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1,K)*qn(2*n+1,K)) * RacWig6j(i,j,m,k,l,n,K)
    end
    return sol
end

#Define G-symbol
function Gsymb(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real,K::Int)
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = Fsymb(i,j,m,k,l,n,K) /(visqrt(m,K)*visqrt(n,K))
    end
    return sol
end


function SuperIndex(AllSpin,K)
	x = K+1
	Size = length(AllSpin)
	ans = Int(2*sum(AllSpin[i]*x^(Size-i) for i in 1:Size)+1)
    return ans
end

function InverseSuperIndex(y,Size)
    yt = y - 1
    j = zeros(0)
    for i in (Size-1):-1:0
        d,r = divrem(yt,x^i)
        push!(j,d/2)
        yt = r
    end
    return j
end

maxK = 10
DataAmpK = Array{Array}(undef,maxK,1)
DataPosK = Array{Array}(undef,maxK,1)
for K in 1:1:maxK
	temp = []
	pos = []
	for i in 0:0.5:K/2, j in 0:0.5:K/2, m in 0:0.5:K/2
		if delta(i,j,m,K) !=0
			for l in 0:0.5:K/2, n in 0:0.5:K/2
				if delta(i,l,n,K) !=0
					for k in 0:0.5:K/2
						if delta(k,j,n,K) !=0 && delta(k,l,m,K) !=0
							if Gsymb(i,j,m,k,l,n,K) !=0
								push!(temp,Gsymb(i,j,m,k,l,n,K))
								push!(pos,SuperIndex([i,j,m,k,l,n],K))
							end
						end
					end
				end
			end
		end
	end
	DataAmpK[K] = temp
	DataPosK[K] = pos
end
