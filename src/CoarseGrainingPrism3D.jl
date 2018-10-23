module CoarseGrainingPrism3D
#####
using LinearAlgebra

#####
#####

### GLOBAL CONSTANT
const K = 2
const x = K+1
const y = K/2

### FUNDAMENTAL FUNCTION
#chop function
export numchop
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

#Define delta_{ijkl} -> coupling rules
export delta
function delta(i::Real,j::Real,k::Real)
    sol = 0
    if i <=(j+k) && j<=(i+k) && k<=(i+j) && i+j+k <= K && 2*(i+j+k)%2 ==0
        sol = 1
    end
    return sol
end

#Define quantum numbers qn (this is real)
export qn
function qn(n::Real)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2)))
    return real(sol)
end

#Define qn factorial
export qnfact
function qnfact(n::Real)
    sol = 1
    for i in 1:n
        sol *= qn(i)
    end
    return sol
end

#Define square root of quantum dimension
export visqrt
function visqrt(i::Real)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1))
    return sol
end

#Define triangle equality
export trian
function trian(i::Real,j::Real,k::Real)
    sol = 0
    if delta(i,j,k) == 1
        sol = delta(i,j,k)*sqrt(qnfact(i+j-k)*qnfact(i-j+k)*qnfact(-i+j+k)/qnfact(i+j+k+1))
    end
    return sol
end

# Define Racah-Wigner six-j symbol
export RacWig6j
function RacWig6j(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sumz = 0
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1)/
                ((qnfact(e-z)*qnfact(f-z)*qnfact(g-z))* (qnfact(z-a)*qnfact(z-b)*qnfact(z-c)*qnfact(z-d)))
        end
        sol = trian(i,j,m)*trian(i,l,n)*trian(k,j,n)*trian(k,l,m)*sumz
    end
    return sol
end

#Define F-symbol
export Fsymb
function Fsymb(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1)*qn(2*n+1)) * RacWig6j(i,j,m,k,l,n)
    end
    return sol
end

#Define G-symbol
export Gsymb
function Gsymb(i::Real,j::Real,m::Real,k::Real,l::Real,n::Real)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = Fsymb(i,j,m,k,l,n) /(visqrt(m)*visqrt(n))
    end
    return sol
end

### DATA
# Define amplitude for a prism function of 12 edges : usual 9 + 3 diag
export Prism
function Prism(ja1::Float64,jb1::Float64,jd1::Float64,je1::Float64,jf1::Float64,ja2::Float64,jb2::Float64,
                jd2::Float64,je2::Float64,jf2::Float64,jc::Float64,jg::Float64)
    sol = 0
    if delta(jd1,jd2,je1) != 0 && delta(jd1,jb1,jg) != 0 && delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0 && delta(jd1,jc,jb2) != 0 && delta(jf1,jd2,jb2) != 0 && delta(jf1,jc,je1) != 0 && delta(jd1,ja2,jf2) != 0 && delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0
        dims = visqrt(ja1)*visqrt(jb1)*visqrt(jg)*visqrt(jf1)*visqrt(jc)*visqrt(jb2)*visqrt(je2)*visqrt(ja2)*visqrt(jf2)*visqrt(jd1)*visqrt(jd2)*visqrt(je1)
        sol =  dims*Gsymb(jd1,jd2,je1,ja1,jb1,jg) * Gsymb(jd1,jd2,je1,jf1,jc,jb2) * Gsymb(jd1,jc,jb2,je2,ja2,jf2)
    end
    return sol
end

function SuperIndex(AllSpin)
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
#actual data of the non zero amplitude for Prism store as two vectors: amp value and their positions
export PrA
function PrA()
    PrA = zeros(0)
    s = zeros(0)
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y,
                jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        ans = numchop( Prism(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
        if ans != 0
            push!(PrA,ans)
            push!(s,SuperIndex([ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]))
        end
    end
    return PrA,s
end

export PrAOpti
function PrAOpti()
	PrA = zeros(0)
	s = zeros(0)

	sol = 0
	dims = 0
	for jd1 in 0:0.5:y, jd2 in 0:0.5:y, je1  in 0:0.5:y
		if delta(jd1,jd2,je1) != 0
			for jb1 in 0:0.5:y, jg in 0:0.5:y
				if delta(jd1,jb1,jg) != 0
					for ja1 in 0:0.5:y
						if delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0
							for jc in 0:0.5:y, jb2 in 0:0.5:y
								if delta(jd1,jc,jb2) != 0
									for jf1 in 0:0.5:y
										if delta(jf1,jd2,jb2) != 0  && delta(jf1,jc,je1) != 0
											for ja2 in 0:0.5:y, jf2 in 0:0.5:y
												if delta(jd1,ja2,jf2) != 0
													for je2 in 0:0.5:y
														if delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0
															dims = visqrt(ja1)*visqrt(jb1)*visqrt(jg)*visqrt(jf1)*visqrt(jc)*visqrt(jb2)*visqrt(je2)*visqrt(ja2)*visqrt(jf2)*visqrt(jd1)*visqrt(jd2)*visqrt(je1)
													        sol =  numchop(dims*Gsymb(jd1,jd2,je1,ja1,jb1,jg) * Gsymb(jd1,jd2,je1,jf1,jc,jb2) * Gsymb(jd1,jc,jb2,je2,ja2,jf2))
															if sol !=0
																push!(PrA,sol)
												            	push!(s,SuperIndex([ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]))
															end
														end
													end
												end
											end
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
	return PrA, s
end


#OUTDATED prismB : Prism with modified diag
function prismB(ja1::Float64,jb1::Float64,jd1p::Float64,je1::Float64,jf1::Float64,ja2::Float64,jb2::Float64,
		                 jd2p::Float64,je2::Float64,jf2::Float64,jc::Float64,jg::Float64)
    sol = 0
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y
        sol += numchop(Fsymb(ja1,jg,jd2,jb2,jf1,jd2p)*Fsymb(jb1,jg,jd1,ja2,jf2,jd1p)*prismA(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
    end
    return sol
end
export PrB
function PrB()
    PrA = zeros(0)
    s = zeros(0)
    for ja1 in 0.:0.5:y, jb1 in 0.:0.5:y, jd1 in 0.:0.5:y, je1 in 0.:0.5:y, jf1 in 0.:0.5:y, ja2 in 0.:0.5:y,
                jb2 in 0.:0.5:y, jd2 in 0.:0.5:y, je2 in 0.:0.5:y, jf2 in 0.:0.5:y, jc in 0.:0.5:y, jg in 0.:0.5:y
        ans = numchop( prismB(ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg))
        if ans != 0
            push!(PrA,ans)
            push!(s,SuperIndex([ja1,jb1,jd1,je1,jf1,ja2,jb2,jd2,je2,jf2,jc,jg]))
        end
    end
    return PrA,s
end

##### coarse-graining algorithm
# create function returning M_{AB}^C, where C fixed spins, A spins of tetra and B spins of pyramid
export BlocksPrism
function BlocksPrism(PrismData,jd1::Real,je1::Real,jd2::Real)
    PrA = PrismData[1]
    sIndex = PrismData[2]

    TempTetra = Array{Array}(undef,0)
    TempPyra = Array{Array}(undef,0)
    TempAmp = Array{Real}(undef,0)
    PosAmp = 1
    for i in sIndex
        (ja1,jb1,jd1t,je1t,jf1,ja2,jb2,jd2t,je2,jf2,jc,jg) = InverseSuperIndex(i,12)
        if jd1t == jd1 && je1t == je1 && jd2t == jd2
            push!(TempAmp,PrA[PosAmp])
            push!(TempTetra,[ja1 jb1 jg])
            push!(TempPyra,[jf1 ja2 jb2 je2 jf2 jc])
        end
        PosAmp += 1
    end
    PosTetra = unique(TempTetra)
    PosPyra = unique(TempPyra)
    SizeA = length(PosTetra)
    SizeB = length(PosPyra)
    mat = zeros(SizeA,SizeB)
    for i in 1:length(TempTetra)
        qa = findall(x -> x == TempTetra[i], PosTetra)
        qb = findall(x -> x == TempPyra[i], PosPyra)
        mat[qa[1],qb[1]] = TempAmp[i]
    end

    if SizeA == 0
        mat = 0
    end
    return mat, PosTetra, PosPyra
end

export svdPrism
function svdPrism(PrismData,jd1::Real,je1::Real,jd2::Real)# parameters d1,d2,e1
    mat, Aindex, Bindex = BlocksPrism(PrismData,jd1,je1,jd2)

    U, s, V = svd(mat)
    V = adjoint(V)

    return U, numchop(s), V, Aindex, Bindex
end

#HYP : only take the firsts TRUNC singular values
export svdPrismTruncated
function svdPrismTruncated(PrismData,jd1::Real,je1::Real,jd2::Real,trunc)
	U, s, V, Aindex, Bindex = svdPrism(PrismData,jd1::Real,je1::Real,jd2::Real)

	TruncU = U[:,1:trunc]*sqrt.(s[1:trunc])
	TruncV = transpose(sqrt.(s[1:trunc]))*V[1:trunc,:]

	return TruncU, TruncV, Aindex, Bindex, s[1:trunc]
end

export FullsvdPrismTruncated
function FullsvdPrismTruncated(PrismData,trunc)
	FullU = Array{Array}(undef,0)
	FullV = Array{Array}(undef,0)
	FullAindex = Array{Array}(undef,0)
	FullBindex = Array{Array}(undef,0)
	FullCindex = Array{Array}(undef,0)
	Fulls = Array{Array}(undef,0)
	for jd1 in 0:0.5:y, je1 in 0:0.5:y, jd2 in 0:0.5:y
		if delta(jd1,je1,jd2) == 1
			temp = svdPrismTruncated(PrismData,jd1,je1,jd2,trunc)
			push!(FullU,temp[1])
			push!(FullV,temp[2])
			push!(FullAindex,temp[3])
			push!(FullBindex,temp[4])
			push!(FullCindex,[jd1,je1,jd2])
			push!(Fulls,temp[5])
		end
	end
	return FullU, FullV, FullAindex, FullBindex, FullCindex, Fulls
end

export PrismUVTruncated
function PrismUVTruncated(PrismData,trunc)
	FullU, FullV, FullAindex, FullBindex, FullCindex,s = FullsvdPrismTruncated(PrismData,trunc)

	AmpUV = Array{Real}(undef,0)
	PosUV12 = Array{Real}(undef,0)
	PosUVFace1 = Array{Real}(undef,0)
	PosUVFace2 = Array{Real}(undef,0)
	SizeC = length(FullCindex)
	for i in 1:SizeC, j in 1:SizeC
		(dA,bA,aA) = FullCindex[i]
		(dB,bB,aB) = FullCindex[j]
		if dA == dB
			TempA = FullAindex[i]
			TempB = FullBindex[j]
			for iA in 1:length(TempA), jB in 1:length(TempB)
				(j1A,uA,vA) = TempA[iA]
				(j2B,uB,j3B,j4B,vB,j5B) = TempB[jB]
				if uA == uB && vA == vB
					push!(AmpUV,FullU[i][iA]*FullV[j][jB])
					push!(PosUV12,SuperIndex([vA,j2B,j5B]))
					push!(PosUVFace1,SuperIndex([bA,j1A,j4B,j3B,uA]))
					push!(PosUVFace2,SuperIndex([bA,aA,bB,aB,dA]))
				end
			end
		end
	end
	return AmpUV,PosUV12,PosUVFace1,PosUVFace2, s
end

export PrismEff
function PrismEff(PrismData,trunc)
	AmpUV,PosUV12,PosUVFace1,PosUVFace2,s = PrismUVTruncated(PrismData,trunc)
	FullAmp = Array{Real}(undef,0)
	BigSIndex = Array{Real}(undef,0)
	for i in 1:length(PosUVFace1)
		# spins of first prism
		(bA,j1A,j4B,j3B,uA) = InverseSuperIndex(PosUVFace1[i],5) #shared spins
		(vA,j2B,j5B) = InverseSuperIndex(PosUV12[i],3)
		(bA,aA,bB,aB,dA) = InverseSuperIndex(PosUVFace2[i],5)

		pos = findall(x->x == PosUVFace1[i],PosUVFace2)
		for j in pos
			# spins of second prism
			#(bA,j1A,j4B,j3B,uA)
			(vAt,j2Bt,j5Bt) = InverseSuperIndex(PosUV12[j],3)
			(bA,j1At,j4Bt,j3Bt,uAt) = InverseSuperIndex(PosUVFace1[j],5)
			#super index for 18 spins prism
			TempIndex = SuperIndex([bA,aA,bB,aB,dA,j1At,j4Bt,j3Bt,uAt,j1A,j4B,j3B,vA,vAt,j2B,j2Bt,j5B,j5Bt])
			# FullTot = vcat(FullTot,[TempIndex AmpUV[i]*AmpUV[j]])
			push!(FullAmp,AmpUV[i]*AmpUV[j])
			push!(BigSIndex,TempIndex)
		end
	end
	FullTot = hcat(BigSIndex,FullAmp)
	p = sortperm(FullTot[:,1])
	FullTot = FullTot[p,:]
	NewFullTot = zeros(length(unique(FullTot[:,1])),2)
	pos = 1
	cnt = 1
	for i in 2:size(FullTot,1)
		if FullTot[i,1] == FullTot[i-1,1]
			FullTot[i-cnt,2] += FullTot[i,2]
			cnt += 1
		else
			NewFullTot[pos,:] = FullTot[i-cnt,:]
			cnt = 1
			pos +=1
		end
	end
	return NewFullTot
end

export EmbeddingMap1
function EmbeddingMap1(PrismData,trunc,truncEmbedding1)
	FullAmp = PrismEff(PrismData,trunc)
	Amp = FullAmp[:,2]
	sIndex = FullAmp[:,1]

	#first embedding map, SVD on j1
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in sIndex
		(β,α1,α2,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,18)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j1)
		push!(TempRest,SuperIndex([β,α1,α2,d,α1t,α2t,α3t,dt,q,j2,u,ut,v,vt,q1,q1t]))
		pos +=1
	end

	PosLine = unique(TempLine)
	PosRest = unique(TempRest)
	SizeL = length(PosLine)
	SizeR = length(PosRest)
	mat = zeros(SizeR,SizeL)
	for i in 1:length(TempLine)
        qb = findall(x -> x == TempLine[i], PosLine)
        qa = findall(x -> x == TempRest[i], PosRest)
        mat[qa[1],qb[1]] = TempAmp[i]
    end
    if SizeL == 0
        mat = 0
    end

	U, s, V = svd(mat)
	V = adjoint(V)
	sEmb1 = numchop(s)
	AmpEmb1 = U[:,1:truncEmbedding1]*sqrt.(sEmb1[1:truncEmbedding1])
	TruncVEmb1 = transpose(sqrt.(sEmb1[1:truncEmbedding1]))*V[1:truncEmbedding1,:]

	# second embedding, on j2,  should be equivalent to the first one. multiplying by the V
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in PosRest
		(β,α1,α2,d,α1t,α2t,α3t,dt,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,17)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j2)
		push!(TempRest,SuperIndex([β,α1,α2,d,α1t,α2t,α3t,dt,q,u,ut,v,vt,q1,q1t]))
		pos +=1
	end

	PosLine = unique(TempLine)
	PosRest = unique(TempRest)
	SizeL = length(PosLine)
	SizeR = length(PosRest)
	mat = zeros(SizeR,SizeL)
	for i in 1:length(TempLine)
		qb = findall(x -> x == TempLine[i], PosLine)
		qa = findall(x -> x == TempRest[i], PosRest)
		mat[qa[1],qb[1]] = TempAmp[i]
	end
	if SizeL == 0
		mat = 0
	end

	U, s, V = svd(mat)
	V = adjoint(V)
	sEmb2 = numchop(s)
	AmpEmb2 = U[:,1:truncEmbedding1]*sqrt.(sEmb2[1:truncEmbedding1])
	TruncVEmb2 = transpose(sqrt.(sEmb2[1:truncEmbedding1]))*V[1:truncEmbedding1,:]


	return AmpEmb2, PosRest, TruncVEmb1,TruncVEmb2, sEmb1,sEmb2

end

export EmbeddingMap1Commutator
function EmbeddingMap1Commutator(PrismData,trunc,truncEmbedding1)
	FullAmp = PrismEff(PrismData,trunc)
	Amp = FullAmp[:,2]
	sIndex = FullAmp[:,1]

	#first embedding map, SVD on j1
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in sIndex
		(β,α1,α2,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,18)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j2)
		push!(TempRest,SuperIndex([β,α1,α2,d,α1t,α2t,α3t,dt,j1,q,u,ut,v,vt,q1,q1t]))
		pos +=1
	end

	PosLine = unique(TempLine)
	PosRest = unique(TempRest)
	SizeL = length(PosLine)
	SizeR = length(PosRest)
	mat = zeros(SizeR,SizeL)
	for i in 1:length(TempLine)
        qb = findall(x -> x == TempLine[i], PosLine)
        qa = findall(x -> x == TempRest[i], PosRest)
        mat[qa[1],qb[1]] = TempAmp[i]
    end
    if SizeL == 0
        mat = 0
    end

	U, s, V = svd(mat)
	V = adjoint(V)
	sEmb1 = numchop(s)
	AmpEmb1 = U[:,1:truncEmbedding1]*sqrt.(sEmb1[1:truncEmbedding1])
	TruncVEmb1 = transpose(sqrt.(sEmb1[1:truncEmbedding1]))*V[1:truncEmbedding1,:]

	# second embedding, on j2,  should be equivalent to the first one. multiplying by the V
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in PosRest
		(β,α1,α2,d,α1t,α2t,α3t,dt,j1,q,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,17)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j1)
		push!(TempRest,SuperIndex([β,α1,α2,d,α1t,α2t,α3t,dt,q,u,ut,v,vt,q1,q1t]))
		pos +=1
	end

	PosLine = unique(TempLine)
	PosRest = unique(TempRest)
	SizeL = length(PosLine)
	SizeR = length(PosRest)
	mat = zeros(SizeR,SizeL)
	for i in 1:length(TempLine)
		qb = findall(x -> x == TempLine[i], PosLine)
		qa = findall(x -> x == TempRest[i], PosRest)
		mat[qa[1],qb[1]] = TempAmp[i]
	end
	if SizeL == 0
		mat = 0
	end

	U, s, V = svd(mat)
	V = adjoint(V)
	sEmb2 = numchop(s)
	AmpEmb2 = U[:,1:truncEmbedding1]*sqrt.(sEmb2[1:truncEmbedding1])
	TruncVEmb2 = transpose(sqrt.(sEmb2[1:truncEmbedding1]))*V[1:truncEmbedding1,:]


	return AmpEmb2, PosRest, TruncVEmb1,TruncVEmb2, sEmb1,sEmb2

end


end # module
