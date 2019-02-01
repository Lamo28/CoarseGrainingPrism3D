module CoarseGrainingPrism3D
#####
using LinearAlgebra
using Profile
#####
#####

### GLOBAL CONSTANT
const K =3
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
function numchop(x::ComplexF64)
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
	if i <= (j+k) && j<= (i+k) && i+j+k <= K &&2*(i+j+k)%2 == 0
		sol =1
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
function Prism(β::Float64,α1::Float64,β1::Float64,α2::Float64,δ::Float64,a1::Float64,b::Float64,
                    a2::Float64,d::Float64,c1::Float64,c2::Float64,f::Float64)
    sol = 0
    if delta(c2,d,δ) !=0 && delta(c2,α2,a2) !=0 && delta(β,d,a2) !=0 && delta(β,α2,δ) !=0 &&
        delta(f,d,α1) !=0 && delta(f,c1,b) !=0 && delta(a1,d,b) !=0 && delta(a1,c1,α1) !=0 &&
        delta(α1,β1,δ) !=0 && delta(c2,β1,f) !=0

        dims = visqrt(β)*visqrt(α1)*visqrt(β1)*visqrt(α2)*visqrt(δ)*visqrt(a1)*visqrt(b)*visqrt(a2)*visqrt(d)*visqrt(c1)*visqrt(c2)*visqrt(f)
        sol = dims*Gsymb(c2,d,δ,β,α2,a2)*Gsymb(f,d,α1,a1,c1,b)*Gsymb(c2,d,δ,α1,β1,f)

    if imag(sol) == 0
        sol = real(sol)
    end
    end
    return sol
end

export SuperIndex
function SuperIndex(AllSpin)
	Size = length(AllSpin)
	ans = Int(2*sum(AllSpin[i]*x^(Size-i) for i in 1:Size)+1)
    return ans
end

export InverseSuperIndex
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

export PrAOpti
function PrAOpti()
    PrA = Array{Float64}(undef,0)
	s = Array{Int64}(undef,0)

	sol = 0
	dims = 0
    for β in 0.:0.5:y, d in 0.:0.5:y, a2 in 0.:0.5:y
        if delta(β,d,a2) !=0
            for α2 in 0.:0.5:y, δ in 0.:0.5:y
                if  delta(β,α2,δ) !=0
                    for c2 in 0.:0.5:y
                        if delta(c2,d,δ) !=0 && delta(c2,α2,a2) !=0
                            for f in 0.:0.5:y, β1 in 0.:0.5:y
                                if delta(c2,β1,f) !=0
                                    for α1 in 0.:0.5:y
                                        if delta(f,d,α1) !=0 && delta(α1,β1,δ) !=0
                                            for c1 in 0.:0.5:y, a1 in 0.:0.5:y
                                                if delta(a1,c1,α1) !=0
                                                    for b in 0.:0.5:y
                                                        if delta(f,c1,b) !=0 && delta(a1,d,b) !=0
                                                            dims = visqrt(β)*visqrt(α1)*visqrt(β1)*visqrt(α2)*visqrt(δ)*visqrt(a1)*visqrt(b)*visqrt(a2)*visqrt(d)*visqrt(c1)*visqrt(c2)*visqrt(f)
                                                            sol = numchop(dims*Gsymb(c2,d,δ,β,α2,a2)*Gsymb(f,d,α1,a1,c1,b)*Gsymb(c2,d,δ,α1,β1,f))
                                                            if sol !=0
                                                                push!(PrA,sol)
                                                                push!(s,SuperIndex([β,α1,β1,α2,δ,a1,b,a2,d,c1,c2,f]))
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


##### coarse-graining algorithm
# create function returning M_{AB}^C, where C fixed spins, A spins of tetra and B spins of pyramid
export BlocksPrism
function BlocksPrism(PrismData,c2::Float64,d::Float64,δ::Float64)
    PrA = PrismData[1]
    sIndex = PrismData[2]

    TempTetra = Array{Array}(undef,0)
    TempPyra = Array{Array}(undef,0)
    TempAmp = Array{Real}(undef,0)
    PosAmp = 1
    for i in sIndex
        (β,α1,β1,α2,δCut,a1,b,a2,dCut,c1,c2Cut,f) = InverseSuperIndex(i,12)
        if δCut == δ && dCut == d && c2Cut == c2
            push!(TempAmp,PrA[PosAmp])
            push!(TempTetra,[β α2 a2])
            push!(TempPyra,[c1 a1 α1 β1 b f])
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
function svdPrism(PrismData,c2::Float64,d::Float64,δ::Float64)# parameters c2,d,δ
    mat, TetraIndex, PyraIndex = BlocksPrism(PrismData,c2,d,δ)

    U, s, V = svd(mat)
    V = adjoint(V)

    return U, numchop(s), V, TetraIndex, PyraIndex
end

#HYP : only take the firsts TRUNC singular values
export svdPrismTruncated
function svdPrismTruncated(PrismData,c2::Float64,d::Float64,δ::Float64,trunc)
	U, s, V, TetraIndex, PyraIndex = svdPrism(PrismData,c2,d,δ)

	TruncU = U[:,1:trunc]*sqrt.(s[1:trunc])
	TruncV = transpose(sqrt.(s[1:trunc]))*V[1:trunc,:]

	if real(sqrt(visqrt(c2)*visqrt(d)*visqrt(δ))) == 0
		TruncU = im*sqrt(visqrt(c2)*visqrt(d)*visqrt(δ))*TruncU
		TruncV = -im*sqrt(visqrt(c2)*visqrt(d)*visqrt(δ))*TruncV
	else
		TruncU = sqrt(visqrt(c2)*visqrt(d)*visqrt(δ))*TruncU
		TruncV = sqrt(visqrt(c2)*visqrt(d)*visqrt(δ))*TruncV
	end
    if length(TetraIndex) != 0
        (β,α2,a2) = TetraIndex[1]
        if sign(TruncU[1,1]) != sign(AmpTetraInPrism(c2,d,δ,β,α2,a2))
             TruncU = -(1)*TruncU
             TruncV = -(1)*TruncV
        end
    end

	return real(TruncU), real(TruncV), TetraIndex, PyraIndex, s[1:trunc]
end

export FullsvdPrismTruncated
function FullsvdPrismTruncated(PrismData,trunc)
	FullU = Array{Array}(undef,0)
	FullV = Array{Array}(undef,0)
	FullTetraindex = Array{Array}(undef,0)
	FullPyraindex = Array{Array}(undef,0)
	FullCutindex = Array{Array}(undef,0)
	Fulls = Array{Array}(undef,0)
	for c2 in 0:0.5:y, d in 0:0.5:y, δ in 0:0.5:y
		if delta(c2,d,δ) == 1
			temp = svdPrismTruncated(PrismData,c2,d,δ,trunc)
            if temp[1] !=[0]
    			push!(FullU,temp[1])
    			push!(FullV,temp[2])
    			push!(FullTetraindex,temp[3])
    			push!(FullPyraindex,temp[4])
    			push!(FullCutindex,[c2,d,δ])
    			push!(Fulls,temp[5])
            end
		end
	end

    if FullU[1][1] != 1
        FullU = (1/FullU[1][1])*FullU
        FullV = (1/FullV[1][1])*FullV
    end

	return FullU, FullV, FullTetraindex, FullPyraindex, FullCutindex, Fulls
end

export PrismUVTruncated
function PrismUVTruncated(PrismData,trunc)
	FullU, FullV, FullTetraindex, FullPyraindex, FullCutindex,s = FullsvdPrismTruncated(PrismData,trunc)

	AmpUV = Array{Float64}(undef,0)
	PosUVFaceExt = Array{Int64}(undef,0)
	PosUVFaceTrue = Array{Int64}(undef,0)
	PosUVFaceFalse = Array{Int64}(undef,0)
	SizeC = length(FullCutindex)
	for i in 1:SizeC, j in 1:SizeC
		(c2,d,δ) = FullCutindex[i]
		(c2t,dt,δt) = FullCutindex[j]
		#glue along dA/dB
		if d == dt
			TempTetra = FullTetraindex[i]
			TempPyra = FullPyraindex[j]
			for iA in 1:length(TempTetra), jB in 1:length(TempPyra)
				(β,α2,a2) = TempTetra[iA]
				(c1t,a1t,α1t,β1t,bt,ft) = TempPyra[jB]
				#and glue along uA/uB and vA/vB
				if β == bt && a2 == a1t
					push!(AmpUV,FullU[i][iA]*FullV[j][jB]/(visqrt(d)*visqrt(β)*visqrt(a2)))
					push!(PosUVFaceExt,SuperIndex([β1t,ft,β]))
					push!(PosUVFaceTrue,SuperIndex([c2,δt,c2t,δ,d]))
					push!(PosUVFaceFalse,SuperIndex([c2,α1t,c1t,α2,a2]))
				end
			end
		end
	end
	return AmpUV,PosUVFaceExt,PosUVFaceTrue,PosUVFaceFalse, s
end


export PrismEff
function PrismEff(PrismData,trunc)

	AmpUV,PosUVFaceExt,PosUVFaceTrue,PosUVFaceFalse, s = PrismUVTruncated(PrismData,trunc)
	FullAmp = Array{Real}(undef,0)
	BigSIndex = Array{Real}(undef,0)
	for i in 1:length(PosUVFaceTrue)
		# spins of first prism
		(Χ,c1,a2,c2,d) = InverseSuperIndex(PosUVFaceTrue[i],5) #shared spins
		(γ,e,f) = InverseSuperIndex(PosUVFaceExt[i],3)
		(χ,α1,a1,α2,β) = InverseSuperIndex(PosUVFaceFalse[i],5) # was face 2

		pos = findall(x->x == PosUVFaceTrue[i],PosUVFaceFalse)
		for j in pos
			# spins of second prism
			#(bA,j1A,j4B,j3B,uA)
			(γt,et,ft) = InverseSuperIndex(PosUVFaceExt[j],3)
			(Χ,c1t,a2t,c2t,dt) = InverseSuperIndex(PosUVFaceTrue[j],5)
			#super index for 18 spins prism
			TempIndex = SuperIndex([χ,α1,a1,α2,β,c1t,a2t,c2t,dt,γ,e,f,γt,et,ft,c1,c2,a2])
			# FullTot = vcat(FullTot,[TempIndex AmpUV[i]*AmpUV[j]])
			push!(FullAmp,AmpUV[i]*AmpUV[j])
			push!(BigSIndex,TempIndex)
		end
	end

    p = sortperm(BigSIndex)
    BigSIndex = BigSIndex[p]
    FullAmp = FullAmp[p]

    NewSindex = Array{Int64}(undef,0)
    NewAmp = Array{Float64}(undef,0)

    cnt = 0
    for i in 2:length(BigSIndex)
        cnt +=1
        if BigSIndex[i] == BigSIndex[i-1] && i != length(BigSIndex)
            FullAmp[i-cnt] += FullAmp[i]
        else
            push!(NewAmp,FullAmp[i-cnt])
            push!(NewSindex,BigSIndex[i-1])
            cnt = 0
        end
    end
    if BigSIndex[end] == BigSIndex[end-1]
        NewAmp[end] += FullAmp[end]
    else
        push!(NewAmp,FullAmp[end])
        push!(NewSindex,BigSIndex[end])
    end
    return NewAmp, NewSindex
end

####
"embedding map"
####

# 2-2 move for diag
# WARNING: Move done with Fsym, wich need to be either redefine at every step ? or always keep
# initial one ? but then we know something weird might happen -> do test !! IN THE FOLLOWING
# DONE WITH INITIAL Fsymb

export Move2_2
function Move2_2(SetAmpAndIndex,SizeSetSpins,PosMove,PosTrianglesOfMove)
	# SetIndexAndAmp: full set of amplitude and SuperIndex pos of non zero
	# SizeSetSpins: number of edges in polyhedra associated to SetIndexAndAmp
	# PosMove: Position of the edge/spin that we want to flip
	# PosTrianglesOfMove: Positions of the OTHER 4 edges, forming with PosMove the two ini triangles
	# 				that is PosTrianglesOfMove = [j1,j2 , h1,h2] and j1j2Posmove / h1h2Posmove ini triangles
	#				AND j1,h1 / j2,h2 share a vertex
	Amp,sIndex = SetAmpAndIndex

	TempAmp = Array{Real}(undef,0)
	TempsIndex = Array{Real}(undef,0)
	for i in 1:length(sIndex)
		SetSpins = InverseSuperIndex(sIndex[i],SizeSetSpins)
		for qMove in 0:0.5:y
            TempFsymb = Fsymb(SetSpins[PosTrianglesOfMove[1]],SetSpins[PosTrianglesOfMove[2]],SetSpins[PosMove],SetSpins[PosTrianglesOfMove[4]],SetSpins[PosTrianglesOfMove[3]],qMove)
            if numchop(TempFsymb) != 0
                TempSetSpins = copy(SetSpins)
				push!(TempsIndex,SuperIndex(insert!(deleteat!(TempSetSpins,PosMove),PosMove,qMove)))
				push!(TempAmp,TempFsymb*Amp[i])
			end
		end
	end

    p = sortperm(TempsIndex)
    TempsIndex = TempsIndex[p]
    TempAmp = TempAmp[p]

    NewsIndex = Array{Int64}(undef,0)
	NewData = Array{Float64}(undef,0)

	cnt = 0
	for i in 2:length(TempsIndex)
        cnt +=1
		if TempsIndex[i] == TempsIndex[i-1] && i != length(TempsIndex)
			TempAmp[i-cnt] += TempAmp[i]
		else
            if numchop(TempAmp[i-cnt]) != 0
                push!(NewData,TempAmp[i-cnt])
                push!(NewsIndex,TempsIndex[i-1])
            end
		    cnt = 0
		end
	end
    if TempsIndex[end] == TempsIndex[end-1]
        NewData[end] += TempAmp[end]
    else
        push!(NewData,TempAmp[end])
        push!(NewsIndex,TempsIndex[end])
    end

	return NewData,NewsIndex
end


export Move3_1
function Move3_1(SetAmpAndIndex,SizeSetSpins,PosTrianglesOfMove)
    # SetIndexAndAmp: full set of amplitude and SuperIndex pos of non zero
	# SizeSetSpins: number of edges in polyhedra associated to SetIndexAndAmp
	# PosTrianglesOfMove: Position of the 6 edges of the move, of the form
    #                   [j1 j2 j3 h1 h2 h3] where ji are disappearing spins, h1 h2 h3 keeping spins
    #                                triangles are: j1j3h1, j1j2h2, j2j3h3
	Amp,sIndex = SetAmpAndIndex

    PositionDelete = Array{Bool}(undef,0)
    for i in 1:SizeSetSpins
        if i in PosTrianglesOfMove[1:3]
            push!(PositionDelete,true)
        else
            push!(PositionDelete,false)
        end
    end

    TempAmp = Array{Real}(undef,0)
	TempsIndex = Array{Real}(undef,0)
    for i in 1:length(sIndex)
        SetSpins = InverseSuperIndex(sIndex[i],SizeSetSpins)
        TempDimFactor = (visqrt(SetSpins[PosTrianglesOfMove[1]])*visqrt(SetSpins[PosTrianglesOfMove[3]]))/visqrt(SetSpins[PosTrianglesOfMove[4]])
        TempFsymb =Fsymb(SetSpins[PosTrianglesOfMove[6]],SetSpins[PosTrianglesOfMove[5]],SetSpins[PosTrianglesOfMove[4]],SetSpins[PosTrianglesOfMove[1]],SetSpins[PosTrianglesOfMove[3]],SetSpins[PosTrianglesOfMove[2]])
        if numchop(TempFsymb) != 0
            TempSetSpins = copy(SetSpins)
            push!(TempsIndex,SuperIndex(deleteat!(TempSetSpins,PositionDelete)))
            push!(TempAmp,(1/(TempDimFactor*TempFsymb))*Amp[i])
        end
    end

   p = sortperm(TempsIndex)
   TempsIndex = TempsIndex[p]
   TempAmp = TempAmp[p]

   NewsIndex = Array{Int64}(undef,0)
   NewData = Array{Real}(undef,0)
   cntlast = 0
   cnt = 0
   for i in 2:length(TempsIndex)
       cnt+=1
       if TempsIndex[i] == TempsIndex[i-1] && i != length(TempsIndex)
           #TempAmp[i-cnt] += TempAmp[i]
       else
           if numchop(TempAmp[i-cnt]) != 0
               #push!(NewData,(1/cnt)*TempAmp[i-cnt])
               push!(NewData,TempAmp[i-cnt])
               push!(NewsIndex,TempsIndex[i-1])
               if i == length(TempsIndex)
                   cntlast=cnt
               end
           end
           cnt = 0
        end
    end
    if TempsIndex[end] == TempsIndex[end-1]
        #NewData[end] = cntlast*NewData[end]
        #NewData[end] += TempAmp[end]
        #NewData[end] = (1/(cntlast+1))*NewData[end]
    else
        push!(NewData,TempAmp[end])
        push!(NewsIndex,TempsIndex[end])
    end

    if NewData[1] != 1
        NewData = (1/NewData[1])*NewData
    end

	return NewData,NewsIndex
end


#### Embedding of PrismEff to space of Prism again
export EmbPrismEff
function EmbPrismEff(PrismData,trunc)
    SetAmpAndIndex = PrismEff(PrismData,trunc)

    #Move2_2 for big diagonal
    #Move2_2_BigDiag = Move2_2(SetAmpAndIndex,18,18,[10,11,14,15])
    Move2_2FirstTriangle = Move2_2(SetAmpAndIndex,18,16,[2,10,6,13])

    # Move 2_2 first triangle
    #Move2_2FirstTriangle = Move2_2(Move2_2_BigDiag,18,16,[2,10,6,13])
    Move2_2_BigDiag = Move2_2(Move2_2FirstTriangle,18,18,[10,11,14,15])

    # Move3_1 first triangle
    #Move3_1FirstTriangle = Move3_1(Move2_2FirstTriangle,18,[10,13,14,18,16,7])
    Move3_1FirstTriangle = Move3_1(Move2_2_BigDiag,18,[10,13,14,18,16,7])

    # Move 2_2 second triangle: SHOULD OR NOT ?????
    Move2_2SecondTriangle = Move2_2(Move3_1FirstTriangle,15,14,[4,11,8,12])

    # Move3_1 second triangle: SHOULD OR NOT ?????
    Move3_1SecondTriangle = Move3_1(Move2_2SecondTriangle,15,[10,11,12,15,3,14])

    EffPrism = Move3_1SecondTriangle

    return EffPrism
end


##### test function



#### data from initial tetra in prism
export AmpTetraInPrism
function AmpTetraInPrism(j4::Real,j5::Real,j6::Real,j1::Real,j2::Real,j3::Real)
	ans = 0
	if delta(j4,j5,j6)==1 && delta(j4,j2,j3)==1 && delta(j1,j5,j3)==1 && delta(j1,j2,j6)==1
		dims = visqrt(j1)*visqrt(j2)*visqrt(j3)*visqrt(j4)*visqrt(j5)*visqrt(j6)
		ans = dims*Gsymb(j4,j5,j6,j1,j2,j3)
	end
	return ans
end

export FullAmpTetraInPrism
function FullAmpTetraInPrism()
    Amp = Array{Real}(undef,0)
    Index = Array{Real}(undef,0)
    for j1 in 0:0.5:y, j2 in 0:0.5:y, j3 in 0:0.5:y, j4 in 0:0.5:y, j5 in 0:0.5:y, j6 in 0:0.5:y
        TempAmp = numchop(AmpTetraInPrism(j1,j2,j3,j4,j5,j6))
        if TempAmp!=0
            push!(Amp,TempAmp)
            push!(Index,SuperIndex([j1,j2,j3,j4,j5,j6]))
        end
    end
    return Amp,Index
end

export FullAmpPyramidInPrism
function FullAmpPyramidInPrism()
	Amp = Array{Real}(undef,0)
	Index = Array{Real}(undef,0)
	for j1 in 0:0.5:y, j2 in 0:0.5:y, j3 in 0:0.5:y, j4 in 0:0.5:y, j5 in 0:0.5:y, j6 in 0:0.5:y,
		a in 0:0.5:y, b in 0:0.5:y, c in 0:0.5:y
		TempAmp1 = numchop(AmpTetraInPrism(a,b,c,j1,j2,j3))
		TempAmp2 = numchop(AmpTetraInPrism(a,b,c,j4,j5,j6))
		if TempAmp1!=0 && TempAmp2!=0
			push!(Amp,TempAmp1*TempAmp2/(visqrt(a)*visqrt(b)*visqrt(c)))
			push!(Index,SuperIndex([a,b,c,j1,j2,j3,j4,j5,j6]))
		end
	end
	return Amp,Index
end



end # module
