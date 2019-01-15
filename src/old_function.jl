export BlocksPrism
function BlocksPrism(PrismData,c2::Float64,d::Float64,δ::Float64)
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

#### SVD embedding map, by extracting the edges we don't want to keep

export EmbeddingMap1
function EmbeddingMap1(PrismData,trunc,truncEmbedding1)
	sIndex,Amp = PrismEff(PrismData,trunc)

	#first embedding map, SVD on j1
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in sIndex
		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,18)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j1)
		push!(TempRest,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,q,j2,u,ut,v,vt,q1,q1t]))
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
		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,17)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j2)
		push!(TempRest,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,q,u,ut,v,vt,q1,q1t]))
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
	sIndex,Amp = PrismEff(PrismData,trunc)

	#first embedding map, SVD on j1
	TempLine = Array{Real}(undef,0)
	TempRest = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)
	pos = 1
	for i in sIndex
		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,18)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j2)
		push!(TempRest,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,u,ut,v,vt,q1,q1t]))
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
		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,17)
		push!(TempAmp,Amp[pos])
		push!(TempLine,j1)
		push!(TempRest,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,q,u,ut,v,vt,q1,q1t]))
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



function Move2_2diag(PrismData,trunc)
	sIndex,Amp = PrismEff(PrismData,trunc)

	NewData = Array{Real}(undef,0)
	NewsIndex = Array{Real}(undef,0)

	TempIndexWithoutqWithqMove = Array{Real}(undef,0)
	TempAmp = Array{Real}(undef,0)

	for i in 1:length(sIndex)
		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(sIndex[i],18)
		for qMove in 0:0.5:y
			if delta(qMove,q1,ut)==1 && delta(qMove,q1t,v)==1 && delta(qMove,α3,dt) ==1
				push!(TempIndexWithoutqWithqMove,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,j2,u,ut,v,vt,q1,q1t,qMove]))
				push!(TempAmp,Fsymb(ut,q1t,q,v,q1,qMove)*Amp[i])
			end
		end
	end

	UniqueTempIndexWithoutqWithqMove = unique(TempIndexWithoutqWithqMove)
	for i in 1:length(UniqueTempIndexWithoutqWithqMove)
		pos = findall(j-> j == UniqueTempIndexWithoutqWithqMove[i],TempIndexWithoutqWithqMove)
		push!(NewData,sum(TempAmp[pos]))
	end

#	p = sortperm(vec(TempIndexWithoutqWithqMove))
#	TempIndexWithoutqWithqMove = TempIndexWithoutqWithqMove[p]
#	TempPos = TempPos[p]
#
#	UniqueTempIndexWithoutqWithqMove = unique(TempIndexWithoutqWithqMove)
#	CountUniqueTempIndexWithoutqWithqMove = Array{Real}(undef,length(UniqueTempIndexWithoutqWithqMove))
#	for i in 1:length(UniqueTempIndexWithoutqWithqMove)
#		CountUniqueTempIndexWithoutqWithqMove[i] = count(j->j == UniqueTempIndexWithoutqWithqMove[i],TempIndexWithoutqWithqMove)
#	end

	# Comparaison = 0
	# for i in TempIndexWithoutqWithqMove
	# 	if Comparaison == i
	# 		(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,j2,u,ut,v,vt,q1,q1t,qMove) = InverseSuperIndex(i,18)
	# 		sol += Fsymb(ut,q1t,qa,v,q1,qMove)*Amp[j]
	# 	end
	# 	Comparaison = i
	# end
	#
	# for i in UniqueTempIndexWithoutqWithqMove
	# 	sol = 0
	# 	cnt = 1
	# 	if TempIndexWithoutqWithqMove[cnt] == i
	# 		sol += Fsymb(ut,q1t,qa,v,q1,qMove)*Amp[j]


	# for qMove in 0:0.5:y, i in sIndex
	# 	(β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,q,j2,u,ut,v,vt,q1,q1t) = InverseSuperIndex(i,18)
	# 	if delta(qMove,q1,ut)==1 && delta(qMove,q1t,v)==1 && delta(qMove,α3,dt) ==1
	# 		sol = 0
	# 		for j in 1:length(sIndex)
	# 			(βa,α1a,α2a,α3a,da,α1ta,α2ta,α3ta,dta,j1a,qa,j2a,ua,uta,va,vta,q1a,q1ta) = InverseSuperIndex(sIndex[j],18)
	# 			if β==βa && α1a==α1 && α2a==α2 && α3==α3a && da==d && α1ta==α1t && α2ta==α2t && α3ta==α3ta && dta==dt && j1a==j1 && j2a==j2 && ua==u && uta==ut && va==v && vta==vt && q1a==q1 && q1ta==q1t
	# 				sol += Fsymb(ut,q1t,qa,v,q1,qMove)*Amp[j]
	# 			end
	# 		end
	# 		push!(NewData,sol)
	# 		push!(NewsIndex,SuperIndex([β,α1,α2,α3,d,α1t,α2t,α3t,dt,j1,qMove,j2,u,ut,v,vt,q1,q1t]))
	# 	end
	# end

	return NewData,UniqueTempIndexWithoutqWithqMove
end


export PrismUVTruncated
function PrismUVTruncated(PrismData,trunc)
	FullU, FullV, FullTetraindex, FullPyraindex, FullCindex,s = FullsvdPrismTruncated(PrismData,trunc)

	AmpUV = Array{Float64}(undef,0)
	PosUVFaceExt = Array{Int64}(undef,0)
	PosUVFaceTrue = Array{Int64}(undef,0)
	PosUVFaceFalse = Array{Int64}(undef,0)
	SizeC = length(FullCindex)
	for i in 1:SizeC, j in 1:SizeC
		(dA,bA,aA) = FullCindex[i]
		(dB,bB,aB) = FullCindex[j]
		#glue along dA/dB
		if dA == dB
			TempA = FullTetraindex[i]
			TempB = FullPyraindex[j]
			for iA in 1:length(TempA), jB in 1:length(TempB)
				(j1A,uA,vA) = TempA[iA]
				(j2B,uB,j3B,j4B,vB,j5B) = TempB[jB]
				#and glue along uA/uB and vA/vB
				if uA == uB && vA == vB
					push!(AmpUV,FullU[i][iA]*FullV[j][jB]/(visqrt(dA)*visqrt(uA)*visqrt(vB)))
					push!(PosUV12,SuperIndex([vA,j2B,j5B]))
					push!(PosUVFace1,SuperIndex([bA,j1A,j4B,j3B,uA]))
					push!(PosUVFace2,SuperIndex([bA,aA,bB,aB,dA]))
				end
			end
		end
	end
	return AmpUV,PosUV12,PosUVFace1,PosUVFace2, s
end
