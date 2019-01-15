using CoarseGrainingPrism3D
println("data")
@time data = PrAOpti()
@show length(data[1])

# @time FullsvdPrismTruncated(data,1)
# @time tempPrismUV = PrismUVTruncated(data,1)
# @show length(tempPrismUV[1])
#
# @show sum(numchop((sort(data[1])-sort(tempPrismUV[1]))).==0)

# println("prism eff")
# @time dataEff = PrismEff(data,1)
# @show length(dataEff[1])
println("test 2-2 move on initial prism")
println("one time")
@time Move2_2(data,12,12,[10,7,3,11])
#using Profile
#@profile Move2_2(data,12,12,[10,7,3,11])
#Profile.print()
@time dataIniPrismMove2_2 = Move2_2(data,12,12,[10,7,3,11])
@show sum(numchop(sort(data[1])-sort(dataIniPrismMove2_2[1])).==0)
@time dataIniPrismMove2_2_2nd_possibility = Move2_2(data,12,5,[2,3,1,4])

@time dataIniPrismMove2_2_second = Move2_2(dataIniPrismMove2_2,12,12,[3,10,11,7])

@show length(dataIniPrismMove2_2[1])
@show length(dataIniPrismMove2_2_2nd_possibility[1])
@show length(dataIniPrismMove2_2_second[1])


println("3-1 move test")
println("stupid test, all spins to 0")
@show Move3_1([1,1],6,[1 2 3 4 5 6])

println("time test for Prism Eff")
@time dataeff = PrismEff(data,1)
@show length(dataeff[1])

println("time test for move")
@time Move2_2_BigDiag = Move2_2(dataeff,18,17,[10,11,14,15])
@time Move2_2FirstTriangle = Move2_2(Move2_2_BigDiag,18,16,[2,10,6,13])
@time Move3_1FirstTriangle = Move3_1(Move2_2FirstTriangle,18,[10,13,14,17,16,7])
@show Move3_1FirstTriangle
#@time Move2_2SecondTriangle = Move2_2(Move3_1FirstTriangle,15,15,[4,11,8,12])
#@time Move3_1SecondTriangle = Move3_1(Move2_2SecondTriangle,15,[10,11,12,14,3,15])



println("test Embedding")
#@show length(unique(DataEmbPrismEff[2]))
#@show length(DataEmbPrismEff[1])

#@show (DataEmbPrismEff[1][1],DataEmbPrismEff[2][1])

#println("comparing after embedding and initial prism")
#@show sum(numchop(sort(DataEmbPrismEff[1])-sort(data[1])).==0)
#@show sum((sort(DataEmbPrismEff[2])-sort(data[2])).==0)

#@show sort(DataEmbPrismEff[1])
# @show sort(data[1])
# println("2-2 move for big diag")

# @time dataEffMove2_2 = Move2_2(dataEff,18,11,[14,18,17,15])
# @show length(dataEffMove2_2[1])
# @show sum(dataEffMove2_2[1].!=0)

# println("amp of tetrahedra in Prism")
# @time dataSVD = FullsvdPrismTruncated(data,1)
# @show sum(length(dataSVD[1][i]) for i in 1:length(dataSVD[1]))
# dataTetraFromPrism = Array{Real}(undef,0)
# for i in 1:length(dataSVD[1])
#     for j in dataSVD[1][i]
#         if j!=0
#             push!(dataTetraFromPrism,j)
#         end
#     end
# end
# dataPyramidFromPrism = Array{Real}(undef,0)
# for i in 1:length(dataSVD[2])
#     for j in dataSVD[2][i]
#         if j!=0
#             push!(dataPyramidFromPrism,j)
#         end
#     end
# end
#
#
# println("function for AmpTetraInPrism")
# @time dataTetraInPrism = FullAmpTetraInPrism()
# @show length(dataTetraInPrism[1])
#
# println("difference between tetrahedron from prism and teterahedron")
# DiffTetra = numchop((sort(dataTetraFromPrism) - sort(dataTetraInPrism[1])))
# @show sum(DiffTetra.!=0)
# for i in 1:length(DiffTetra)
#     if DiffTetra[i] !=0
#         println(sort(dataTetraFromPrism)[i],"\t",sort(dataTetraInPrism[1])[i])
#     end
# end
#
#
# println("function for Amp of pyramid")
# @time dataPyramidInPrism = FullAmpPyramidInPrism()
# @show length(dataPyramidInPrism[1])
# @show length(dataPyramidFromPrism)
#
# println("difference between pyramid from prism and pyramid")
# DiffPyramid = numchop(sort(dataPyramidFromPrism) - sort(dataPyramidInPrism[1]))
# # for i in 1:length(DiffPyramid)
# #     if DiffPyramid[i] != 0
# #         println(sort(dataPyramidFromPrism)[i] ,"\t" ,sort(dataPyramidInPrism[1])[i])
# #     end
# # end
# @show sum(DiffPyramid.!=0)
#
# println("prismUV and its diff with PrA")
# @time dataUV = PrismUVTruncated(data,1)
# @show sum(numchop(sort(dataUV[1]) - sort(data[1])).!=0)

# println("prism eff")
# @time dataEff = PrismEff(data,1)
#
# println("2-2 move for big diag")






# println("amp of tetrahedra in Prism")
# @time dataSVD = FullsvdPrismTruncated(data,1)
# @show sum(length(dataSVD[1][i]) for i in 1:length(dataSVD[1]))
# dataTetraFromPrism = Array{Real}(undef,0)
# for i in 1:length(dataSVD[1])
#     for j in dataSVD[1][i]
#         push!(dataTetraFromPrism,j)
#     end
# end
#
# println("function for AmpTetraInPrism")
# @time dataTetraInPrism = FullAmpTetraInPrism()
# @show length(dataTetraInPrism[1])
#
# println("is there a diff ???")
# @show sum(numchop((sort(dataTetraFromPrism) - sort(dataTetraInPrism[1]))).!=0)
#
# println("is spin config the same for would-be symmetric edges ?")
# Spin = Array{Real}(undef,0)
# Spint = Array{Real}(undef,0)
# for i in dataEff[:,1]
#     spins = InverseSuperIndex(i,18)
#     push!(Spin,spins[13])
#     push!(Spint,spins[14])
# end
#
# @show sum((sort(Spin)-sort(Spint)).!=0)
#
# println("embedding stuff")
# @time temp = EmbeddingMap1(data,1,1)
# @time tempC = EmbeddingMap1Commutator(data,1,1)
#
# println("1 then 2")
#
# @show temp[5]
# @show temp[6]
# @show length(temp[1])
#
# println("2 then 1")
# @show tempC[5]
# @show tempC[6]
# @show length(tempC[1])
#
# println("diff")
# @show sum(numchop(sort(temp[1]) - sort(tempC[1])).!=0)
# println("amp of tetrahedra in Prism")
# @time dataSVD = FullsvdPrismTruncated(data,1)
# @show sum(length(dataSVD[1][i]) for i in 1:length(dataSVD[1]))
# dataTetraFromPrism = Array{Real}(undef,0)
# for i in 1:length(dataSVD[1])
#     for j in dataSVD[1][i]
#         push!(dataTetraFromPrism,j)
#     end
# end
#
# println("function for AmpTetraInPrism")
# @time dataTetraInPrism = FullAmpTetraInPrism()
# @show length(dataTetraInPrism[1])
#
# println("is there a diff ???")
# @show sum(numchop((sort(dataTetraFromPrism) - sort(dataTetraInPrism[1]))).!=0)
#
