using CoarseGrainingPrism3D
println("data")
@time data = PrAOpti()
@show length(data[1])

println("amp of tetrahedra in Prism")
@time dataSVD = FullsvdPrismTruncated(data,1)
@show sum(length(dataSVD[1][i]) for i in 1:length(dataSVD[1]))
dataTetraFromPrism = Array{Real}(undef,0)
for i in 1:length(dataSVD[1])
    for j in dataSVD[1][i]
        if j!=0
            push!(dataTetraFromPrism,j)
        end
    end
end
dataPyramidFromPrism = Array{Real}(undef,0)
for i in 1:length(dataSVD[2])
    for j in dataSVD[2][i]
        if j!=0
            push!(dataPyramidFromPrism,j)
        end
    end
end
@show length(dataPyramidFromPrism)

println("function for AmpTetraInPrism")
@time dataTetraInPrism = FullAmpTetraInPrism()
@show length(dataTetraInPrism[1])

println("difference between tetrahedron from prism and teterahedron")
@show sum(numchop((sort(dataTetraFromPrism) - sort(dataTetraInPrism[1]))).!=0)

println("function for Amp of pyramid")
@time dataPyramidInPrism = FullAmpPyramidInPrism()
@show length(dataPyramidInPrism[1])

println("difference between pyramid from prism and pyramid")
@show sum(numchop((sort(real(dataPyramidFromPrism)) - sort(real(dataPyramidInPrism[1])))).!=0)

println("prismUV and its diff with PrA")
@time dataUV = PrismUVTruncated(data,1)
@show length(dataUV[1])
@show sum(numchop(sort(dataUV[1]) - sort(data[1])).!=0)

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
