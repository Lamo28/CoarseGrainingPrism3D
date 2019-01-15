using CoarseGrainingPrism3D
println("data")
@time data = PrAOpti()
@show length(data[1])

println("prism eff")
@time dataEff = PrismEff(data,1)
@show length(dataEff[1])

println("test 2-2 move on initial prism")
println("one time")
@time Move2_2(data,12,12,[10,7,3,11])
#using Profile
#@profile Move2_2(data,12,12,[10,7,3,11])
#Profile.print()
@time dataIniPrismMove2_2 = Move2_2(data,12,12,[10,7,3,11])
@show length(dataIniPrismMove2_2[2])
# @show sum(numchop(sort(data[1])-sort(dataIniPrismMove2_2[1])).==0)

# println("3-1 move test")
# println("stupid test, all spins to 0")
# @show Move3_1([1,1],6,[1 2 3 4 5 6])
#
# println("time test for move")
# @time Move2_2_BigDiag = Move2_2(dataEff,18,17,[10,11,14,15])
# @time Move2_2FirstTriangle = Move2_2(Move2_2_BigDiag,18,16,[2,10,6,13])
# @time Move3_1FirstTriangle = Move3_1(Move2_2FirstTriangle,18,[10,13,14,17,16,7])
# @time Move2_2SecondTriangle = Move2_2(Move3_1FirstTriangle,15,15,[4,11,8,12])
# @time Move3_1SecondTriangle = Move3_1(Move2_2SecondTriangle,15,[10,11,12,14,3,15])
#
# println("test Embedding")
# @time DataEmbPrismEff = EmbPrismEff(data,1)
# @show length(unique(DataEmbPrismEff[2]))
# @show length(DataEmbPrismEff[1])
#
# println("comparing after embedding and initial prism")
# @show sum(numchop(sort(DataEmbPrismEff[1])-sort(data[1])).==0)
# @show sum((sort(DataEmbPrismEff[2])-sort(data[2])).==0)
