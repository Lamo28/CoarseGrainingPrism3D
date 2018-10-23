using CoarseGrainingPrism3D
println("data")
@time data = PrAOpti()
@show length(data[1])

@time temp = EmbeddingMap1(data,1,1)
@time tempC = EmbeddingMap1Commutator(data,1,1)

println("1 then 2")

@show temp[5]
@show temp[6]
@show length(temp[1])

println("2 then 1")
@show tempC[5]
@show tempC[6]
@show length(tempC[1])

println("diff")
@show sum(numchop(sort(temp[1]) - sort(tempC[1])).!=0)

# println("time computation PrismEff")
# @time dataEff = PrismEff(data,1)
# @show length(dataEff[:,1])
# cnt = 0
# for i in numchop(dataEff[:,2])
#    if i != 0
#       global cnt +=1
#    end
# end
# println(cnt)
