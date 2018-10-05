using CoarseGrainingPrism3D
println("data")
@time data = PrAOpti()
@show length(data[1])

println("time computation Prism UV and number of diff between PrismUV and initial Prism")
# @time temp = FullsvdPrismTruncated(data,1)
# AmpUV = sort(temp[1])
# AmpIni = sort(data[1])
# cnt = 0
# for i in 1:length(AmpIni)
#    if numchop(AmpIni[i]-AmpUV[i]) !=0
#        # println(AmpIni[i]-AmpUV[i])
#        global cnt +=1
#    end
# end
# println(cnt)

println("time computation PrismEff")
@time dataEff = PrismEff(data,1)
@show dataEff
# @show length(dataEff)
