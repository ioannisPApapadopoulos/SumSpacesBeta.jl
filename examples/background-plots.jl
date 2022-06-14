using SumSpaces
using Plots

xx = -4:0.001:4
ewU = ExtendedWeightedChebyshevU()
eT = ExtendedChebyshevT()
## Rainbow
# y = ewU[xx,1]
# p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255),dpi=1000)
# savefig(p, "ewU0.png")

# y = eT[xx,2]
# p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255), dpi=1000)
# savefig(p, "eT1.png")

## Imperial blue
icblue = RGB(0/255,62/255,116/255)
lightgrey = RGB(235/255,238/255,238/255)
coolgrey = RGB(157/255,157/255,157/255)
lightblue = RGB(212/255,239/255,252/255)

for j = 1:20
    y = ewU[xx,j]
    p = plot(xx, y, grid=false, aspect_ratio=:equal, 
                linewidth=12, color=coolgrey, 
                background_color=lightblue,
                legend=false,
                axis=([],false)) #line_z=y,
    savefig(p, "ewU$(j-1).pdf")
end

for j = 1:20
    y = eT[xx,j]
    p = plot(xx, y, grid=false, aspect_ratio=:equal, 
                linewidth=12, color=coolgrey, 
                background_color=lightblue,
                legend=false,
                axis=([],false))
    savefig(p, "eT$(j-1).pdf") #line_z=y, 
end

##################

ewT = ExtendedWeightedChebyshevT()
eT = ExtendedChebyshevT()
ewU = ExtendedWeightedChebyshevU()
eU = ExtendedChebyshevU()

xx = -8:0.01:8
y = ewT[xx, 2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black)
savefig(p, "ewT1.pdf")

xx = -3:0.01:3
y = eT[xx,2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black, size=(1000,300))
savefig(p, "eT1.pdf")

xx = -3:0.01:3
y = ewU[xx, 1]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black)
savefig(p, "ewU0-2.pdf")

xx = -8:0.01:8
y = eU[xx,3]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black, size=(800,300))
savefig(p, "eU0-2.pdf")