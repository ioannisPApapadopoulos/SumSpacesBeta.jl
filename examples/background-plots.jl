using SumSpaces
using Plots

xx = -4:0.001:4
y = ewU[xx,1]
p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255),dpi=1000)
savefig(p, "ewU0.png")

y = eT[xx,2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255), dpi=1000)
savefig(p, "eT1.png")


ewT = ExtendedWeightedChebyshevT()
eT = ExtendedChebyshevT()

xx = -8:0.01:8
y = ewT[xx, 2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black)
savefig(p, "ewT1.pdf")

xx = -3:0.01:3
y = eT[xx,2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, legend=false, axis=([],false), linewidth=4, color=:black, size=(1000,300))
savefig(p, "eT1.pdf")