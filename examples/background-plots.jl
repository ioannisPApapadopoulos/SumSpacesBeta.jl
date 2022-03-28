using Plots

xx = -4:0.001:4
y = ewU[xx,1]
p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255),dpi=1000)
savefig(p, "ewU0.png")

y = eT[xx,2]
p = plot(xx, y, grid=false, aspect_ratio=:equal, linewidth=3, color=:rainbow, line_z=y, background_color=RGB(212/255,239/255,252/255), dpi=1000)
savefig(p, "eT1.png")