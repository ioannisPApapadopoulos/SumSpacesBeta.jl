using SumSpaces
using Plots
using DelimitedFiles
using LaTeXStrings

errors = readdlm("uS-saved/errors-inf.txt") 
N = Int32.(errors[:,2])
errors = errors[:,1]


plot(N, log10.(errors),
title=L"$\mathrm{Error}$",
legend=false,
# ylabel=L"$\Vert u - \mathbf{S}^{\mathbf{I},\!\!\!\!+}_n\!\!\!\!\! \mathbf{u}\Vert_\infty$"
ylabel=L"$\log_{10}(\infty\mathrm{-norm \;\; error})$",
xlabel=L"$n$",
ylim=[log10(1e-15),log10(1e-3)],
xtickfontsize=10, ytickfontsize=10,xlabelfontsize=15,ylabelfontsize=15,
color=:black,
linewidth=2
)

savefig("frac-helmholtz-error.pdf")