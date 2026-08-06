using CairoMakie, LaTeXStrings

# characteristics: ((x1, t1), (x2, t2)), all of slope dt/dx = 5/4
segments = [
    ((0, 10/4), (2, 5)),
    ((0, 5/4), (3, 5)),
    ((0, 0), (4, 5)),
    ((1, 0),       (5, 5)),
    ((2, 0),       (6, 5)),
    ((3, 0),       (7, 5)),
]

fig = Figure(size = (420, 340), figure_padding = 10)
ax = Axis(fig[1, 1]; aspect = DataAspect(), limits = (-1, 5, -1, 3))
hidespines!(ax)
hidedecorations!(ax)

for ((x1, t1), (x2, t2)) in segments
    lines!(ax, [x1, x2], [t1, t2]; color = :black, linewidth = 1)
end

# axes through the origin, arrowheads as triangle markers
lines!(ax, [-1, 4.95], [0, 0]; color = :black, linewidth = 1)
lines!(ax, [0, 0], [-1, 2.95]; color = :black, linewidth = 1)
scatter!(ax, [4.95], [0]; marker = :rtriangle, color = :black, markersize = 10)
scatter!(ax, [0], [2.95]; marker = :utriangle, color = :black, markersize = 10)
text!(ax, 4.95, 0; text = L"x", align = (:right, :top), offset = (-3, -6))
text!(ax, 0, 2.95; text = L"t", align = (:right, :top), offset = (-6, -3))

# marked points
scatter!(ax, [2, 0], [0, 5/4]; color = :black, markersize = 7)
text!(ax, 2, 0;   text = L"u(x_0, 0)", align = (:center, :top),  offset = (0, -8))
text!(ax, 0, 5/4; text = L"u(0, t)",   align = (:right, :center), offset = (-8, 0))

save("chapters/assets/characteristics.svg", fig)
fig
