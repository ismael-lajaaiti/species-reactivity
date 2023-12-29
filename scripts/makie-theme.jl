import Makie: texfont

p_theme = Theme(;
    # fonts = Attributes(
    #     :bold => texfont(:bold),
    #     :bolditalic => texfont(:bolditalic),
    #     :italic => texfont(:italic),
    #     :regular => texfont(:regular),
    # ),
    Axis = Attributes(;
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticklabelsize = 10,
        yticklabelsize = 10,
    ),
    Colorbar = Attributes(; ticklabelsize = 10),
    Legend = Attributes(;
        backgroundcolor = :grey80,
        rowgap = -5,
        labelsize = 10,
        titlesize = 10,
    ),
)

