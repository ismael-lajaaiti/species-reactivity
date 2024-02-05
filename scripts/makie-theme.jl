import ColorSchemes: Accent_3

# Theme for publications.
publication_theme = Theme(;
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

# Useful variables for plotting.
letters = string.(collect('A':'Z')) # To label figure panels.
palette = ColorSchemes.Accent_8

# Save figure with physical dimension.
cm_to_pt = 28.3465 # 1 cm = 28.3465 pt, cf. CairoMakie documentation.
single_column_width = 8.2 # In centimeners, cf. Ecology Letter guidelines.
two_third_page_width = 11
full_page_width = 17.3
width_height_ratio = 1.5
