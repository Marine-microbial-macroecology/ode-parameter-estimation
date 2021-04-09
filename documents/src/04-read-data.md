# Phytoplankton culture data

Here I will show culture data for N-starved cultures of phytoplankton, tracked over about 10 days with DIN (nitrate + nitrite, ammonia), PON (cell quota), and cell numbers, plus estimated growth rate and many other parameters.

```@example Ex1
using DataFrames, CSV, StatsPlots, Query
import Gadfly
```

Read the data from a file. Simplify some of the names

```@example Ex1
liefer = CSV.File("liefer-growth-data.csv") |> DataFrame ; 
# names(liefer)
rename!(liefer, [2 => :days, 5 => :cell_density, 7 => :dilution_factor])
nothing # hide
```

Cultures were diluted after time 0 according to the data in "dilution". Use this factor to scale the cell density and DIN at time 0.

```@example Ex1
liefer |> @group_by(:Species) |> @map(max(_.dilution_factor)) # max doesn't work here?
```

We will use the following columns

* Species
* Replicate
* days (time since start of each experiment, in days; sometimes missing)
* cell_density (cell L-1)
* DIN (µmol L-1)
* N (pg cell-1)

Make some plots to show the data.

```@example Ex1
@df liefer scatter(:days, :N, group = :Species)
liefer |> @filter(_.Species == "Micromonas sp.")  |> @df scatter(:days, :N, group = :Replicate)
@df liefer scatter(:days, log.(:cell_density), group = :Species)
```

Using Gadfly

```@example Ex1
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "N", Gadfly.Geom.subplot_grid(Gadfly.Geom.point))
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "X", Gadfly.Geom.subplot_grid(Gadfly.Geom.point))
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "cell_density", Gadfly.Geom.subplot_grid(Gadfly.Geom.point, free_y_axis=true), 
    Gadfly.Scale.y_log10)
```

