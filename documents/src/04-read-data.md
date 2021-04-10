# Phytoplankton culture data

Here I will show culture data for N-starved cultures of phytoplankton, tracked over about 10 days with DIN (nitrate + nitrite, ammonia), PON (cell quota), and cell numbers, plus estimated growth rate and many other parameters.

```@example Ex1
using DataFrames, CSV, StatsPlots, Query
import Gadfly
```

Read the data from a file. Simplify some of the names. Change missing dilution factors to 1, and compute cumulative dilution factor by multiplication.

```@example Ex1
liefer = CSV.File("liefer-growth-data.csv") |> DataFrame ; 
# names(liefer)
rename!(liefer, [2 => :days, 5 => :cell_density, 7 => :dilution_factor])
replace!(liefer.dilution_factor, missing => 1.0)
nothing # hide
```

Make a grouped data view of the dataframe and work on that to manipulate one species and replicate at a time. Cultures were diluted after time 0 according to the data in "dilution". Use this factor to scale the cell density and DIN to simulate the values at time 0, just after dilution.

```@example Ex1
# liefer = @pipe groupby(liefer, [:Species, :Replicate]) |> transform(_, :dilution_factor => cumprod)
# liefer = @pipe groupby(liefer, [:Species, :Replicate]) |> 
#               transform(_, [:cell_density, :dilution_factor_cumprod ] => ./ , renamecols = false)
# rename!(liefer, 44 => :cells)
#rename!(liefer, Dict(:cell_density_dilution_factor_cumprod => "cells"))
# liefer.cells
# liefer = @pipe groupby(liefer, [:Species, :Replicate]) |> 
#                transform(_, [:DIN, :dilution_factor_cumprod ] => ./ , renamecols = false)
# rename!(liefer, Dict(:DIN_dilution_factor_cumprod => "DIN_scaled"))
gdf = groupby(liefer, [:Species, :Replicate])
transform!(gdf, :dilution_factor => cumprod)
transform!(gdf, :dilution_factor_cumprod => ( (x) -> x ./ minimum(x)) => :dilution_factor_min)
transform!(gdf, [:cell_density, :dilution_factor_min] => ( ./ ) => :cells)
```

At time 0, the cultures were washed to remove most of the nutrients, so we don't have a good estimate of DIN at time 0. We can approximate it by computing total mass at the second sample point, or just convert them to missing values.

```@example Ex1
# transform!(gdf, [:DIN, :dilution_factor_cumprod] => ( .* ) => :DIN_scaled)
transform!(liefer, [:days, :DIN] => ByRow((x,y) -> x === 0 ? missing : y) => :DIN_corr)
transform!(liefer, [:DIN_corr, :N, :cells] => ( (R,Q,X) -> R .+ Q .* X ) => :mass)
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
liefer |> @filter(_.Species == "Micromonas sp.")  |> @df scatter(:days, :DIN_corr, group = :Replicate)
liefer |> @filter(_.Species == "Micromonas sp.")  |> @df scatter(:days, :mass, group = :Replicate)
@df liefer scatter(:days, log.(:cells), group = :Species)
```

Using Gadfly

```@example Ex1
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "N", Gadfly.Geom.subplot_grid(Gadfly.Geom.point,
    free_y_axis = true))
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "DIN_corr", Gadfly.Geom.subplot_grid(Gadfly.Geom.point,
    free_y_axis = true))
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "cells", Gadfly.Geom.subplot_grid(Gadfly.Geom.point,
    free_y_axis=true), 
    Gadfly.Scale.y_log10)
Gadfly.plot(liefer, 
    ygroup = "Species", x = "days", y = "mass", Gadfly.Geom.subplot_grid(Gadfly.Geom.point,
    free_y_axis = true))
```

Create subsets of the data, one for each species, and just the variables t, R, Q, X to be used in the 
model in the next section.

