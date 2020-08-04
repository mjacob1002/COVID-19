using DataFrames, CSV

county_data = CSV.read(joinpath(Pkg.dir("DataFrames"),
"/Users/evaxueyaoguo/Desktop/NYT_us-counties.csv"));

group_by_county = groupby(county_data, :fips)

#group_by_state_and_county = groupby(county_data, [:state, :county])
