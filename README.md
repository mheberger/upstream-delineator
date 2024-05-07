Notes on the Upstream Tech project.

# Timeline

Before the end of July 2024. 

# Contact

Alden Keefe Sampson, alden@upstream.tech


# Scope

Modify delineator.py scripts with 4 new features:

1. Return sub-basins upstream of the watershed outlet, rather than a single, merged polygon 
representing the entire watershed. 

2. Input will be 1 or more main outlet points, and 0 or more internal points for each outlet. 

3. Using MERIT-Basins, do the delineation and sub-basin splitting at internal points. 

4. Return a data structure which contains: 
   a) a list of sub-basin polygons 
   b) a list of sub-basin outlet points 
   c) a graph representation of which sub-basins flow into each other
   d) potentially other bonus info like river linestring geometries and the whole drainage 
   area geometry for the outlet

5. Automatically merge adjacent unit catchments to increase their size and reduce their
 number (preserving the sub-basins at internal outlets.)

6. List coordinates of sub-basin *outlets* (not centroids). 

* Not in scope for now: delineation with USGS data (NHD, HUCs). Current scope includes delineation based on MERIT-Hydro and MERIT-Basins only. 
  if Upstream desires to use USGS data in the future (it is more accurate, but only available for continental US), we can explore this fall, 
  say after September. I have experience using these data, via the NHD API and the 


# Notes

See if I can reproduce the approximate level of detail in this figure: 
https://www.upstream.tech/hydroforecast

For an outlet at Folsom Lake,  38.707, -121.157

The merging of sub-basins turned out to be a somewhat complex task. 

The river network can be represented as an "acyclic directed graph" or a "tree graph."

At each step, needs to be "pruned," using a "linear branch collapser" pruner.

https://www.researchgate.net/publication/273508304_BiNChE_A_web_tool_and_library_for_chemical_enrichment_analysis_based_on_the_ChEBI_ontology

In the end, it seems rather simple. However, it took me a fair bit of experimenting and coding to 
conclude that this was the optimal solution. 

Consider introducing Shreve stream order. Can be a useful concept. 

4/24 Vision starting to emerge on how the project can work. I don't know a priori what size basins they want. It will be good to 
have a set of options. I can create a set of successively larger basins, up 1 level, then up 2, then 3. Then the delineator routine can simply be run on those.

Next step: get the river polylines for each layer, so I can extract the basin outlets. Add this as fields in the sub-basin attribute table? 
Or I can store it in the river reaches attribute table. 


 
