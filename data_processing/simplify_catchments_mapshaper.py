"""
I discovered that my "simplified" versions of the unit catchments had an issue, 
namely that some small features had been removed during the simplification. 
(I had used QGIS and the Simplify Tool, having selected Visvalingam.)

Alternatively, mapshaper has an option that it will never remove features. 
I had previously only used this via the web interface, thinking it would be
complicated to install, since it is written in javascript and requires node.js,
which I had never used. 

It turned out to be pretty straightforward to install and use. 

However, I could not get mapshaper to run from a batch file. I don't know why. 
Rather than spending time trying to figure it out, I just used this script, 
which uses Python's `subprocess` library. This worked like a charm.

"""

import subprocess


def simplify_shapefile(input_shapefile, output_shapefile, simplify_percentage):
    # Command to run Mapshaper to simplify the shapefile
    command = f"mapshaper -i {input_shapefile} -simplify {simplify_percentage}% keep-shapes -o format=shapefile {output_shapefile}"
    print(command)
    # Run the command using subprocess
    subprocess.run(command, shell=True)


basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

for basin in basins:

    input_shapefile  = rf"C:\Data\GIS\MERITBasins\catchments\src\cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"
    output_shapefile = rf"C:\Data\GIS\MERITBasins\catchments\tmp\cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"
    
    simplify_percentage = 22  # I found that 22% gave reasonable results, and dramatically reduces the file sizes. 

    simplify_shapefile(input_shapefile, output_shapefile, simplify_percentage)
