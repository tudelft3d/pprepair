## What is pprepair?

pprepair (*p*lanar *p*artition repair) takes a set of polygons and ensures that they form a valid planar partition, made of valid polygons and having no gaps or overlaps. It can be used as a validator, telling of problems in individual polygons or in the planar partition, and also as an automatic repair tool, outputting a set of polygons that do form a valid planar partition. If you are only interested in repairing individual polygons, have a look at [prepair](https://github.com/tudelft3d/prepair).

## What is a planar partition and why is it important to have valid planar partitions?

Planar partitions are subdivisions of the plane into polygons, and are frequently used in GIS to model various concepts like land use, geology, administrative subdivisions, natural features and cadastral parcels, among many others.

However, the polygons in a planar partition are often created separately, and thus cannot be expected to fit with each other exactly. Other times, these polygons are stored and modified separately, causing different errors and inconsistencies to be introduced. These come in the form of invalid polygons, gaps, overlaps and disconnected polygons.

When software that expects a planar partition received one that is not so, it can give erroneous results (in the best case), or fail to give a result at all, often without a clear explanation to the user.

## How does pprepair work?

In short, pprepair creates a constrained triangulation of the polygons, tags each triangle with the polygon that it belongs to, modifies the triangulation to ensure that only one tag is present in each, and reconstructs the polygons from the triangulation. 

Many more details are available in Ken Arroyo Ohori's MSc thesis [here](https://3d.bk.tudelft.nl/ken/files/10_msc_thesis.pdf).

If you use pprepair for your research, please cite this publication:

> Arroyo Ohori, Ken, Ledoux, Hugo and Meijers, Martijn (2012). Validation and automatic repair of planar partitions using a constrained triangulation. *Photogrammetrie, Fernerkundung, Geoinformation (PFG)*. 5:613–630. [ [PDF] ](https://3d.bk.tudelft.nl/ken/files/12_pfg.pdf) [ [DOI] ](http://dx.doi.org/10.1127/1432-8364/2012/0143)

## How do I use pprepair?

pprepair is a command-line program, which we provide as source code. It compiles on Linux and macOS; the build system also supports Windows with MSVC.

To compile pprepair, you first need to install the free libraries [CGAL](http://www.cgal.org) (version 6 or newer) and [GDAL](http://www.gdal.org) (version 3 or newer), as well as [CMake](http://www.cmake.org) 3.20 or newer. Under macOS, you can install pprepair directly using [Homebrew](http://brew.sh/) by typing `brew install tudelft3d/software/pprepair`.

Then, from the root of the source tree, run an out-of-source build:

    $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    $ cmake --build build --parallel
    $ ctest --test-dir build

On Windows, run the same commands from a Developer Command Prompt. With a multi-config generator such as Visual Studio, use `cmake --build build --config Release` instead.

After successful compilation, run:

    $ ./build/pprepair -i inputfile -o outputfile -fix

You can get all the (advanced) options by running pprepair with no arguments, and the version with `--version`:

    $ ./build/pprepair --version
    $ ./build/pprepair

Input may be any vector format supported by GDAL. The output driver is selected from the file extension; supported extensions are `.csv`, `.dxf`, `.gdb`, `.geojson`/`.json`, `.gml`, `.gpkg`, `.kml`, `.parquet` and `.shp`. GeoParquet support requires GDAL's Parquet driver, which is provided by the `libgdal-arrow-parquet` package on conda-forge.

## Can I use pprepair in my organisation? What is the license of pprepair?

pprepair is available under the [GPLv3](http://www.gnu.org/copyleft/gpl.html) licence, which allows you to use, copy and modify the software freely. However, if you incorporate pprepair in your software, you must distribute the source code of your software, as well as any modifications made to pprepair, under the GPLv3 as well.

## Help! pprepair is crashing or giving weird results

This can be due to several reasons. 99% of the time it can be solved by:
  - If your data set has very large polygons, try passing the -bd flag
  - If your data set has points that are **very** close together, try uncommenting line 26 in definitions/CGALDefinitions.h (#define EXACT_CONSTRUCTIONS)

If your problem persists, please report it [here](https://github.com/tudelft3d/pprepair/issues?state=open).
