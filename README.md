# Appearance Mimicking Surfaces
This repo is an libigl style implementation of the [Appearance Mimicking Surfaces](https://cims.nyu.edu/gcl/papers/mimicking-2014.pdf) paper.

> **To get started:** 
> 
>     git clone --recursive http://github.com/[username]/appearance-mimicking-surfaces.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the demo from inside the `build/` by running

    ./ams [path to mesh.obj] [path to matcap image]
The example uses [MatCaps](https://libigl.github.io/tutorial/#matcaps) to style the mesh and optionally uses matlab to process long-running computations for large meshes.

## Background
See `report/entry.md`.

