# Usage

## Unzip STL file
$ gzip -d constant/triSurface/puyo_LOW.stl.gz  

## Source OpenFOAM v1706 etc/bashrc and generate mesh
$ source ~/OpenFOAM/OpenFOAM-v1706/etc/bashrc  
$ mkdir constant/polyMesh  
$ cp blockMeshDict constant/polyMesh  
$ blockMesh  
$ snappyHexMesh -overwrite  

## Source foam-extend-4.0 (hf-foam-4.0 merged) etc/bashrc and run the case
$ source ~/foam/foam-extend-4.0/etc/bashrc  
$ ./Allrun  

## Cleanup
$ ./Allclean

# Original STL Data
http://honda-3d.com/
