package require vtk
package require vtkinteraction
package require vtktesting

# Demonstrate how to use structured grid blanking with an image. There are two
# techniques demonstrated: one uses an image to perform the blanking;
# the other uses scalar values to do the same thing. Both images should
# be identical.
#

# get the interactor ui

# create pipeline - start by extracting a single plane from the grid
#

vtkMultiBlockPLOT3DReader pl3d
pl3d SetXYZFileName "$VTK_DATA_ROOT/Data/combxyz.bin"
pl3d SetQFileName "$VTK_DATA_ROOT/Data/combq.bin"
pl3d SetScalarFunctionNumber 100
pl3d SetVectorFunctionNumber 202
pl3d Update
set output [[pl3d GetOutput] GetBlock 0]

vtkExtractGrid plane
plane SetInputData $output
plane SetVOI 0 57 0 33 0 0
plane Update

# Create some data to use for the (image) blanking
#
vtkImageData blankImage

#vtkType.h has definitions for vtk datatypes VTK_INT, VTK_FLOAT, etc. that
#don't get wrapped in Tcl.
set VTK_UNSIGNED_CHAR 3

blankImage SetDimensions 57 33 1
blankImage AllocateScalars $VTK_UNSIGNED_CHAR 1
[[blankImage GetPointData] GetScalars] SetName "blankScalars"

set blanking [[blankImage GetPointData] GetScalars]
set numBlanks [expr 57*33]
for {set i 0} {$i<$numBlanks} {incr i} {
   $blanking SetComponent $i 0 1
}

# Manually blank out areas corresponding to dilution holes
$blanking SetComponent 318 0 0
$blanking SetComponent 945 0 0
$blanking SetComponent 1572 0 0
$blanking SetComponent 641 0 0
$blanking SetComponent 1553 0 0

# The first blanking technique uses the image to set the blanking values
#
vtkBlankStructuredGridWithImage blankIt
blankIt SetInputConnection [plane GetOutputPort]
blankIt SetBlankingInputData blankImage

vtkStructuredGridGeometryFilter blankedPlane
blankedPlane SetInputConnection [blankIt GetOutputPort]
blankedPlane SetExtent 0 100 0 100 0 0

vtkPolyDataMapper planeMapper
planeMapper SetInputConnection [blankedPlane GetOutputPort]
planeMapper SetScalarRange 0.197813 0.710419

vtkActor planeActor
planeActor SetMapper planeMapper

# The second blanking technique uses grid data values to create the blanking.
# Here we borrow the image data and threshold on that.
#
vtkStructuredGrid anotherGrid
anotherGrid CopyStructure [plane GetOutput]
[anotherGrid GetPointData] SetScalars [[blankImage GetPointData] GetScalars]

vtkBlankStructuredGrid blankGrid
blankGrid SetInputData anotherGrid
blankGrid SetArrayName "blankScalars"
blankGrid SetMinBlankingValue -0.5
blankGrid SetMaxBlankingValue  0.5

vtkStructuredGridGeometryFilter blankedPlane2
blankedPlane2 SetInputConnection [blankGrid GetOutputPort]
blankedPlane2 SetExtent 0 100 0 100 0 0

vtkPolyDataMapper planeMapper2
planeMapper2 SetInputConnection [blankedPlane2 GetOutputPort]
planeMapper2 SetScalarRange 0.197813 0.710419

vtkActor planeActor2
planeActor2 SetMapper planeMapper2

# An outline around the data
#
vtkStructuredGridOutlineFilter outline
outline SetInputData $output

vtkPolyDataMapper outlineMapper
outlineMapper SetInputConnection [outline GetOutputPort]

vtkActor outlineActor
outlineActor SetMapper outlineMapper
eval [outlineActor GetProperty] SetColor $black

vtkPolyDataMapper outlineMapper2
outlineMapper2 SetInputConnection [outline GetOutputPort]

vtkActor outlineActor2
outlineActor2 SetMapper outlineMapper2
eval [outlineActor2 GetProperty] SetColor $black

# create planes
# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
ren1 SetViewport 0 0 0.5 1

vtkRenderer ren2
ren2 SetViewport 0.5 0 1 1

vtkRenderWindow renWin
renWin AddRenderer ren1
renWin AddRenderer ren2

vtkRenderWindowInteractor iren
iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 AddActor planeActor
ren2 AddActor outlineActor2
ren2 AddActor planeActor2

ren1 SetBackground 1 1 1
ren2 SetBackground 1 1 1

renWin SetSize 500 250

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange 3.95297 50
$cam1 SetFocalPoint 8.88908 0.595038 29.3342
$cam1 SetPosition -12.3332 31.7479 41.2387
$cam1 SetViewUp 0.060772 -0.319905 0.945498
ren2 SetActiveCamera [ren1 GetActiveCamera]

# render the image
#
renWin Render

iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize

# prevent the tk window from showing up then start the event loop
wm withdraw .
