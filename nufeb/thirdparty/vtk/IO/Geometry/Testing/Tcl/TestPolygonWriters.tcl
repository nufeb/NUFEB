package require vtk
package require vtkinteraction

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# read data
#
vtkPolyDataReader input
  input SetFileName "$VTK_DATA_ROOT/Data/brainImageSmooth.vtk"

#
# generate vectors
vtkCleanPolyData clean
  clean SetInputConnection [input GetOutputPort]

vtkWindowedSincPolyDataFilter smooth
  smooth SetInputConnection [clean GetOutputPort]
  smooth GenerateErrorVectorsOn
  smooth GenerateErrorScalarsOn
  smooth Update

vtkPolyDataMapper mapper
  mapper SetInputConnection [smooth GetOutputPort]
  eval mapper SetScalarRange [[smooth GetOutput] GetScalarRange]

vtkActor brain
  brain SetMapper mapper


# Add the actors to the renderer, set the background and size
#
ren1 AddActor brain

renWin SetSize 320 240

[ren1 GetActiveCamera] SetPosition 149.653 -65.3464 96.0401
[ren1 GetActiveCamera] SetFocalPoint 146.003 22.3839 0.260541
[ren1 GetActiveCamera] SetViewAngle 30
[ren1 GetActiveCamera] SetViewUp -0.255578 -0.717754 -0.647695
[ren1 GetActiveCamera] SetClippingRange 79.2526 194.052

iren Initialize
renWin Render

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .

#
# If the current directory is writable, then test the witers
#

if {[catch {set channel [open "test.tmp" "w"]}] == 0 } {
   close $channel
   file delete -force "test.tmp"

   #
   #
   # test the writers
   vtkDataSetWriter dsw
   dsw SetInputConnection [smooth GetOutputPort]
   dsw SetFileName "brain.dsw"
   dsw Write
   file delete -force "brain.dsw"

   vtkPolyDataWriter pdw
   pdw SetInputConnection [smooth GetOutputPort]
   pdw SetFileName "brain.pdw"
   pdw Write
   file delete -force "brain.pdw"

   if { [info command "vtkIVWriter"] != "" } {
      vtkIVWriter iv
      iv SetInputConnection [smooth GetOutputPort]
      iv SetFileName "brain.iv"
      iv Write
      file delete -force "brain.iv"
   }

   #
   # the next writers only handle triangles
   vtkTriangleFilter triangles
   triangles SetInputConnection [smooth GetOutputPort]

   if { [info command "vtkIVWriter"] != "" } {
      vtkIVWriter iv2
      iv2 SetInputConnection [triangles GetOutputPort]
      iv2 SetFileName "brain2.iv"
      iv2 Write
      file delete -force "brain2.iv"
   }

   if { [info command "vtkIVWriter"] != "" } {
      vtkExtractEdges edges
      edges SetInputConnection [triangles GetOutputPort]
      vtkIVWriter iv3
      iv3 SetInputConnection [edges GetOutputPort]
      iv3 SetFileName "brain3.iv"
      iv3 Write
      file delete -force "brain3.iv"
   }

   vtkBYUWriter byu
   byu SetGeometryFileName "brain.g"
   byu SetScalarFileName "brain.s"
   byu SetDisplacementFileName "brain.d"
   byu SetInputConnection [triangles GetOutputPort]
   byu Write
   file delete -force "brain.g"
   file delete -force "brain.s"
   file delete -force "brain.d"

   vtkMCubesWriter mcubes
   mcubes SetInputConnection [triangles GetOutputPort]
   mcubes SetFileName "brain.tri"
   mcubes SetLimitsFileName "brain.lim"
   mcubes Write
   file delete -force "brain.lim"
   file delete -force "brain.tri"

   vtkSTLWriter stl
   stl SetInputConnection [triangles GetOutputPort]
   stl SetFileName "brain.stl"
   stl Write
   file delete -force "brain.stl"

   vtkSTLWriter stlBinary
   stlBinary SetInputConnection [triangles GetOutputPort]
   stlBinary SetFileName "brainBinary.stl"
   stlBinary SetFileType 2
   stlBinary Write
   file delete -force "brainBinary.stl"

#   vtkCGMWriter cgm
#   cgm SetInputConnection [triangles GetOutputPort]
#   cgm SetFileName "brain.cgm"
#   cgm Write
#   file delete -force "brain.cgm"
}



