package require vtk
package require vtkinteraction
package require vtktesting

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin SetMultiSamples 0
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create pipeline
#
vtkSLCReader slc
  slc SetFileName "$VTK_DATA_ROOT/Data/nut.slc"
  slc Update

set types  "Char UnsignedChar Short UnsignedShort Int UnsignedInt Long UnsignedLong Float Double"
set i 3
foreach vtkType $types {
  vtkImageClip clip$vtkType
    clip$vtkType SetInputConnection [slc GetOutputPort]
    clip$vtkType SetOutputWholeExtent -1000 1000 -1000 1000 $i $i
  incr i 2
  vtkImageCast castTo$vtkType
    castTo$vtkType SetOutputScalarTypeTo$vtkType
    castTo$vtkType SetInputConnection [clip$vtkType GetOutputPort]
    castTo$vtkType ClampOverflowOn

  vtkContourFilter iso$vtkType
    iso$vtkType SetInputConnection [castTo$vtkType GetOutputPort]
    iso$vtkType GenerateValues 1 30 30

  vtkPolyDataMapper iso${vtkType}Mapper
    iso${vtkType}Mapper SetInputConnection [iso$vtkType GetOutputPort]
    iso${vtkType}Mapper SetColorModeToMapScalars

  vtkActor iso${vtkType}Actor
    iso${vtkType}Actor SetMapper iso${vtkType}Mapper
    ren1 AddActor iso${vtkType}Actor
}

vtkOutlineFilter outline
  outline SetInputConnection [slc GetOutputPort]
vtkPolyDataMapper outlineMapper
  outlineMapper SetInputConnection [outline GetOutputPort]
vtkActor outlineActor
  outlineActor SetMapper outlineMapper
  outlineActor VisibilityOff

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 ResetCamera
[ren1 GetActiveCamera] SetViewAngle 30
[ren1 GetActiveCamera] Elevation 20
[ren1 GetActiveCamera] Azimuth 20
[ren1 GetActiveCamera] Zoom 1.5
ren1 ResetCameraClippingRange

ren1 SetBackground 0.9 .9 .9
renWin SetSize 200 200
iren Initialize

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .


