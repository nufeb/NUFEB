package require vtk
package require vtkinteraction

# create a rendering window and renderer
vtkRenderer ren1

vtkRenderWindow renWin
renWin AddRenderer ren1
renWin StereoCapableWindowOn

vtkRenderWindowInteractor iren
iren SetRenderWindow renWin

vtkGenericEnSightReader reader
# Make sure all algorithms use the composite data pipeline
vtkCompositeDataPipeline cdp
reader SetDefaultExecutivePrototype cdp
reader SetCaseFileName "$VTK_DATA_ROOT/Data/EnSight/blow1_bin.case"
reader SetTimeValue 1

vtkGeometryFilter geom
geom SetInputConnection [reader GetOutputPort]

vtkHierarchicalPolyDataMapper mapper
mapper SetInputConnection [geom GetOutputPort]
mapper SetColorModeToMapScalars
mapper SetScalarModeToUsePointFieldData
mapper ColorByArrayComponent "displacement" 0
mapper SetScalarRange 0 2.08

vtkActor actor
actor SetMapper mapper

# assign our actor to the renderer
ren1 AddActor actor

# enable user interface interactor
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize

[ren1 GetActiveCamera] SetPosition 99.3932 17.6571 -22.6071
[ren1 GetActiveCamera] SetFocalPoint 3.5 12 1.5
[ren1 GetActiveCamera] SetViewAngle 30
[ren1 GetActiveCamera] SetViewUp 0.239617 -0.01054 0.97081
ren1 ResetCameraClippingRange

renWin Render

# prevent the tk window from showing up then start the event loop
wm withdraw .

reader SetDefaultExecutivePrototype {}
