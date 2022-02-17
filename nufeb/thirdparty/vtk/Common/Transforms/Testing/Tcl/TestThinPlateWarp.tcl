package require vtk

# warp an image with a thin plate spline

# first, create an image to warp
vtkImageGridSource imageGrid
imageGrid SetGridSpacing 16 16 0
imageGrid SetGridOrigin 0 0 0
imageGrid SetDataExtent 0 255 0 255 0 0
imageGrid SetDataScalarTypeToUnsignedChar

vtkLookupTable table
table SetTableRange 0 1
table SetValueRange 1.0 0.0
table SetSaturationRange 0.0 0.0
table SetHueRange 0.0 0.0
table SetAlphaRange 0.0 1.0
table Build

vtkImageMapToColors alpha
alpha SetInputConnection [imageGrid GetOutputPort]
alpha SetLookupTable table

vtkBMPReader reader1
reader1 SetFileName "$VTK_DATA_ROOT/Data/masonry.bmp"

vtkImageBlend blend
blend AddInputConnection 0 [reader1 GetOutputPort]
blend AddInputConnection 0 [alpha GetOutputPort]

# next, create a ThinPlateSpline transform

vtkPoints p1
p1 SetNumberOfPoints 8
p1 SetPoint 0 0 0 0
p1 SetPoint 1 0 255 0
p1 SetPoint 2 255 0 0
p1 SetPoint 3 255 255 0
p1 SetPoint 4 96 96 0
p1 SetPoint 5 96 159 0
p1 SetPoint 6 159 159 0
p1 SetPoint 7 159 96 0

vtkPoints p2
p2 SetNumberOfPoints 8
p2 SetPoint 0 0 0 0
p2 SetPoint 1 0 255 0
p2 SetPoint 2 255 0 0
p2 SetPoint 3 255 255 0
p2 SetPoint 4 96 159 0
p2 SetPoint 5 159 159 0
p2 SetPoint 6 159 96 0
p2 SetPoint 7 96 96 0

vtkThinPlateSplineTransform transform
  transform SetSourceLandmarks p2
  transform SetTargetLandmarks p1
  transform SetBasisToR2LogR
# you must invert the transform before passing it to vtkImageReslice
  transform Inverse

vtkImageReslice reslice
  reslice SetInputConnection [blend GetOutputPort]
  reslice SetResliceTransform transform
  reslice SetInterpolationModeToLinear

# set the window/level to 255.0/127.5 to view full range
vtkImageViewer viewer
viewer SetInputConnection [reslice GetOutputPort]
viewer SetColorWindow 255.0
viewer SetColorLevel 127.5
viewer SetZSlice 0
viewer Render




