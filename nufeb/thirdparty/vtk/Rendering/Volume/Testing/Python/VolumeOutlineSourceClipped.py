#!/usr/bin/env python
import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

# Simple volume rendering example.
reader = vtk.vtkSLCReader()
reader.SetFileName(VTK_DATA_ROOT + "/Data/sphere.slc")

# Create transfer functions for opacity and color
opacityTransferFunction = vtk.vtkPiecewiseFunction()
opacityTransferFunction.AddPoint(0, 0.0)
opacityTransferFunction.AddPoint(30, 0.0)
opacityTransferFunction.AddPoint(80, 0.5)
opacityTransferFunction.AddPoint(255, 0.5)

colorTransferFunction = vtk.vtkColorTransferFunction()
colorTransferFunction.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
colorTransferFunction.AddRGBPoint(64.0, 1.0, 0.0, 0.0)
colorTransferFunction.AddRGBPoint(128.0, 0.0, 0.0, 1.0)
colorTransferFunction.AddRGBPoint(192.0, 0.0, 1.0, 0.0)
colorTransferFunction.AddRGBPoint(255.0, 0.0, 0.2, 0.0)

# Create properties, mappers, volume actors, and ray cast function
volumeProperty = vtk.vtkVolumeProperty()
volumeProperty.SetColor(colorTransferFunction)
volumeProperty.SetScalarOpacity(opacityTransferFunction)
volumeProperty.SetInterpolationTypeToLinear()
volumeProperty.ShadeOn()

ren1 = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren1)

renWin.SetSize(300, 300)

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
renWin.SetMultiSamples(0)
ren1.SetBackground(0.1, 0.2, 0.4)

# Translate the volume to center it at (0,0,0)
userTrans = vtk.vtkTransform()
userTrans.PostMultiply()
userTrans.Identity()
userTrans.Translate(-24.5, -24.5, -24.5)

# Clipping planes are in world coords
plane1 = vtk.vtkPlane()
plane1.SetOrigin(-24, 0, 0)
plane1.SetNormal(1, 0, 0)

plane2 = vtk.vtkPlane()
plane2.SetOrigin(24, 0, 0)
plane2.SetNormal(-1, 0, 0)

plane3 = vtk.vtkPlane()
plane3.SetOrigin(0, -15, 0)
plane3.SetNormal(0.163176, 0.925417, -0.342020)

plane4 = vtk.vtkPlane()
plane4.SetOrigin(0, 24, 0)
plane4.SetNormal(0, -1, 0)

plane5 = vtk.vtkPlane()
plane5.SetOrigin(0, 0, -24)
plane5.SetNormal(0, 0, 1)

plane6 = vtk.vtkPlane()
plane6.SetOrigin(0, 0, 24)
plane6.SetNormal(0, 0, -1)

clippingPlanes = vtk.vtkPlaneCollection()
clippingPlanes.AddItem(plane1)
clippingPlanes.AddItem(plane2)
clippingPlanes.AddItem(plane3)
clippingPlanes.AddItem(plane4)
clippingPlanes.AddItem(plane5)
clippingPlanes.AddItem(plane6)

# Cropping planes are in data coords
volumeMapper1 = vtk.vtkFixedPointVolumeRayCastMapper()
volumeMapper1.SetInputConnection(reader.GetOutputPort())
volumeMapper1.CroppingOn()
volumeMapper1.SetCroppingRegionPlanes(16, 33, 16, 33, 16, 33)
volumeMapper1.SetClippingPlanes(clippingPlanes)

volume1 = vtk.vtkVolume()
volume1.SetMapper(volumeMapper1)
volume1.SetProperty(volumeProperty)

outline1 = vtk.vtkVolumeOutlineSource()
outline1.SetVolumeMapper(volumeMapper1)
outline1.GenerateFacesOn()
outline1.GenerateScalarsOn()

preTrans1 = vtk.vtkTransformPolyDataFilter()
preTrans1.SetInputConnection(outline1.GetOutputPort())
preTrans1.SetTransform(userTrans)

outlineClip1 = vtk.vtkClipClosedSurface()
outlineClip1.SetInputConnection(preTrans1.GetOutputPort())
outlineClip1.SetClippingPlanes(clippingPlanes)
outlineClip1.GenerateFacesOff()
outlineClip1.GenerateOutlineOn()
outlineClip1.SetScalarModeToColors()
outlineClip1.SetClipColor(1, 1, 0)
outlineClip1.SetActivePlaneId(2)
outlineClip1.SetActivePlaneColor(0, 1, 0)

postTrans1 = vtk.vtkTransformPolyDataFilter()
postTrans1.SetInputConnection(outlineClip1.GetOutputPort())
postTrans1.SetTransform(userTrans.GetInverse())

outlineMapper1 = vtk.vtkDataSetMapper()
outlineMapper1.SetInputConnection(postTrans1.GetOutputPort())

outlineActor1 = vtk.vtkActor()
outlineActor1.SetMapper(outlineMapper1)

volume1.SetUserTransform(userTrans)

outlineActor1.SetUserTransform(userTrans)

ren1.AddViewProp(outlineActor1)
ren1.AddViewProp(volume1)

volumeMapper1.SetCroppingRegionFlagsToFence()

outline1.GenerateScalarsOn()

ren1.ResetCamera()
ren1.GetActiveCamera().Zoom(1.35)

renWin.Render()

iren.Initialize()
#iren.Start()
