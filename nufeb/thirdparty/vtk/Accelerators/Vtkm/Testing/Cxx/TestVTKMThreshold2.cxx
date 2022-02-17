/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestThreshold.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkmThreshold.h"
#include "vtkElevationFilter.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRTAnalyticSource.h"


int TestVTKMThreshold2(int argc, char *argv[])
{
  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  renWin->AddRenderer(ren.GetPointer());
  iren->SetRenderWindow(renWin.GetPointer());

  //---------------------------------------------------
  // Test using different thresholding methods
  //---------------------------------------------------
  vtkNew<vtkRTAnalyticSource> source;

  vtkNew<vtkElevationFilter> elevation;
  elevation->SetInputConnection(source->GetOutputPort());
  elevation->SetScalarRange(0.0, 1.0);
  elevation->SetLowPoint(-10.0, -10.0, -10.0);
  elevation->SetHighPoint(10.0, 10.0, 10.0);

  vtkNew<vtkmThreshold> threshold;
  threshold->SetInputConnection(elevation->GetOutputPort());
  threshold->SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RTData");

  double L=100;
  double U=200;
  threshold->ThresholdBetween(L,U);
  threshold->SetAllScalars(0);
  threshold->Update();

  threshold->UseContinuousCellRangeOn();
  threshold->Update();

  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetInputConnection(threshold->GetOutputPort());

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(surface->GetOutputPort());
  mapper->ScalarVisibilityOn();
  mapper->SetScalarModeToUsePointFieldData();
  mapper->SelectColorArray("Elevation");
  mapper->SetScalarRange(0.0, 1.0);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());

  ren->AddActor(actor.GetPointer());
  ren->ResetCamera();
  renWin->Render();

  int retVal = vtkRegressionTestImage(renWin.GetPointer());
  if(retVal == vtkRegressionTester::DO_INTERACTOR)
  {
  iren->Start();
  retVal = vtkRegressionTester::PASSED;
  }
  return (!retVal);
}
