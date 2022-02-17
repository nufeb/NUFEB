/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkShadowMapPass.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// These internal classes are shared by vtkShadowMapBakerPass and
// vtkShadowMapPass.

#ifndef vtkShadowMapPassInternal_h
#define vtkShadowMapPassInternal_h

#include "vtkCamera.h"
#include "vtkTextureObject.h"

#include <vector>
class vtkShadowMapBakerPassTextures
{
public:
  std::vector<vtkSmartPointer<vtkTextureObject> > Vector;
};

class vtkShadowMapBakerPassLightCameras
{
public:
  std::vector<vtkSmartPointer<vtkCamera> > Vector;
};

#endif
// VTK-HeaderTest-Exclude: vtkShadowMapPassInternal.h
