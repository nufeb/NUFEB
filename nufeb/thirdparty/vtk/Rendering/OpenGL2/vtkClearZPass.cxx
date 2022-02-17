/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkClearZPass.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkClearZPass.h"
#include "vtkObjectFactory.h"
#include <cassert>
#include "vtkRenderState.h"
#include "vtkRenderer.h"
#include "vtk_glew.h"

vtkStandardNewMacro(vtkClearZPass);

// ----------------------------------------------------------------------------
vtkClearZPass::vtkClearZPass()
{
  this->Depth=1.0;
}

// ----------------------------------------------------------------------------
vtkClearZPass::~vtkClearZPass()
{
}

// ----------------------------------------------------------------------------
void vtkClearZPass::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Depth:" << this->Depth << endl;
}

// ----------------------------------------------------------------------------
// Description:
// Perform rendering according to a render state \p s.
// \pre s_exists: s!=0
void vtkClearZPass::Render(const vtkRenderState *s)
{
  assert("pre: s_exists" && s!=0);
  (void)s;
  this->NumberOfRenderedProps=0;

  glDepthMask(GL_TRUE);
#if GL_ES_VERSION_3_0 != 1
  glClearDepth(this->Depth);
#else
  glClearDepthf(static_cast<GLclampf>(this->Depth));
#endif
  glClear(GL_DEPTH_BUFFER_BIT);
}
