/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: freebodyjoint.cpp                                       *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/

#ifndef FREEBODYJOINT_H
#define FREEBODYJOINT_H

#include "joint.h"


class FreeBodyJoint : public Joint{
public: 
  FreeBodyJoint();
  ~FreeBodyJoint();
  
  JointType GetType();
  bool ReadInJointData(std::istream& in);
  void WriteOutJointData(std::ostream& out);
  void ComputeLocalTransform();
  Matrix GetForward_sP();
  Matrix GetBackward_sP();
  void UpdateForward_sP( Matrix& sP);
  void UpdateBackward_sP( Matrix& sP);
  void ForwardKinematics();
  void BackwardKinematics();
};

#endif
