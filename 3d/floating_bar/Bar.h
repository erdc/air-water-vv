#ifndef BOX_H
#define BOX_H
//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//
//   Demo code about
//
//     - collisions and contacts  with DEM
//
//   CHRONO
//   ------
//   Multibody dinamics engine
//
///////////////////////////////////////////////////

#include "chrono/physics/ChSystemDEM.h"
#include "chrono/physics/ChContactContainerDEM.h"
#include "chrono/solver/ChSolverDEM.h"

#include "chrono_irrlicht/ChIrrApp.h"

#include <irrlicht.h>

//using namespace std;
// Use the namespace of Chrono
using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;

void AddWall(std::shared_ptr<ChBody> body, const ChVector<>& dim, const ChVector<>& loc) {
    body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);
    auto box = std::make_shared<ChBoxShape>();
    box->GetBoxGeometry().Size = dim;
    box->GetBoxGeometry().Pos = loc;
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    body->AddAsset(box);
}

class cppChRigidBar
{
 public:
 cppChRigidBar(double glass_thickness,
	       double* gravityIn,
	       double* bar_center,
               double* bar_dim,
               double* L,
               double mass,
               double* inertia,
               double* free_x,
               double* free_r):
    gravity(gravityIn),
    time_step(0.00001),
    out_step(2000 * time_step),
    barId(100),
    mass(mass),
    pos(bar_center[0],bar_center[1],bar_center[2]),
    rot(1, 0, 0, 0),
    init_vel(0, 0, 0),
    binId(200),
    thickness(glass_thickness),
    free_x(free_x[0], free_x[1], free_x[2]),
    free_r(free_r[0], free_r[1], free_r[2])
  {
    // The following two lines are optional, since they are the default options. They are added for future reference,
    // i.e. when needed to change those models.
    msystem.SetContactForceModel(ChSystemDEM::ContactForceModel::Hertz);
    msystem.SetAdhesionForceModel(ChSystemDEM::AdhesionForceModel::Constant);
    
    msystem.Set_G_acc(ChVector<>(gravity[0],
                                 gravity[1],
                                 gravity[2]));
    
    // Create a material (will be used by both objects)
    material = std::make_shared<ChMaterialSurfaceDEM>();
    material->SetRestitution(0.1f);
    material->SetFriction(0.4f);
    material->SetAdhesion(0);  // Magnitude of the adhesion in Constant adhesion model
    
    // Create the falling bar
    bar = std::make_shared<ChBody>(ChMaterialSurfaceBase::DEM);
    
    bar->SetIdentifier(barId);
    bar->SetMass(mass);
    bar->SetPos(pos);
    bar->SetRot(rot);
    bar->SetPos_dt(init_vel);
    bar->SetBodyFixed(false);
    bar->SetMaterialSurface(material);
    
    bar->SetCollide(true);
    
    bar->GetCollisionModel()->ClearModel();
    bar->GetCollisionModel()->AddBox(0.5*bar_dim[0],
                                     0.5*bar_dim[1],
                                     0.5*bar_dim[2]);
    bar->GetCollisionModel()->BuildModel();
    bar->SetInertiaXX(ChVector<>(inertia[0],
                                 inertia[4],
                                 inertia[8]));
    
    barShape = std::make_shared<ChBoxShape>();
    barShape->GetBoxGeometry().Size = ChVector<>(0.5*bar_dim[0],
						 0.5*bar_dim[1],
						 0.5*bar_dim[2]);
    bar->AddAsset(barShape);
    
    mtexture = std::make_shared<ChTexture>();
    mtexture->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    bar->AddAsset(mtexture);
    
    msystem.AddBody(bar);
  
    // Create container
    bin = std::make_shared<ChBody>(ChMaterialSurfaceBase::DEM);
    
    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);
    bin->SetMaterialSurface(material);
    
    bin->GetCollisionModel()->ClearModel();
    AddWall(bin, ChVector<>(0.5*L[0], 0.5*L[1], 0.5*thickness), ChVector<>(0.5*L[0], 0.5*L[1], -thickness));
    AddWall(bin, ChVector<>(0.5*thickness, 0.5*L[1], 0.5*L[2]), ChVector<>(- 0.5*thickness, 0.5*L[1], 0.5*L[2]));
    AddWall(bin, ChVector<>(0.5*thickness, 0.5*L[1], 0.5*L[2]), ChVector<>( L[0] + 0.5*thickness, 0.5*L[1], 0.5*L[2]));
    AddWall(bin, ChVector<>(0.5*L[0], 0.5*thickness, 0.5*L[1]), ChVector<>(0.5*L[0], -0.5*thickness, 0.5*L[2]));
    AddWall(bin, ChVector<>(0.5*L[0], 0.5*thickness, 0.5*L[2]), ChVector<>(0.5*L[0],  L[1] + 0.5*thickness, 0.5*L[2]));
    bin->GetCollisionModel()->BuildModel();
    
    msystem.AddBody(bin);
    
    /* // Complete asset construction */
    /* application.AssetBindAll(); */
    /* application.AssetUpdateAll(); */
    
    // The soft-real-time cycle
    double time = 0.0;
    double out_time = 0.0;
  }
  
  double* gravity;
  double time_step;
  double out_step;
  int barId;
  double radius;
  double mass;
  ChVector<> pos;
  ChQuaternion<> rot;
  ChVector<> init_vel;
  int binId;
  double thickness;
  ChSystemDEM msystem;
  std::shared_ptr<ChMaterialSurfaceDEM> material;
  std::shared_ptr<ChBody> bar;
  std::shared_ptr<ChBoxShape> barShape;
  std::shared_ptr<ChTexture> mtexture;
  std::shared_ptr<ChBody> bin;
  double time = 0.0;
  double out_time = 0.0;
  ChVector<double> pos_last;
  ChMatrix33<double> A_last, A;
  ChVector<double> free_x, free_r;
  //setInertia_xx and xy
  double hx(double* x, double t)
  {
    ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, A_last);
    ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, A);
    return xNew(0) - x[0];
  }
  double hy(double* x, double t)
  {
    ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, A_last);
    ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, A);
    return xNew(1) - x[1];
  }
  double hz(double* x, double t)
  {
    ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, A_last);
    ChVector<double> xNew = ChTransform<double>::TransformLocalToParent(local, pos, A);
    return xNew(2) - x[2];
  }
  void step(double* force, double* torque, double dt)
  {
    pos_last = bar->GetPos();
    A_last = bar->GetA();
    bar->Empty_forces_accumulators();
    bar->Accumulate_force(ChVector<double>(force[0]*free_x(0),
                                           force[1]*free_x(1),
                                           force[2]*free_x(2)),
                          pos_last,
                          false);
    bar->Accumulate_torque(ChVector<double>(torque[0]*free_r(0),
                                            torque[1]*free_r(1),
                                            torque[2]*free_r(2)),
                           false);
    msystem.DoStepDynamics(dt);
    pos = bar->GetPos();
    A = bar->GetA();
  }
};

cppChRigidBar* newChRigidBar(double glass_thickness,
			     double* gravityIn,
			     double* bar_center,
			     double* bar_dim,
			     double* L,
			     double mass,
			     double* inertia,
			     double* free_x,
					     double* free_r)
{
  return new cppChRigidBar(glass_thickness,
			   gravityIn,
			   bar_center,
			   bar_dim,
			   L,
			   mass,
			   inertia,
			   free_x,
			   free_r);
}

#endif
