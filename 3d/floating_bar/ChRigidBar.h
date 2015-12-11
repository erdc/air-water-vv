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
#include "chrono/lcp/ChLcpSolverDEM.h"

#include "chrono_irrlicht/ChIrrApp.h"

#include <irrlicht.h>

// Use the namespace of Chrono
using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;

void AddWall(ChSharedPtr<ChBody> body, const ChVector<>& dim, const ChVector<>& loc) {
    body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);

    ChSharedPtr<ChBoxShape> box(new ChBoxShape);
    box->GetBoxGeometry().Size = dim;
    box->GetBoxGeometry().Pos = loc;
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    body->AddAsset(box);
}

class cppChRigidBar
{
 public:
 cppChRigidBar(double* bar_center,
               double* bar_dim,
               double* L,
               double* gravityIn,
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
    width(2),
    length(2),
    height(1),
    thickness(0.1),
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
    material = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
    material->SetRestitution(0.1f);
    material->SetFriction(0.4f);
    material->SetAdhesion(0);  // Magnitude of the adhesion in Constant adhesion model
    
    // Create the falling bar
    bar = ChSharedPtr<ChBody>(new ChBody(ChMaterialSurfaceBase::DEM));
    
    bar->SetIdentifier(barId);
    bar->SetMass(mass);
    bar->SetPos(pos);
    bar->SetRot(rot);
    bar->SetPos_dt(init_vel);
    bar->SetBodyFixed(false);
    bar->SetMaterialSurface(material);
    
    bar->SetCollide(true);
    
    bar->GetCollisionModel()->ClearModel();
    bar->GetCollisionModel()->AddBox(bar_dim[0],
                                     bar_dim[1],
                                     bar_dim[2]);
    bar->GetCollisionModel()->BuildModel();
    bar->SetInertiaXX(ChVector<>(inertia[0],
                                 inertia[4],
                                 inertia[8]));
    
    barShape = ChSharedPtr<ChBoxShape>(new ChBoxShape);
    barShape->GetBoxGeometry().SetLengths(ChVector<>(bar_dim[0],
                                                     bar_dim[1],
                                                     bar_dim[2]));
    bar->AddAsset(barShape);
    
    mtexture = ChSharedPtr<ChTexture>(new ChTexture);
    mtexture->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    bar->AddAsset(mtexture);
    
    msystem.AddBody(bar);
  
    // Create container
    bin = ChSharedPtr<ChBody>(new ChBody(ChMaterialSurfaceBase::DEM));
    
    bin->SetIdentifier(binId);
    bin->SetMass(1);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(ChQuaternion<>(1, 0, 0, 0));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);
    bin->SetMaterialSurface(material);
    
    bin->GetCollisionModel()->ClearModel();
    AddWall(bin, ChVector<>(width, thickness, length), ChVector<>(0, 0, 0));
    // AddWall(bin, ChVector<>(thickness, height, length), ChVector<>(-width + thickness, height, 0));
    // AddWall(bin, ChVector<>(thickness, height, length), ChVector<>(width - thickness, height, 0));
    // AddWall(bin, ChVector<>(width, height, thickness), ChVector<>(0, height, -length + thickness));
    // AddWall(bin, ChVector<>(width, height, thickness), ChVector<>(0, height, length - thickness));
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
  double width;
  double length;
  double height;
  double thickness;
  ChSystemDEM msystem;
  ChSharedPtr<ChMaterialSurfaceDEM> material;
  ChSharedPtr<ChBody> bar;
  ChSharedPtr<ChBoxShape> barShape;
  ChSharedPtr<ChTexture> mtexture;
  ChSharedPtr<ChBody> bin;
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
    //ChBody.h
    //on the body
    //    setPosition()//
    //bar->SetPos(ChVector<>(dummy_pos));
    //bar->SetPos_dt(ChVector<>(dummy_vel));
    //  setLinearVel()
    //  setTorque
    //  setForce
    bar->Empty_forces_accumulators();
    ChVector<> Fstar;
    bar->Accumulate_force(ChVector<double>(force[0]*free_x(0),
                                           force[1]*free_x(1),
                                           force[2]*free_x(2)),
                          ChVector<>(0,0,0),true);
    bar->Accumulate_torque(ChVector<double>(torque[0]*free_r(0),
                                            torque[1]*free_r(1),
                                            torque[2]*free_r(2)),
                           true);
    msystem.DoStepDynamics(dt);
    pos = bar->GetPos();
    A = bar->GetA();
  }
};

cppChRigidBar* newChRigidBar(double* gravityIn,
                             double* bar_center,
                             double* bar_dim,
                             double* L,
                             double mass,
                             double* inertia,
                             double* free_x,
                             double* free_r)
{
  return new cppChRigidBar(gravityIn,
                           bar_center,
                           bar_dim,
                           L,
                           mass,
                           inertia,
                           free_x,
                           free_r);
}
