#include "chrono/physics/ChSystemDEM.h"
#include <memory>
#include <iostream>
#include <fstream>
using namespace chrono;
using namespace std;




class cppSystem {
 public:
  ChSystem system;
  double* gravity;
  cppSystem(double* gravity);
  void step(double dt);
};

class cppRigidBody {
 public:
  ChVector<> free_x;
  ChVector<> free_r;
  ChVector<> pos;
  ChVector<> pos_last;
  ChVector<> vel;
  ChVector<> vel_last;
  ChVector<> acc;
  ChVector<> acc_last;
  ChVector<> angvel;
  ChVector<> angvel_last;
  ChVector<> angacc;
  ChVector<> angacc_last;
  ChMatrix33<double> rotm;
  ChMatrix33<double> rotm_last;
  ChQuaternion<double> rotq;
  ChQuaternion<double> rotq_last;
  ChVector<> F;
  ChVector<> F_last;
  ChVector<> M;
  ChVector<> M_last;
  double mass;
  /* ChVector <> inertia; */
  double* inertia;
  std::shared_ptr<ChBody> body;
  cppSystem* system;
  cppRigidBody(cppSystem* system,
               double* pos,
               double* rotq,
               double mass,
               double* inertia,
               double* free_x,
               double* free_r);
  double hx(double* x, double t);
  double hy(double* x, double t);
  double hz(double* x, double t);
  void prestep(double* force, double* torque, double dt);
  void poststep();
  void setRotation(double* quat);
  void setPosition(double* quat);
  void setConstraints(double* free_x, double* free_y);
  void setInertiaXX(double* inertia);
  void addSpring(double stiffness,
                 double damping,
                 double* fairlead,
                 double* anchor,
                 double rest_length);
};

cppSystem::cppSystem(double* gravity):
gravity(gravity)
{
  system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
}

void cppSystem::step(double dt)
{
  double dt2 = dt/20;
  for (int i = 0; i < 20; ++i) {
    system.DoStepDynamics(dt2);
  }
  std::vector< std::shared_ptr< ChBody > > * list =	system.Get_bodylist();
  std::shared_ptr<ChBody> bod = list->front();
  ChVector< double > acc = bod->GetPos_dtdt();
  ChVector< double > torque = bod->Get_Xtorque();
}

cppRigidBody::cppRigidBody(cppSystem* system,
                           double* posin,
                           double* rotin,
                           double mass,
                           double* inertia,
                           double* free_xin,
                           double* free_rin):
system(system),
  mass(mass),
  inertia(inertia),
  free_x(free_xin[0], free_xin[1], free_xin[2]),
  free_r(free_rin[0], free_rin[1], free_rin[2])
{

  body = std::make_shared<ChBody>();
  // add body to system
  system->system.AddBody(body);
  // basic attributes of body
  pos = ChVector<>(posin[0], posin[1], posin[2]);
  rotq = ChQuaternion<>(rotin[0], rotin[1], rotin[2], rotin[3]);
  body->SetPos(pos);
  body->SetRot(rotq);
  body->SetInertiaXX(ChVector<>(1.,
                                1.,
                                inertia[2]));  // careful division by zero!
  rotm = body->GetA();
  rotm_last = body->GetA();
  pos = body->GetPos();
  pos_last = body->GetPos();
  body->SetMass(mass);
}

double cppRigidBody::hx(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
  return xNew(0) - x[0];
}

double cppRigidBody::hy(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
  return xNew(1) - x[1];
}

double cppRigidBody::hz(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
  return xNew(2) - x[2];
}

void cppRigidBody::prestep(double* force, double* torque, double dt)
{
  /* step to call before running chrono system step */
  pos_last = body->GetPos();
  vel_last = body->GetPos_dt();
  acc_last = body->GetPos_dtdt();
  rotm_last = body->GetA();
  rotq_last = body->GetRot();
  angacc_last = body->GetWvel_loc();
  angvel_last = body->GetWvel_loc();
  F_last = body->Get_Xforce();
  M_last = body->Get_Xtorque();
  // apply external forces
  body->Empty_forces_accumulators();
  // calculate opposite force of gravity if free_x is 0
  double forceG[3]={0.,0.,0.};
  if (free_x(0) == 0) {forceG[0] = -system->system.Get_G_acc()(0)*body->GetMass();}
  if (free_x(1) == 0) {forceG[1] = -system->system.Get_G_acc()(1)*body->GetMass();}
  if (free_x(2) == 0) {forceG[2] = -system->system.Get_G_acc()(2)*body->GetMass();}
  body->Accumulate_force(ChVector<double>(forceG[0]+force[0]*free_x(0),
                                          forceG[1]+force[1]*free_x(1),
                                          forceG[2]+force[2]*free_x(2)),
                         pos_last,
                         false);
  body->Accumulate_torque(ChVector<double>(torque[0]*free_r(0),
                                           torque[1]*free_r(1),
                                           torque[2]*free_r(2)),
                          true);
}

void cppRigidBody::poststep()
{
  pos = body->GetPos();
  vel = body->GetPos_dt();
  acc = body->GetPos_dtdt();
  rotm = body->GetA();
  rotq = body->GetRot();
  angacc = body->GetWvel_loc();
  angvel = body->GetWvel_loc();
  F = body->Get_Xforce();
  M = body->Get_Xtorque();
}

void cppRigidBody::setPosition(double* position){
  body->SetPos(ChVector<>(position[0], position[1], position[2]));
}

void cppRigidBody::setRotation(double* quat) {
  body->SetRot(ChQuaternion<double>(quat[0], quat[1], quat[2], quat[3]));
}

void cppRigidBody::setConstraints(double* free_x_in, double* free_r_in){
  free_x = ChVector<>(free_x_in[0], free_x_in[1], free_x_in[2]);
  free_r = ChVector<>(free_r_in[0], free_r_in[1], free_r_in[2]);
}

void cppRigidBody::setInertiaXX(double* inertia){
  body->SetInertiaXX(ChVector<>(inertia[0], inertia[1], inertia[2]));
}


void cppRigidBody::addSpring(double stiffness,
                             double damping,
                             double* fairlead,
                             double* anchor,
                             double rest_length)
{
  std::shared_ptr<ChLinkSpring> spring = std::make_shared<ChLinkSpring>();
  std::shared_ptr<ChBody> anchor_body = std::make_shared<ChBody>();
  anchor_body->SetPos(ChVector<>(anchor[0], anchor[1], anchor[2]));
  anchor_body->SetBodyFixed(true);
  system->system.AddBody(anchor_body);
  spring->Initialize(body,
                     anchor_body,
                     true, // true for pos relative to bodies
                     ChVector<>(fairlead[0], fairlead[1], fairlead[2]),
                     ChVector<>(0.,0.,0.),
                     false,  // true for auto rest length (distance between body1 and body2)
                     rest_length);
  spring->Set_SpringK(stiffness);
  spring->Set_SpringR(damping);
  system->system.AddLink(spring);
}

cppSystem * newSystem(double* gravity)
{
  return new cppSystem(gravity);
}



cppRigidBody * newRigidBody(cppSystem* system,
                            double* position,
                            double* rotq,
                            double mass,
                            double* inertia,
                            double* free_x,
                            double* free_r)
{
  return new cppRigidBody(system,
                          position,
                          rotq,
                          mass,
                          inertia,
                          free_x,
                          free_r);
}
