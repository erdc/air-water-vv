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
  std::vector< ChSharedPtr< ChBody > > * list =	system.Get_bodylist();
  ChSharedPtr<ChBody> bod = list->front();
  ChVector< double > acc = bod->GetPos_dtdt();
  ChVector< double > torque = bod->Get_Xtorque();
  ofstream myfile;
  myfile.open ("acceleration.txt", std::ios_base::app);
  myfile << system.GetChTime() << ",";
  myfile << torque(0) << "," << torque(1) << "," << torque(2) << ",";
  myfile << acc(0) << "," << acc(1) << "," << acc(2) << "\n";
  myfile.close();
}
  
class cppRigidBody {
 public:
  ChVector<> free_x;
  ChVector<> free_r;
  ChVector<> pos;
  ChVector<> pos_last;
  ChQuaternion<double> rot;
  ChQuaternion<double> rot_last;
  double mass;
  double* inertia;
  ChMatrix33<double> a_last;
  ChMatrix33<double> a;
  ChSharedPtr<ChBody> body;
  cppSystem* system;
  cppRigidBody(cppSystem* system,
               double* pos,
               double* rot,
               double mass,
               double* inertia,
               double* free_x,
               double* free_r);
  double hx(double* x, double t);
  double hy(double* x, double t);
  double hz(double* x, double t);
  void prestep(double* force, double* torque, double dt);
  void poststep();
};

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

  body = ChSharedPtr<ChBody>(new ChBody());
  // add body to system
  system->system.AddBody(body);
  // basic attributes of body
  pos = ChVector<>(posin[0], posin[1], posin[2]);
  rot = ChQuaternion<>(rotin[0], rotin[1], rotin[2], rotin[3]);
  body->SetPos(pos);
  body->SetRot(rot);
  body->SetInertiaXX(ChVector<>(1.,
                                1.,
                                inertia[2]));  // careful division by zero!
  a = body->GetA();
  a_last = body->GetA();
  pos = body->GetPos(); 
  pos_last = body->GetPos(); 
  body->SetMass(mass);
}

double cppRigidBody::hx(double* x, double t)
{
  a = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, a_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, a);
  return xNew(0) - x[0];
}
double cppRigidBody::hy(double* x, double t)
{
  a = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, a_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, a);
  return xNew(1) - x[1];
}
double cppRigidBody::hz(double* x, double t)
{
  a = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, a_last);
  ChVector<double> xNew = ChTransform<double>::TransformLocalToParent(local, pos, a);
  return xNew(2) - x[2];
}

void cppRigidBody::prestep(double* force, double* torque, double dt)
{
  pos_last = body->GetPos();
  a_last = body->GetA();
  pos_last = body->GetPos();
  rot_last = body->GetRot();
  body->Empty_forces_accumulators();
  body->Accumulate_force(ChVector<double>(force[0]*free_x(0),
                                          force[1]*free_x(1),
                                          force[2]*free_x(2)),
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
  a = body->GetA();
}


cppSystem * newSystem(double* gravity)
{
  return new cppSystem(gravity);
}



cppRigidBody * newRigidBody(cppSystem* system,
                            double* position,
                            double* rot,
                            double mass,
                            double* inertia,
                            double* free_x,
                            double* free_r)
{
  return new cppRigidBody(system,
                          position,
                          rot,
                          mass,
                          inertia,
                          free_x,
                          free_r);

}
