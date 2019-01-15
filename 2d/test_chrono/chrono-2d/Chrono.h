// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Milad Rakhsha
// =============================================================================
//
// This file contains the chrono model to set up a FSI simulation.
// =============================================================================
#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChParticlesClones.h"
#include "chrono/core/ChFileutils.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"

#include <typeinfo>

#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveBinary.h"
#include "chrono/serialization/ChArchiveAsciiDump.h"
#include "chrono/serialization/ChArchiveJSON.h"

using namespace chrono;
using namespace chrono::collision;
std::string out_dir = "outputs/";
double time_step = 1e-3;

class cppMBDModel
{
  private:
    double container_center[3]; //(x,y,z)
    double container_dims[3];   //(height,width,thickness)
    double wallThickness = 0.1;

    double gravity[3];

    ChSystemNSC my_system;

    double time_step;
    int outframe = 0;
    double time = 0.0;

    double friction_const = 0.01;

    std::vector<std::shared_ptr<ChBodyEasySphere> > pt_body;

    ChState x;
    ChStateDelta v;
    double T;
    ChStateDelta a;
    ChVectorDynamic<> L;

  public:
    int num_particles;
    double *diam_;
    double *pos_;
    double *vel_;
    double *angular_vel_;

  public:
    ~cppMBDModel(){
        free(diam_);
        free(pos_);
        free(vel_);
        free(angular_vel_);
    }
    cppMBDModel(double m_timeStep,
                double *m_container_center,
                double *m_container_dims,
                double *m_gravity,
                int nParticle,
                double *ball_center,
                double *ball_radius,
                double *ball_density) : time_step(m_timeStep)
    {
        memcpy(container_center, m_container_center, 3 * sizeof(double));
        memcpy(container_dims, m_container_dims, 3 * sizeof(double));
        memcpy(gravity, m_gravity, 3 * sizeof(double));


        //printf("container_dims -first= %f,%f,%f\n", container_dims[0], container_dims[1], container_dims[2]);

        my_system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
        //std::cout << "Set_G_acc " << gravity[0] << " " << gravity[1] << " " << gravity[2] << " " << std::endl;

        const std::string rmCmd = (std::string("rm ") + out_dir + std::string("/*"));
        system(rmCmd.c_str());

        if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0)
        {
            double a = 0;
            std::cout << "Error creating directory " << out_dir << std::endl;
            std::cin >> a;
        }

        num_particles = nParticle;
        diam_ = (double *)malloc(num_particles * 1 * sizeof(double));
        pos_ = (double *)malloc(num_particles * 3 * sizeof(double));
        vel_ = (double *)malloc(num_particles * 3 * sizeof(double));
        angular_vel_ = (double *)malloc(num_particles * 3 * sizeof(double));

        //// Create ground body used for link all balls
        auto Ground = std::make_shared<ChBody>();
        Ground->SetBodyFixed(true);
        Ground->SetPos(ChVector<>(0.0, 0.0, -10.0));//so there is no effect to balls
        Ground->SetRot(ChQuaternion<>(1.0, 0.0, 0.0, 0.0));
        my_system.Add(Ground);
        
        //// add particles
        std::shared_ptr<ChLinkLockPointPlane> constraint_xy;

        for (int i=0;i<num_particles;++i)
        {
            diam_[i] = 2*ball_radius[i];
            auto msphereBody = std::make_shared<ChBodyEasySphere>(ball_radius[i],        // radius size
                                                                  ball_density[i], // density
                                                                  true,                // collide enable?
                                                                  false);              // visualization?
            ChVector<> mpos = ChVector<>(ball_center[3*i+0],ball_center[3*i+1],ball_center[3*i+2]);
            msphereBody->SetPos(mpos);
            ChVector<> mvel = ChVector<>(0.0,0.0,0.0);
            ChVector<> mangular_vel = ChVector<>(0.0,0.0,0.0);
            msphereBody->SetPos_dt(mvel);
            msphereBody->SetWvel_loc(mangular_vel);
            //printf("initialized particle %d (%f,%f,%f)\n", i, mpos.x(), mpos.y(), mpos.z());
            msphereBody->GetMaterialSurfaceNSC()->SetFriction(friction_const);

            double mass = ball_radius[i]*ball_radius[i]*M_PI*ball_density[i];
            msphereBody->SetMass(mass);
            ChVector<> inertia_of_disk = ChVector<>(0.25*mass*ball_radius[i]*ball_radius[i],//for 2D object
                                                    0.25*mass*ball_radius[i]*ball_radius[i],//for 2D object
                                                    0.5* mass*ball_radius[i]*ball_radius[i]);
            msphereBody->SetInertiaXX(inertia_of_disk);
            msphereBody->SetInertiaXY(ChVector<>(0.0,0.0,0.0));

            my_system.Add(msphereBody);
            pt_body.push_back(msphereBody);

            pos_[i * 3 + 0] = mpos.x();
            pos_[i * 3 + 1] = mpos.y();
            pos_[i * 3 + 2] = mpos.z();
            vel_[i * 3 + 0] = mvel.x();
            vel_[i * 3 + 1] = mvel.y();
            vel_[i * 3 + 2] = mvel.z();
            angular_vel_[i * 3 + 0] = mangular_vel.x();
            angular_vel_[i * 3 + 1] = mangular_vel.y();
            angular_vel_[i * 3 + 2] = mangular_vel.z();

            constraint_xy = std::make_shared<ChLinkLockPointPlane>();
            constraint_xy->Initialize(msphereBody, Ground, ChCoordsys<>(msphereBody->GetPos()));
            my_system.Add(constraint_xy);
        }

        // auto mat = std::make_shared<ChMaterialSurfaceNSC>();
        // mat->SetFriction(friction_const);
        // The inconsistency here is to resolve the chrono default container creation. Chrono container is at z=dimZ/2 not 0

        //container_center[2]-container_dims[2]/2 is because of utils::CreateBoxContainer
        // ChVector<> center(container_center[0], container_center[1], container_center[2] - container_dims[2] / 2);
        // ChVector<> boxDim(container_dims[0] / 2, container_dims[1] / 2, container_dims[2] / 2);
        //printf("container_dims = %f,%f,%f\n", container_dims[0], container_dims[1], container_dims[2]);
        //printf("container_center = %f,%f,%f\n", container_center[0], container_center[1], container_center[2]);
        //printf("boxDim = %f,%f,%f\n", boxDim.x(), boxDim.y(), boxDim.z());
        //printf("center = %f,%f,%f\n", center.x(), center.y(), center.z());

        // utils::CreateBoxContainer(&my_system, 0,
        //                           mat,//material
        //                           boxDim, wallThickness,
        //                           center, Q_from_AngAxis(0, VECT_Y),
        //                           true, false,//collide, y_up
        //                           true, true);//overlap, closed
        add_2d_walls_without_front_back(my_system,container_center,container_dims,friction_const,ball_radius[0]);
        writeThisFrame();


        my_system.Setup();

    }

    void add_2d_walls_without_front_back(ChSystemNSC& mphysicalSystem, double *box_center, double *box_dim, double friction_const, double ball_radius)
    {

        auto left_wall = std::make_shared<ChBodyEasyBox>(1.0, box_dim[1], box_dim[2], //width is 1
                                                        1000,//density. Fix anyway
                                                        true,//collision
                                                        false);//visual
        left_wall->SetPos(ChVector<>(box_center[0]-box_dim[0]*0.5-0.5, box_center[1], box_center[2]));
        left_wall->SetBodyFixed(true);
        mphysicalSystem.Add(left_wall);
        left_wall->GetMaterialSurfaceNSC()->SetFriction(friction_const);

        auto right_wall = std::make_shared<ChBodyEasyBox>(1.0, box_dim[1], box_dim[2], //width is 1
                                                        1000,//density. Fix anyway
                                                        true,//collision
                                                        false);//visual
        right_wall->SetPos(ChVector<>(box_center[0]+box_dim[0]*0.5+0.5, box_center[1], box_center[2]));
        right_wall->SetBodyFixed(true);
        mphysicalSystem.Add(right_wall);
        right_wall->GetMaterialSurfaceNSC()->SetFriction(friction_const);

        
        auto bottom_wall = std::make_shared<ChBodyEasyBox>(box_dim[0], 1.0, box_dim[2], //width is 1
                                                        1000,//density. Fix anyway
                                                        true,//collision
                                                        false);//visual
        bottom_wall->SetPos(ChVector<>(box_center[0], box_center[1]-box_dim[1]*0.5-0.5, box_center[2]));
        bottom_wall->SetBodyFixed(true);
        mphysicalSystem.Add(bottom_wall);
        bottom_wall->GetMaterialSurfaceNSC()->SetFriction(friction_const);

        
        auto top_wall = std::make_shared<ChBodyEasyBox>(box_dim[0], 1.0, box_dim[2], //width is 1
                                                        1000,//density. Fix anyway
                                                        true,//collision
                                                        false);//visual
        top_wall->SetPos(ChVector<>(box_center[0], box_center[1]+box_dim[1]*0.5+0.5, box_center[2]));
        top_wall->SetBodyFixed(true);
        mphysicalSystem.Add(top_wall);
        top_wall->GetMaterialSurfaceNSC()->SetFriction(friction_const);

        // auto front_wall = std::make_shared<ChBodyEasyBox>(box_dim[0], box_dim[1], 1, //width is 1
        //                                                 1000,//density. Fix anyway
        //                                                 true,//collision
        //                                                 false);//visual
        // front_wall->SetPos(ChVector<>(box_center[0], box_center[1], box_center[2]+ball_radius+0.5+1e-6));//use ball radius
        // front_wall->SetBodyFixed(true);
        // mphysicalSystem.Add(front_wall);
        // front_wall->GetMaterialSurfaceNSC()->SetFriction(0.0);

        
        // auto back_wall = std::make_shared<ChBodyEasyBox>(box_dim[0], box_dim[1], 1, //width is 1
        //                                                 1000,//density. Fix anyway
        //                                                 true,//collision
        //                                                 false);//visual
        // back_wall->SetPos(ChVector<>(box_center[0], box_center[1], box_center[2]-ball_radius-0.5-1e-6));
        // back_wall->SetBodyFixed(true);
        // mphysicalSystem.Add(back_wall);
        // back_wall->GetMaterialSurfaceNSC()->SetFriction(0.0);

    }
    void reset_timer()
    {
        my_system.ResetTimers();
    }
    void save_state_to_file()
    {
        std::string jsonfile = out_dir + "/foo_archive.json";
        ChStreamOutAsciiFile mfileo(jsonfile.c_str());
        ChArchiveOutJSON marchiveout(mfileo);
        my_system.ArchiveOUT(marchiveout);
    }
    void get_state_from_file()
    {
        std::string jsonfile = out_dir + "/foo_archive.json";
        ChStreamInAsciiFile mfilei(jsonfile.c_str());
        ChArchiveInJSON marchivein(mfilei);
        my_system.ArchiveIN(marchivein);
    }
    void save_state()
    {
        // save_state_to_file();
        // Create state
        x.Reset(my_system.GetNcoords_x(), &my_system);//call Setup first
        v.Reset(my_system.GetNcoords_v(), &my_system);
        a.Reset(my_system.GetNcoords_a(), &my_system);
        L.Reset(my_system.GetNconstr());
        my_system.StateSetup(x,v,a);
        std::cout<<"size of x "<<x.GetLength()<<","<<x.GetRows()<<","<<x.GetColumns()<<std::endl;
        std::cout<<"size of v "<<v.GetLength()<<","<<v.GetRows()<<","<<v.GetColumns()<<std::endl;
        std::cout<<"size of z "<<a.GetLength()<<","<<a.GetRows()<<","<<a.GetColumns()<<std::endl;

        my_system.StateGather(x,v,T);
        my_system.StateGatherAcceleration(a);
        my_system.StateGatherReactions(L);
    }
    void reset_state()
    {
        // get_state_from_file();
        my_system.StateScatter(x,v,T);
        my_system.StateScatterAcceleration(a);
        my_system.StateScatterReactions(L);
    }

    void set_state(double *_pos, double *_vel, double *_ang_vel)
    {
        for (int i = 0; i < num_particles; i++)
        {
            auto body = pt_body[i];
            body->SetPos(ChVector<double>(_pos[3 * i + 0], _pos[3 * i + 1], _pos[3 * i + 2]));
            body->SetPos_dt(ChVector<double>(_vel[3 * i + 0], _vel[3 * i + 1], _vel[3 * i + 2]));
            body->SetWvel_par(ChVector<double>(_ang_vel[3 * i + 0], _ang_vel[3 * i + 1], _ang_vel[3 * i + 2]));

            body->Update(true);
            body->SyncCollisionModels();
        }
    }

    void get_state(double *_pos, double *_vel, double *_ang_vel)
    {
        for (int i = 0; i < num_particles; i++)
        {
            auto body = pt_body[i];
            auto position = body->GetPos();
            auto velocity = body->GetPos_dt();
            auto angular_velocity = body->GetWvel_par();
            _pos[3 * i + 0] = position.x();
            _pos[3 * i + 1] = position.y();
            _pos[3 * i + 2] = position.z();
            _vel[3 * i + 0] = velocity.x();
            _vel[3 * i + 1] = velocity.y();
            _vel[3 * i + 2] = velocity.z();
            _ang_vel[3 * i + 0] = angular_velocity.x();
            _ang_vel[3 * i + 1] = angular_velocity.y();
            _ang_vel[3 * i + 2] = angular_velocity.z();
        }
    }
    double *step(double *forces, double *torques, double dt)
    {

        for (int i = 0; i < num_particles; i++)
        {
            ChVector<double> myForce = ChVector<double>(forces[3 * i + 0], forces[3 * i + 1], forces[3 * i + 2]);
            ChVector<double> myTorque = ChVector<double>(torques[3 * i + 0], torques[3 * i + 1], torques[3 * i + 2]);
            // No big difference
            // ChVector<double> myForce = ChVector<double>(forces[3 * i + 0], forces[3 * i + 1], 0.0);//2D simulation
            // ChVector<double> myTorque = ChVector<double>(0.0, 0.0, torques[3 * i + 2]);//2D simulation

            // auto body = my_system.Get_bodylist()->at(i+1);//since 1st body is ground
            auto body = pt_body[i];
            body->Empty_forces_accumulators();//empty force and torque in
            ChVector<> pos = body->GetPos();
            ChVector<> W = body->GetMass() * my_system.Get_G_acc();//since 1st body is ground
            body->Accumulate_force(myForce, pos, 0);//input force vector in the global coord
            body->Accumulate_torque(myTorque, 0);//input force vector in the global coord
            //printf("Add force=%f,%f,%f, mg=%f,%f,%f to particle %d\n", myForce.x(), myForce.y(), myForce.z(),
            //       W.x(), W.y(), W.z(), i);
        }

        int sync = std::round(dt / time_step);
        if (sync >= 1)
        {
            //printf("%d * DoStepChronoSystem with dt= %f\n", sync, dt / sync);
            for (int t = 0; t < sync; t++)
            {
                my_system.DoStepDynamics(dt / sync);////////////////////////solving system
                time += dt / sync;
            }
        }
        else
        {
            //printf("DoStepChronoSystem with dt= %f\n ", dt);
            my_system.DoStepDynamics(dt);
        }
        //printf("===================================================\n");
        //printf("chrono time : %f\n", my_system.GetChTime());
        for (int i = 0; i < num_particles; i++)
        {
            ChVector<> pos = pt_body[i]->GetPos();
            pos_[i * 3 + 0] = pos.x();
            pos_[i * 3 + 1] = pos.y();
            pos_[i * 3 + 2] = pos.z();
            ChVector<> vel = pt_body[i]->GetPos_dt();
            vel_[i * 3 + 0] = vel.x();
            vel_[i * 3 + 1] = vel.y();
            vel_[i * 3 + 2] = vel.z();
            //return the coordinate of angular velocity vector in the absolute system
            ChVector<> angular_vel = pt_body[i]->  GetWvel_par();//return//since 1st body is ground
//            ChVector<> angular_vel = my_system.Get_bodylist()->at(i+1)->  GetWvel_loc();

            angular_vel_[i * 3 + 0] = angular_vel.x();
            angular_vel_[i * 3 + 1] = angular_vel.y();
            angular_vel_[i * 3 + 2] = angular_vel.z();
            //printf("Chrono.h file: pos of particle %d = %f,%f,%f\n", i, pos.x(), pos.y(), pos.z());
            //printf("Chrono.h file: vel of particle %d = %f,%f,%f\n", i, vel.x(), vel.y(), vel.z());
            //printf("Chrono.h file: angular vel of particle %d = %f,%f,%f\n", i, angular_vel.x(), angular_vel.y(), angular_vel.z());
        }
    }
//
//    void calc_d_N_IBM(double *pos, double *output)
//    {
//        // double output[4];
//        double r = 1e6;
//        int p = 0;
//
//        for (int i = 0; i < num_particles; i++)
//        {
//            double x = pos_[i * 3 + 0] - pos[0];
//            double y = pos_[i * 3 + 1] - pos[1];
//            double z = pos_[i * 3 + 2] - pos[2];
//            double r0 = diam_[i] / 2;
//            if (x * x + y * y + z * z < (r + r0) * (r + r0))
//            {
//                p = i;
//                output[0] = std::pow(x * x + y * y + z * z, 0.5) - r0;
//                output[1] = (pos[0] - pos_[i * 3 + 0]) / (r + r0);
//                output[2] = (pos[1] - pos_[i * 3 + 1]) / (r + r0);
//                output[3] = (pos[2] - pos_[i * 3 + 2]) / (r + r0);
//                r = std::pow(x * x + y * y + z * z, 0.5) - r0;
//            }
//        }
//        // printf("point= %f,%f,%f, is close to: %f,%f,%f, normal =%f,%f,%f\n", pos[0], pos[1], pos[2],
//        //        particle_pos_diameter[p * 3 + 0], particle_pos_diameter[p * 3 + 1], particle_pos_diameter[p * 3 + 2],
//        //        output[0], output[1], output[2]);
//    }


    void writeThisFrame()
    {
        const std::string nameFluid =
            out_dir + std::string("particles") + std::to_string(outframe) + std::string(".csv");
        const std::string bodies =
            out_dir + std::string("bodies") + std::to_string(outframe) + std::string(".csv");

        std::ofstream fileParticles;
        fileParticles.open(nameFluid);
        std::stringstream ssFluidParticles;
        ssFluidParticles << "x,y,z,vx,vy,vz,wx,wy,wz,d\n";
        for (int j = 0; j < num_particles; j++)
        {
            ChVector<> mypos = pt_body[j]->GetPos();//since 1st body is ground
            ChVector<> myvel = pt_body[j]->GetPos_dt();
            ChVector<> myangular_vel = pt_body[j]->GetWvel_par();

            ssFluidParticles << mypos.x() << ", " << mypos.y() << ", " << mypos.z() << ", "
                             << myvel.x() << ", " << myvel.y() << ", " << myvel.z() << ", "
                             << myangular_vel.x() << ", " << myangular_vel.y() << ", " << myangular_vel.z() << ", "
                             << diam_[j]
                             << std::endl;
        }
        fileParticles << ssFluidParticles.str();
        fileParticles.close();
        utils::WriteBodies(&my_system,
                           bodies,
                           false,
                           true,
                           ",");
        outframe++;
    }
};

cppMBDModel *newMBDModel(double m_timeStep,
                         double *m_container_center,
                         double *m_container_dims,
                         double *m_gravity,
                         int nParticle,
                         double* ball_center,
                         double* ball_radius,
                         double* ball_density)
{
    return new cppMBDModel(m_timeStep, m_container_center, m_container_dims, m_gravity,
                           nParticle, ball_center, ball_radius, ball_density);
};
