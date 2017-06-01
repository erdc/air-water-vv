#ifndef BAR_H
#define BAR_H

#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemDEM.h"
#include "chrono/physics/ChContactContainerDEM.h"
#include "chrono/solver/ChSolverDEM.h"
#include "chrono_fea/ChElementCableANCF.h"
#include "chrono_fea/ChBuilderBeam.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_irrlicht/ChIrrApp.h"

#include <irrlicht.h>

using namespace chrono;
using namespace irr;
using namespace fea;
using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;

void AddWall(std::shared_ptr<ChBody> body, const ChVector<>& dim, const ChVector<>& loc)
{
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
    bar = std::make_shared<ChBody>();//ChMaterialSurfaceBase::DEM);
    
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
  
    // Create the falling bar
    bar2 = std::make_shared<ChBody>();//ChMaterialSurfaceBase::DEM);
    
    bar2->SetIdentifier(barId+1);
    bar2->SetMass(mass/2.0);
    bar2->SetPos(pos - ChVector<>(bar_dim[0]/2.0 + bar_dim[0]/4.0,0.0,0.0));
    bar2->SetRot(rot);
    bar2->SetPos_dt(init_vel);
    bar2->SetBodyFixed(false);
    bar2->SetMaterialSurface(material);
    
    bar2->SetCollide(true);
    
    bar2->GetCollisionModel()->ClearModel();
    bar2->GetCollisionModel()->AddBox(0.25*bar_dim[0],
    				      0.5*bar_dim[1],
    				      0.5*bar_dim[2]);
    bar2->GetCollisionModel()->BuildModel();
    bar2->SetInertiaXX(ChVector<>(inertia[0]/2.,
				  inertia[4]/2.,
				  inertia[8]/2.));
    
    bar2Shape = std::make_shared<ChBoxShape>();
    bar2Shape->GetBoxGeometry().Size = ChVector<>(0.25*bar_dim[0],
						  0.5*bar_dim[1],
						  0.5*bar_dim[2]);
    bar2->AddAsset(bar2Shape);
    
    mtexture = std::make_shared<ChTexture>();
    mtexture->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    bar2->AddAsset(mtexture);
    
    msystem.AddBody(bar2);


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
    //AddWall(bin, ChVector<>(0.5*L[0], 0.5*L[1], 0.5*thickness), ChVector<>(0.5*L[0], 0.5*L[1], -thickness));
    AddWall(bin, ChVector<>(0.5*thickness, 0.5*L[1], 0.5*L[2]), ChVector<>(- 0.5*thickness, 0.5*L[1], 0.5*L[2]));
    AddWall(bin, ChVector<>(0.5*thickness, 0.5*L[1], 0.5*L[2]), ChVector<>( L[0] + 0.5*thickness, 0.5*L[1], 0.5*L[2]));
    //AddWall(bin, ChVector<>(0.5*L[0], 0.5*thickness, 0.5*L[1]), ChVector<>(0.5*L[0], -0.5*thickness, 0.5*L[2]));
    //AddWall(bin, ChVector<>(0.5*L[0], 0.5*thickness, 0.5*L[2]), ChVector<>(0.5*L[0],  L[1] + 0.5*thickness, 0.5*L[2]));
    bin->GetCollisionModel()->BuildModel();
    
    msystem.AddBody(bin);
    
    /* // Complete asset construction */
    /* application.AssetBindAll(); */
    /* application.AssetUpdateAll(); */
    
    // The soft-real-time cycle
    double time = 0.0;
    double out_time = 0.0;
    // 2. Create the mesh that will contain the finite elements, and add it to the system

    auto mesh = std::make_shared<ChMesh>();
    auto mesh2 = std::make_shared<ChMesh>();

    msystem.Add(mesh);
    msystem.Add(mesh2);


    // 3. Create a material for the beam finite elements.

    //    Note that each FEA element type requires some corresponding
    //    type of material. Here we will use ChElementCableANCF elements:
    //    they use a material of type ChBeamSectionCable, so let's do

    auto beam_material = std::make_shared<ChBeamSectionCable>();
    beam_material->SetDiameter(0.01);
    beam_material->SetYoungModulus(0.01e9);
    beam_material->SetBeamRaleyghDamping(0.01);


    // 4. Create the nodes

    //    - We use a simple for() loop to create nodes along the cable.
    //    - Nodes for ChElementCableANCF must be of ChNodeFEAxyzD class;
    //      i.e. each node has 6 coordinates: {position, direction}, where
    //      direction is the tangent to the cable.
    //    - Each node must be added to the mesh, ex.  mesh.Add(my_node)
    //    - To make things easier in the following, we store node pointers
    //      into an optional 'beam_nodes' array, i.e. a std::vector<>, later we
    //      can use such array for easy creation of elements between the nodes.

    std::vector<std::shared_ptr<ChNodeFEAxyzD> > beam_nodes;

    ChVector<> bar_attach(bar_center[0]+0.5*bar_dim[0],
			  bar_center[1]+0.5*bar_dim[1],
			  bar_center[2]+0.5*bar_dim[2]);
    ChVector<> wall_attach(L[0],
			   L[1],
			   L[2]);
    auto d = wall_attach - bar_attach;
    double length = d.Length();
    int N_nodes = 16;
    double dl = length/(N_nodes-1.0);
    for (int in = 0; in < N_nodes; ++in) {
        // i-th node position
      ChVector<> position(bar_attach+(in*dl/length)*d);
      
      // i-th node direction
      ChVector<> direction=d/length;
      
      // create the node
      auto node = std::make_shared<ChNodeFEAxyzD>(position, direction);
      
      // add it to mesh
      mesh->AddNode(node);

      // add it to the auxiliary beam_nodes
      beam_nodes.push_back(node);
    }


    // 5. Create the elements

    //    - We use a simple for() loop to create elements between the
    //      nodes that we already created.
    //    - Each element must be set with the ChBeamSectionCable material
    //      that we already created
    //    - Each element must be added to the mesh, ex.  mesh.Add(my_element)

    for (int ie = 0; ie < N_nodes - 1; ++ie) {
        // create the element
        auto element = std::make_shared<ChElementCableANCF>();

        // set the connected nodes (pick two consecutive nodes in our beam_nodes container)
        element->SetNodes(beam_nodes[ie], beam_nodes[ie + 1]);

        // set the material
        element->SetSection(beam_material);

        // add it to mesh
        mesh->AddElement(element);
    }


    // 6. Add constraints

    //    - Constraints can be applied to FEA nodes
    //    - For the ChNodeFEAxyzD there are specific constraints that
    //      can be used to connect them to a ChBody, namely
    //      ChLinkPointFrame and ChLinkDirFrame
    //    - To attach one end of the beam to the ground, we need a
    //      'truss' ChBody that is fixed.
    //    - Note. An alternative, only when the node must be fixed 
    //      to absolute reference, is not using constraints, and just
    //      use: beam_nodes[0]->SetFixed(true);  (but would fix also dir)


    // 7. Make the finite elements visible in the 3D view

    //   - FEA fisualization can be managed via an easy
    //     ChVisualizationFEAmesh helper class.
    //     (Alternatively you could bypass this and output .dat
    //     files at each step, ex. for VTK or Matalb postprocessing)
    //   - This will automatically update a triangle mesh (a ChTriangleMeshShape
    //     asset that is internally managed) by setting proper
    //     coordinates and vertex colours as in the FEA elements.
    //   - Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    //     postprocessor that can handle a coloured ChTriangleMeshShape).
    //   - Do not forget AddAsset() at the end!

    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(mesh.get()));
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ANCF_BEAM_AX);
    mvisualizebeamA->SetColorscaleMinMax(-0.005, 0.005);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    mesh->AddAsset(mvisualizebeamA);

    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*(mesh.get()));
    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);  // E_GLYPH_NODE_CSYS for ChNodeFEAxyzrot
    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.005);
    mvisualizebeamC->SetZbufferHide(false);
    mesh->AddAsset(mvisualizebeamC);

    // 4. Create the nodes

    //    - We use a simple for() loop to create nodes along the cable.
    //    - Nodes for ChElementCableANCF must be of ChNodeFEAxyzD class;
    //      i.e. each node has 6 coordinates: {position, direction}, where
    //      direction is the tangent to the cable.
    //    - Each node must be added to the mesh, ex.  mesh.Add(my_node)
    //    - To make things easier in the following, we store node pointers
    //      into an optional 'beam_nodes' array, i.e. a std::vector<>, later we
    //      can use such array for easy creation of elements between the nodes.

    std::vector<std::shared_ptr<ChNodeFEAxyzD> > beam_nodes2;

    ChVector<> bar_attach2(bar_center[0]-0.5*bar_dim[0],
			  bar_center[1]-0.5*bar_dim[1],
			  bar_center[2]+0.5*bar_dim[2]);
    ChVector<> wall_attach2(0,0,L[2]);
    auto d2 = wall_attach2 - bar_attach2;
    double length2 = d2.Length();
    double dl2 = length2/(N_nodes-1.0);
    for (int in = 0; in < N_nodes; ++in) {
        // i-th node position
      ChVector<> position(bar_attach2+(in*dl2/length2)*d2);
      
      // i-th node direction
      ChVector<> direction=d2/length2;
      
      // create the node
      auto node = std::make_shared<ChNodeFEAxyzD>(position, direction);
      
      // add it to mesh
      mesh2->AddNode(node);

      // add it to the auxiliary beam_nodes
      beam_nodes2.push_back(node);
    }


    // 5. Create the elements

    //    - We use a simple for() loop to create elements between the
    //      nodes that we already created.
    //    - Each element must be set with the ChBeamSectionCable material
    //      that we already created
    //    - Each element must be added to the mesh, ex.  mesh.Add(my_element)

    for (int ie = 0; ie < N_nodes - 1; ++ie) {
        // create the element
        auto element = std::make_shared<ChElementCableANCF>();

        // set the connected nodes (pick two consecutive nodes in our beam_nodes container)
        element->SetNodes(beam_nodes2[ie], beam_nodes2[ie + 1]);

        // set the material
        element->SetSection(beam_material);

        // add it to mesh
        mesh2->AddElement(element);
    }


    // 6. Add constraints

    //    - Constraints can be applied to FEA nodes
    //    - For the ChNodeFEAxyzD there are specific constraints that
    //      can be used to connect them to a ChBody, namely
    //      ChLinkPointFrame and ChLinkDirFrame
    //    - To attach one end of the beam to the ground, we need a
    //      'truss' ChBody that is fixed.
    //    - Note. An alternative, only when the node must be fixed 
    //      to absolute reference, is not using constraints, and just
    //      use: beam_nodes[0]->SetFixed(true);  (but would fix also dir)

    //auto truss = std::make_shared<ChBody>();
    //truss->SetBodyFixed(true);
    //system.Add(truss);

    //pin bars
    
    // Define two quaternions representing:
    // - a rotation of -90 degrees around x (z2y)
    ChQuaternion<> z2y;
    z2y.Q_from_AngAxis( CH_C_PI / 2, ChVector<>(1, 0, 0));
    // Revolute joint between ground and crank.
    // The rotational axis of a revolute joint is along the Z axis of the
    // specified joint coordinate frame.  Here, we apply the 'z2y' rotation to
    // align it with the Y axis of the global reference frame.
    auto revolute_ground_crank = std::make_shared<ChLinkRevolute>();
    revolute_ground_crank->SetName("revolute_ground_crank");
    revolute_ground_crank->Initialize(bar2,
				      bar,
				      true,
				      ChFrame<>(ChVector<>( 0.25*bar_dim[0],
							   -0.45*bar_dim[1],
							   -0.5*bar_dim[2]),
						z2y),
				      ChFrame<>(ChVector<>(-0.5*bar_dim[0],
							   -0.45*bar_dim[1],
							   -0.5*bar_dim[2]),z2y));
    msystem.AddLink(revolute_ground_crank);

    auto revolute_ground_crank2 = std::make_shared<ChLinkRevolute>();
    revolute_ground_crank2->SetName("revolute_ground_crank2");
    revolute_ground_crank2->Initialize(bar2,
				       bar,
				       true,
				       ChFrame<>(ChVector<>( 0.25*bar_dim[0],
							     0.45*bar_dim[1],
							    -0.5*bar_dim[2]),z2y),
				       ChFrame<>(ChVector<>(-0.5*bar_dim[0],
							     0.45*bar_dim[1],
							    -0.5*bar_dim[2]),z2y));
    msystem.AddLink(revolute_ground_crank2);

    // lock an end of the wire to the truss
    auto constraint_pos = std::make_shared<ChLinkPointFrame>();
    constraint_pos->Initialize(beam_nodes[N_nodes-1], bin);
    msystem.Add(constraint_pos);

    auto constraint_cyl = std::make_shared<ChLinkPointFrame>();
    constraint_cyl->Initialize(beam_nodes[0], bar);
    msystem.Add(constraint_cyl);

    // lock an end of the wire to the truss
    auto constraint_pos2 = std::make_shared<ChLinkPointFrame>();
    constraint_pos2->Initialize(beam_nodes2[N_nodes-1], bin);
    msystem.Add(constraint_pos2);

    auto constraint_cyl2 = std::make_shared<ChLinkPointFrame>();
    constraint_cyl2->Initialize(beam_nodes2[0], bar);
    msystem.Add(constraint_cyl2);

    // 7. Make the finite elements visible in the 3D view

    //   - FEA fisualization can be managed via an easy
    //     ChVisualizationFEAmesh helper class.
    //     (Alternatively you could bypass this and output .dat
    //     files at each step, ex. for VTK or Matalb postprocessing)
    //   - This will automatically update a triangle mesh (a ChTriangleMeshShape
    //     asset that is internally managed) by setting proper
    //     coordinates and vertex colours as in the FEA elements.
    //   - Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    //     postprocessor that can handle a coloured ChTriangleMeshShape).
    //   - Do not forget AddAsset() at the end!

    auto mvisualizebeamA2 = std::make_shared<ChVisualizationFEAmesh>(*(mesh2.get()));
    mvisualizebeamA2->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ANCF_BEAM_AX);
    mvisualizebeamA2->SetColorscaleMinMax(-0.005, 0.005);
    mvisualizebeamA2->SetSmoothFaces(true);
    mvisualizebeamA2->SetWireframe(false);
    mesh2->AddAsset(mvisualizebeamA2);

    auto mvisualizebeamC2 = std::make_shared<ChVisualizationFEAmesh>(*(mesh2.get()));
    mvisualizebeamC2->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);  // E_GLYPH_NODE_CSYS for ChNodeFEAxyzrot
    mvisualizebeamC2->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC2->SetSymbolsThickness(0.006);
    mvisualizebeamC2->SetSymbolsScale(0.005);
    mvisualizebeamC2->SetZbufferHide(false);
    mesh2->AddAsset(mvisualizebeamC2);

    // 8. Configure the solver and timestepper

    //    - the default SOLVER_SOR of Chrono is not able to manage stiffness matrices
    //      as required by FEA! we must switch to a different solver.
    //    - We pick the SOLVER_MINRES solver and we configure it.
    //    - Note that if you build the MKL module, you could use the more precise MKL solver.

    // Change solver
    msystem.SetSolverType(ChSystem::SOLVER_MINRES);
    msystem.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    msystem.SetMaxItersSolverSpeed(200);
    msystem.SetMaxItersSolverStab(200);
    msystem.SetTolForce(1e-10);

    // Change integrator:
    msystem.SetIntegrationType(chrono::ChSystem::INT_EULER_IMPLICIT_LINEARIZED);  // default: fast, 1st order
    //msystem.SetIntegrationType(chrono::ChSystem::INT_HHT);  // precise, slower, might iterate each step
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
  std::shared_ptr<ChBody> bar, bar2;
  std::shared_ptr<ChBoxShape> barShape, bar2Shape;
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
