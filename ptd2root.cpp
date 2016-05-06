// - Implementation of PTD2Root
// Example on access to data in 'brio' files from flsimulate:
// Uses the 'PTD' data bank 
// 1) Access all available data in PTD data banks
// 2) Write data to a flat TTree ROOT file

// Ourselves
#include "ptd2root.h"

// Standard Library
// Third Party
// - A
// This Project
// Macro which automatically implements the interface needed
// to enable the module to be loaded at runtime
// The first argument is the typename
// The second is the string key used to access the module in pipeline
// scripts. This must be globally unique.

DPP_MODULE_REGISTRATION_IMPLEMENT(PTD2Root,"PTD2Root");

// Construct
PTD2Root::PTD2Root() : dpp::base_module()
{ 
  filename_output_="default.root";
}

// Destruct
PTD2Root::~PTD2Root() {
  // MUST reset module at destruction
  this->reset();
}

// Initialize
void PTD2Root::initialize(const datatools::properties& myConfig,
			  datatools::service_manager& flServices,
			  dpp::module_handle_dict_type& /*moduleDict*/) {
  
  // Throw logic exception if we've already initialized this instance
  DT_THROW_IF(this->is_initialized(),
	      std::logic_error,
	      "PTD2Root already initialized");
  // Extract the filename_out key from the supplied config, if
  // the key exists. datatools::properties throws an exception if
  // the key isn't in the config, so catch this if thrown and don't do
  // anything
  try {
    myConfig.fetch("filename_out",this->filename_output_);
  } catch (std::logic_error& e) {
  }

  // Look for services
  if (flServices.has("geometry"));
  {
    const geomtools::geometry_service& GS = flServices.get<geomtools::geometry_service> ("geometry");

    // initialize geometry manager
    //    std::cout << "Initialize geo manager " << std::endl;
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }

  std::cout << "In INIT: create TFile " << std::endl;
  // Next all root file output here

  hfile_ = new TFile(filename_output_.c_str(),"RECREATE","Output file of Simulation data");
  hfile_->cd();

  tree_ = new TTree("PTD","PTD");
  tree_->SetDirectory(hfile_);
  
  // header data
  tree_->Branch("header.runnumber",&header_.runnumber_);
  tree_->Branch("header.eventnumber",&header_.eventnumber_);
  tree_->Branch("header.date",&header_.date_);
  tree_->Branch("header.runtype",&header_.runtype_);
  tree_->Branch("header.simulated",&header_.simulated_);

  // particle data
  tree_->Branch("particle.nofparticles",&particle_.nofparticles_);
  tree_->Branch("particle.nofgammaonly",&particle_.nofgammas_);
  tree_->Branch("particle.particleID",&particle_.particle_id_);
  tree_->Branch("particle.charge",&particle_.charge_);
  tree_->Branch("particle.vertex_type",&particle_.vertex_type_);
  tree_->Branch("particle.vertex_x",&particle_.vertex_x_);
  tree_->Branch("particle.vertex_y",&particle_.vertex_y_);
  tree_->Branch("particle.vertex_z",&particle_.vertex_z_);
  tree_->Branch("particle.foil_dir_x",&particle_.foil_dir_x_);
  tree_->Branch("particle.foil_dir_y",&particle_.foil_dir_y_);
  tree_->Branch("particle.foil_dir_z",&particle_.foil_dir_z_);
  tree_->Branch("particle.traj_length",&particle_.traj_length_);
  tree_->Branch("particle.traj_cl_del",&particle_.traj_cluster_delayed_);
  tree_->Branch("particle.traj_cl_del_time",&particle_.traj_cluster_delayed_time_);
  tree_->Branch("particle.calo_associated",&particle_.calo_associated_);
  tree_->Branch("particle.calo_type",&particle_.calo_type_);
  tree_->Branch("particle.calo_energy",&particle_.calo_energy_);
  tree_->Branch("particle.calo_sigma_energy",&particle_.calo_sigma_energy_);
  tree_->Branch("particle.calo_time",&particle_.calo_time_);
  tree_->Branch("particle.calo_sigma_time",&particle_.calo_sigma_time_);
  tree_->Branch("particle.calo_side",&particle_.calo_side_);
  tree_->Branch("particle.calo_column",&particle_.calo_column_);
  tree_->Branch("particle.calo_row",&particle_.calo_row_);
  tree_->Branch("particle.calo_wall",&particle_.calo_wall_);
  tree_->Branch("particle.calo_loc_x",&particle_.calo_loc_x_);
  tree_->Branch("particle.calo_loc_y",&particle_.calo_loc_y_);
  tree_->Branch("particle.calo_loc_z",&particle_.calo_loc_z_);

  this->_set_initialized(true);
}

// Process
dpp::base_module::process_status PTD2Root::process(datatools::things& workItem) {  
  // Local variables
  
  // particle event data
  std::vector<int> particleid;
  std::vector<int> charge;
  std::vector<int> vertextype;
  std::vector<double> vertex_x;
  std::vector<double> vertex_y;
  std::vector<double> vertex_z;
  std::vector<double> foil_dir_x;
  std::vector<double> foil_dir_y;
  std::vector<double> foil_dir_z;
  std::vector<double> traj_length;
  std::vector<int> traj_cl_delayed;
  std::vector<double> traj_cl_delayed_time;
  std::vector<int> caloassociated;
  std::vector<int> calotype;
  std::vector<double> caloenergy;
  std::vector<double> calosigmaenergy;
  std::vector<double> calotime;
  std::vector<double> calosigmatime;
  std::vector<int> caloside;
  std::vector<int> calocolumn;
  std::vector<int> calorow;
  std::vector<int> calowall;
  std::vector<double> calo_loc_x;
  std::vector<double> calo_loc_y;
  std::vector<double> calo_loc_z;
  
  // Access the workItem

  // look for reconstructed data
  if(workItem.has("PTD"))
    {
      const snemo::datamodel::particle_track_data & PTD = workItem.get<snemo::datamodel::particle_track_data>("PTD");

      // Extract particle track data
      if (PTD.has_particles()) {
	particle_.nofparticles_ = PTD.get_number_of_particles();

	for (size_t i=0; i<PTD.get_number_of_particles();++i) {
	  const snemo::datamodel::particle_track & the_particle = PTD.get_particle(i);

	  particleid.push_back(the_particle.get_track_id());
	  charge.push_back(the_particle.get_charge());
	  
	  // first the vertices
	  if (the_particle.has_vertices()) {
	    for (unsigned int i=0; i<the_particle.get_vertices().size();++i) {
	      const geomtools::blur_spot & vertex = the_particle.get_vertices().at(i).get();
	      const geomtools::vector_3d & translation  = vertex.get_placement().get_translation();	
	      vertex_x.push_back(translation.x());
	      vertex_y.push_back(translation.y());
	      vertex_z.push_back(translation.z());
	      if (vertex.has_geom_id()) { // vertex is on calorimeter (calo or xcalo or gveto)
		vertextype.push_back(vertex.get_geom_id().get_type());
	      }
	      else if (translation.x() < 0.01 && translation.x() > -0.01) { // test for zero x foil coordinate
		vertextype.push_back(0);
		// get the line direction at the foil vertex position
		if (the_particle.has_trajectory()) {
		  const snemo::datamodel::tracker_trajectory & the_trajectory = the_particle.get_trajectory();
		  const snemo::datamodel::base_trajectory_pattern & the_base_pattern = the_trajectory.get_pattern();
		  if (the_base_pattern.get_pattern_id()=="line") {
		    const geomtools::line_3d & the_shape = (const geomtools::line_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil_dir_x.push_back(direction.x());
		    foil_dir_y.push_back(direction.y());
		    foil_dir_z.push_back(direction.z());
		  }
		  else {
		    const geomtools::helix_3d & the_shape = (const geomtools::helix_3d&)the_base_pattern.get_shape();
		    geomtools::vector_3d direction = the_shape.get_direction_on_curve(the_shape.get_first());
		    foil_dir_x.push_back(direction.x());
		    foil_dir_y.push_back(direction.y());
		    foil_dir_z.push_back(direction.z());
		  }
		}
	      }
	      else // vertex is on wire (i.e. neither calo, nor foil)
		vertextype.push_back(1);

	    }
	  }

	  // then associated calorimeter hits
	  if (the_particle.has_associated_calorimeter_hits()) {
	    for (unsigned int i=0; i<the_particle.get_associated_calorimeter_hits().size();++i) {
	      const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = the_particle.get_associated_calorimeter_hits().at(i).get();
	      const geomtools::mapping & the_mapping = geometry_manager_->get_mapping();
	      if (! the_mapping.validate_id(calo_hit.get_geom_id())) {
		std::vector<geomtools::geom_id> gids;
		the_mapping.compute_matching_geom_id(calo_hit.get_geom_id(), gids); // front calo block = last entry
		const geomtools::geom_info & info = the_mapping.get_geom_info(gids.back()); // in vector gids
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();	
		calo_loc_x.push_back(loc.x());
		calo_loc_y.push_back(loc.y());
		calo_loc_z.push_back(loc.z());
	      }
	      else {
		const geomtools::geom_info & info = the_mapping.get_geom_info(calo_hit.get_geom_id());
		const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();	
		calo_loc_x.push_back(loc.x());
		calo_loc_y.push_back(loc.y());
		calo_loc_z.push_back(loc.z());
	      }
	      caloassociated.push_back(1);
	      caloenergy.push_back(calo_hit.get_energy());
	      calosigmaenergy.push_back(calo_hit.get_sigma_energy());
	      calotime.push_back(calo_hit.get_time());
	      calosigmatime.push_back(calo_hit.get_sigma_time());
	      calotype.push_back(calo_hit.get_geom_id().get_type());
	      
	      if (calo_hit.get_geom_id ().get_type () == 1302)
		{
		  // CALO
		  caloside.push_back(calo_hit.get_geom_id().get(1));
		  calocolumn.push_back(calo_hit.get_geom_id().get(2));
		  calorow.push_back(calo_hit.get_geom_id().get(3));
		}
	      if (calo_hit.get_geom_id ().get_type () == 1232)
		{
		  // XCALO
		  caloside.push_back(calo_hit.get_geom_id().get(1));
		  calocolumn.push_back(calo_hit.get_geom_id().get(3));
		  calowall.push_back(calo_hit.get_geom_id().get(2));
		  calorow.push_back(calo_hit.get_geom_id().get(0));
		}
	      if (calo_hit.get_geom_id ().get_type () == 1252)
		{
		  // GVETO
		  caloside.push_back(calo_hit.get_geom_id().get(1));
		  calocolumn.push_back(calo_hit.get_geom_id().get(3));
		  calowall.push_back(calo_hit.get_geom_id().get(2));
		}
	    }
	  }
	  else { // not a calo vertex, fill blanks
	    caloassociated.push_back(-1); // non particle
	    caloenergy.push_back(-1.0);
	    calosigmaenergy.push_back(-1.0);
	    calotime.push_back(-1.0);
	    calosigmatime.push_back(-1.0);
	    calotype.push_back(-1);
	    caloside.push_back(-1);
	    calocolumn.push_back(-1);
	    calowall.push_back(-1);
	    calorow.push_back(-1);
	    calo_loc_x.push_back(1.0e3);
	    calo_loc_y.push_back(1.0e4);
	    calo_loc_z.push_back(1.0e4);
	  }
	  
	  
	  // then the trajectory length, always
	  if (the_particle.has_trajectory()) {
	    const snemo::datamodel::tracker_trajectory & the_trajectory = the_particle.get_trajectory();
	    const snemo::datamodel::base_trajectory_pattern & the_base_pattern = the_trajectory.get_pattern();
	    const geomtools::i_shape_1d & the_shape = the_base_pattern.get_shape();
	    traj_length.push_back(the_shape.get_length());
	    
	    const snemo::datamodel::tracker_cluster & the_cluster = the_trajectory.get_cluster();
	    traj_cl_delayed.push_back((int)the_cluster.is_delayed());
	    if (the_cluster.is_delayed()>0)
	      traj_cl_delayed_time.push_back(the_cluster.get_hit(0).get_delayed_time());
	    else
	      traj_cl_delayed_time.push_back(-1.0);
	  }
	  else {
	    traj_length.push_back(-1.0);
	    traj_cl_delayed.push_back(-1);
	    traj_cl_delayed_time.push_back(-1.0);
	  }
	}
      }
      else 
	particle_.nofparticles_ = 0;

      // Extract un-associated calorimeter hits
      if (PTD.has_non_associated_calorimeters ()) {
	particle_.nofgammas_ = PTD.get_non_associated_calorimeters().size();

	for (unsigned int i=0; i<PTD.get_non_associated_calorimeters().size();++i) {
	  const snemo::datamodel::calibrated_calorimeter_hit & un_calo_hit = PTD.get_non_associated_calorimeters().at(i).get();

 	  const geomtools::mapping & the_mapping = geometry_manager_->get_mapping();
	  if (! the_mapping.validate_id(un_calo_hit.get_geom_id())) {
	    std::vector<geomtools::geom_id> gids;
	    the_mapping.compute_matching_geom_id(un_calo_hit.get_geom_id(), gids); // front calo block = last entry
	    const geomtools::geom_info & info = the_mapping.get_geom_info(gids.back()); // in vector gids
	    const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();	
	    calo_loc_x.push_back(loc.x());
	    calo_loc_y.push_back(loc.y());
	    calo_loc_z.push_back(loc.z());
	  }
	  else {
	    const geomtools::geom_info & info = the_mapping.get_geom_info(un_calo_hit.get_geom_id());
	    const geomtools::vector_3d & loc  = info.get_world_placement().get_translation();	
	    calo_loc_x.push_back(loc.x());
	    calo_loc_y.push_back(loc.y());
	    calo_loc_z.push_back(loc.z());
	  }
	  
	  caloassociated.push_back(0);
	  caloenergy.push_back(un_calo_hit.get_energy());
	  calosigmaenergy.push_back(un_calo_hit.get_sigma_energy());
	  calotime.push_back(un_calo_hit.get_time());
	  calosigmatime.push_back(un_calo_hit.get_sigma_time());
	  calotype.push_back(un_calo_hit.get_geom_id().get_type());

	  if (un_calo_hit.get_geom_id ().get_type () == 1302)
	    {
	      // CALO
	      caloside.push_back(un_calo_hit.get_geom_id().get(1));
	      calocolumn.push_back(un_calo_hit.get_geom_id().get(2));
	      calorow.push_back(un_calo_hit.get_geom_id().get(3));
	    }
	  if (un_calo_hit.get_geom_id ().get_type () == 1232)
	    {
	      // XCALO
	      caloside.push_back(un_calo_hit.get_geom_id().get(1));
	      calocolumn.push_back(un_calo_hit.get_geom_id().get(3));
	      calowall.push_back(un_calo_hit.get_geom_id().get(2));
	      calorow.push_back(un_calo_hit.get_geom_id().get(0));
	    }
	  if (un_calo_hit.get_geom_id ().get_type () == 1252)
	    {
	      // GVETO
	      caloside.push_back(un_calo_hit.get_geom_id().get(1));
	      calocolumn.push_back(un_calo_hit.get_geom_id().get(3));
	      calowall.push_back(un_calo_hit.get_geom_id().get(2));
	    }
	  // calo unassociated, append blanks
	  particleid.push_back(-1); // not a particle
	  charge.push_back(-100);
	  vertex_x.push_back(1.e3); // calo only
	  vertex_y.push_back(1.e4);
	  vertex_z.push_back(1.e4);
	  vertex_x.push_back(1.e3); // double fill for vertex
	  vertex_y.push_back(1.e4);
	  vertex_z.push_back(1.e4);
	  vertextype.push_back(-1);
	  vertextype.push_back(-1);
	  traj_length.push_back(-1.0);
	  traj_cl_delayed.push_back(-1);
	  traj_cl_delayed_time.push_back(-1.0);
	}
      }
      else {
	particle_.nofgammas_ = 0;
      }

      particle_.particle_id_ = &particleid;
      particle_.charge_ = &charge;
      particle_.vertex_type_ = &vertextype;
      particle_.vertex_x_ = &vertex_x;
      particle_.vertex_y_ = &vertex_y;
      particle_.vertex_z_ = &vertex_z;
      particle_.foil_dir_x_ = &foil_dir_x;
      particle_.foil_dir_y_ = &foil_dir_y;
      particle_.foil_dir_z_ = &foil_dir_z;
      particle_.traj_length_ = &traj_length;
      particle_.traj_cluster_delayed_ = &traj_cl_delayed;
      particle_.traj_cluster_delayed_time_ = &traj_cl_delayed_time;
      particle_.calo_associated_ = &caloassociated;
      particle_.calo_type_ = &calotype;
      particle_.calo_energy_ = &caloenergy;
      particle_.calo_sigma_energy_ = &calosigmaenergy;
      particle_.calo_time_ = &calotime;
      particle_.calo_sigma_time_ = &calosigmatime;
      particle_.calo_side_ = &caloside;
      particle_.calo_column_ = &calocolumn;
      particle_.calo_row_ = &calorow;
      particle_.calo_wall_ = &calowall;
      particle_.calo_loc_x_ = &calo_loc_x;
      particle_.calo_loc_y_ = &calo_loc_y;
      particle_.calo_loc_z_ = &calo_loc_z;

    }

  // look for event header
  if(workItem.has("EH"))
    {
      const snemo::datamodel::event_header & EH = workItem.get<snemo::datamodel::event_header>("EH");
      //      std::cout << "In process: found EH event header " << std::endl;
      header_.runnumber_ = EH.get_id ().get_run_number ();
      header_.eventnumber_ = EH.get_id ().get_event_number ();
      header_.date_ = 0;
      header_.runtype_ = 0;
      header_.simulated_ = (EH.is_simulated () ? true : false);
    }

  tree_->Fill();
  
  // MUST return a status, see ref dpp::processing_status_flags_type
  return dpp::base_module::PROCESS_OK;
}

// Reset
void PTD2Root::reset() {
  // write the output, finished streaming
  hfile_->cd();
  tree_->Write();
  hfile_->Close(); // 
  std::cout << "In reset: finished conversion, file closed " << std::endl;
  
  // clean up
  delete hfile_;
  filename_output_ = "default.root";
  this->_set_initialized(false);
}
