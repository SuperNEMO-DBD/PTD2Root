//! \file PTD2Root.h
//! \brief User processing module for flreconstruct
//! \details Process a things object and convert data to ROOT file output.
#ifndef PTD2ROOT_HH
#define PTD2ROOT_HH

// Standard Library
#include <vector>
#include <iostream>
#include <string>

// Third Party
#include "TFile.h"
#include "TTree.h"
#include <boost/foreach.hpp>


// - Bayeux
#include "bayeux/dpp/base_module.h"
#include "bayeux/datatools/service_manager.h"
#include "bayeux/geomtools/manager.h"
#include "bayeux/geomtools/geometry_service.h"
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"

// - Falaise
#include "falaise/snemo/datamodels/particle_track_data.h"
#include "falaise/snemo/datamodels/particle_track.h"
#include "falaise/snemo/datamodels/calibrated_calorimeter_hit.h"
#include "falaise/snemo/datamodels/base_trajectory_pattern.h"
#include "falaise/snemo/datamodels/event_header.h"



// This Project
typedef struct HeaderEventStorage{
  int runnumber_;
  int eventnumber_;
  int date_;
  int runtype_;
  bool simulated_;
} headereventstorage;

typedef struct ParticleEventStorage{
  int nofparticles_;
  int nofgammas_;
  std::vector<int>* particle_id_;
  std::vector<int>* charge_;
  std::vector<int>* vertex_type_;
  std::vector<double>* vertex_x_;
  std::vector<double>* vertex_y_;
  std::vector<double>* vertex_z_;
  std::vector<double>* foil_dir_x_;
  std::vector<double>* foil_dir_y_;
  std::vector<double>* foil_dir_z_;
  std::vector<double>* traj_length_;
  std::vector<int>* traj_cluster_delayed_;
  std::vector<double>* traj_cluster_delayed_time_;
  std::vector<int>* calo_associated_;
  std::vector<int>* calo_type_;
  std::vector<double>* calo_energy_;
  std::vector<double>* calo_sigma_energy_;
  std::vector<double>* calo_time_;
  std::vector<double>* calo_sigma_time_;
  std::vector<int>* calo_side_;
  std::vector<int>* calo_column_;
  std::vector<int>* calo_row_;
  std::vector<int>* calo_wall_;
  std::vector<double>* calo_loc_x_;
  std::vector<double>* calo_loc_y_;
  std::vector<double>* calo_loc_z_;
} particleeventstorage;


class PTD2Root : public dpp::base_module {
 public:
  //! Construct module
  PTD2Root();

  //! Destructor
  virtual ~PTD2Root();

  //! Configure the module
  virtual void initialize(const datatools::properties& myConfig,
			  datatools::service_manager& flServices,
			  dpp::module_handle_dict_type& what);

  //! Process supplied data record
  virtual dpp::base_module::process_status process(datatools::things& workItem);

  //! Reset the module
  virtual void reset();

 private:
  // configurable data member
  std::string filename_output_;

  // geometry service
  const geomtools::manager* geometry_manager_; //!< The geometry manager

  // Variables
  TFile* hfile_;
  TTree* tree_;
  
  HeaderEventStorage header_;
  ParticleEventStorage particle_;

  // Macro which automatically creates the interface needed
  // to enable the module to be loaded at runtime
  DPP_MODULE_REGISTRATION_INTERFACE(PTD2Root);
};
#endif // PTD2ROOT_HH
