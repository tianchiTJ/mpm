//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPMBase<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::solve() {
  bool status = true;

  // Get analysis type USL/USF
  if (io_->analysis_type() == "MPMExplicitUSL2D" ||
      io_->analysis_type() == "MPMExplicitUSL3D")
    this->usl_ = true;

  console_->error("Analysis{} {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Phase
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Test if strain energy is needed to compute
  bool strain_energy = false;
  if (analysis_.find("strain_energy") != analysis_.end())
    strain_energy = analysis_["strain_energy"].template get<bool>();

  // Test if reference coordinates is needed to be assigned
  mpm::Index reference_step = 0;
  if (analysis_.find("reference_step") != analysis_.end())
    reference_step = analysis_["reference_step"].template get<mpm::Index>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) status = false;

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) status = false;

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) status = false;

  // Assign material to particles
  // Get particle properties
  auto particle_props = io_->json_object("particle");
  // Material id
  const auto material_id =
      particle_props["material_id"].template get<unsigned>();

  // Get material from list of materials
  auto material = materials_.at(material_id);

  // Iterate over each particle to assign material
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, phase, material));

  // Assign material to particle sets
  if (particle_props["particle_sets"].size() != 0) {
    // Assign material to particles in the specific sets
    bool set_material_status = this->apply_properties_to_particles_sets();
  }

  // Compute mass
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));

  // Add new particle
  // Status of add particle
  bool add_particle = false;
  if (particle_props.find("add_particle") != particle_props.end())
    add_particle =
        particle_props["add_particle"]["add_particle"].template get<bool>();
  // Properties of new particle
  // Start step of adding particle
  mpm::Index apstep_start;
  // End step of adding particle
  mpm::Index apstep_end;
  // Step interval of adding particle
  mpm::Index apstep_inv;
  // Number of new particles in one adding step
  unsigned ap_number;
  // Start id of new particles
  mpm::Index start_id;
  // Material id of new particles
  unsigned new_particle_mid;
  // Initial volume of new particle
  double new_particle_volume;
  // Initial coordinates of new particle
  Eigen::Matrix<double, Tdim, 1> new_particle_coordinates;
  // Initial velocity of new particle
  Eigen::Matrix<double, Tdim, 1> new_particle_velocities;
  // Incremental of coordinates
  Eigen::Matrix<double, Tdim, 1> apcoordinates_inv;
  // Initial stress of new particle
  Eigen::Matrix<double, 6, 1> new_particle_stresses;
  // Counter of new particle
  mpm::Index counter_new_particle = 0;

  if (add_particle) {
    // Add particle properties
    auto add_particle_props = particle_props["add_particle"];
    // Assign time properties
    apstep_start =
        add_particle_props["apstep_start"].template get<mpm::Index>();
    apstep_end = add_particle_props["apstep_end"].template get<mpm::Index>();
    apstep_inv = add_particle_props["apstep_inv"].template get<mpm::Index>();
    // Assign start id of new particles
    start_id = add_particle_props["start_id"].template get<mpm::Index>();
    // Assigne initial volume of new particle
    new_particle_volume =
        add_particle_props["new_particle_volume"].template get<double>();
    // Get material from list of materials
    new_particle_mid =
        add_particle_props["new_particle_mid"].template get<unsigned>();
    // Assign number of new particles in one adding step
    ap_number = add_particle_props["ap_number"].template get<unsigned>();
    // Assigne initial and incremental coordinates of new particle
    if (add_particle_props.at("new_particle_coordinates").is_array() &&
        add_particle_props.at("new_particle_coordinates").size() == Tdim &&
        add_particle_props.at("apcoordinates_inv").is_array() &&
        add_particle_props.at("apcoordinates_inv").size() == Tdim) {
      for (unsigned i = 0; i < Tdim; ++i) {
        new_particle_coordinates[i] =
            add_particle_props.at("new_particle_coordinates").at(i);

        apcoordinates_inv[i] = add_particle_props.at("apcoordinates_inv").at(i);
      }
    } else {
      throw std::runtime_error(
          "Specified coordinates of the new particle dimension is invalid");
    }
    // Assigne initial velocities of new particle
    if (add_particle_props.at("new_particle_velocities").is_array() &&
        add_particle_props.at("new_particle_velocities").size() == Tdim) {
      for (unsigned i = 0; i < Tdim; ++i) {
        new_particle_velocities[i] =
            add_particle_props.at("new_particle_velocities").at(i);
      }
    } else {
      throw std::runtime_error(
          "Specified velocities of the new particle dimension is invalid");
    }
    // Assigne initial stress of new particle
    if (add_particle_props.at("new_particle_stress").is_array()) {
      new_particle_stresses[0] =
          add_particle_props.at("new_particle_stress").at(0);
      new_particle_stresses[1] =
          add_particle_props.at("new_particle_stress").at(1);
      new_particle_stresses[3] =
          add_particle_props.at("new_particle_stress").at(3);
      if (Tdim == 3) {
        new_particle_stresses[2] =
            add_particle_props.at("new_particle_stress").at(2);
        new_particle_stresses[4] =
            add_particle_props.at("new_particle_stress").at(4);
        new_particle_stresses[5] =
            add_particle_props.at("new_particle_stress").at(5);
      }
    } else {
      throw std::runtime_error(
          "Specified stress of the new particle dimension is invalid");
    }
  }

  // Check point resume
  if (resume) {
    mpm::Index resume_step = analysis_["resume"]["step"];
    if (resume_step > apstep_start) {
      mpm::Index end_step = resume_step;
      if (apstep_end < resume_step) end_step = apstep_end;
      for (unsigned j = 0;
           j < floor((end_step - apstep_start) / apstep_inv) + 1; ++j) {
        for (unsigned i = 0; i < ap_number; ++i) {
          // New particle id
          mpm::Index new_particle_id = start_id + counter_new_particle;
          // Add new particle
          this->add_new_particle(
              new_particle_id,
              (new_particle_coordinates + apcoordinates_inv * i),
              new_particle_velocities, new_particle_volume,
              new_particle_stresses);
          // Assign material to new particle
          mesh_->assign_new_particle_material(new_particle_id, phase,
                                              materials_.at(new_particle_mid));
          // Counter number of new particle
          counter_new_particle++;
        }
      }
    }
    this->checkpoint_resume();
  }

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    // Add new particle
    if (add_particle && step_ >= apstep_start && step_ < apstep_end &&
        (((step_ - apstep_start) % apstep_inv) == 0)) {
      for (unsigned i = 0; i < ap_number; ++i) {
        // New particle id
        mpm::Index new_particle_id = start_id + counter_new_particle;
        // Add new particle
        this->add_new_particle(
            new_particle_id, (new_particle_coordinates + apcoordinates_inv * i),
            new_particle_velocities, new_particle_volume,
            new_particle_stresses);
        // Assign material to new particle
        mesh_->assign_new_particle_material(new_particle_id, phase,
                                            materials_.at(new_particle_mid));
        // Counter number of new particle
        counter_new_particle++;
      }
    }

    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for initialising nodes and cells
    task_group.run([&] {
      // Apply change material step
      bool change_material_status =
          this->apply_change_material_step(step_, false);

      // Apply remove step
      bool remove_status = mesh_->apply_remove_step(step_);

      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

      // mesh_->find_active_nodes();
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, phase));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, phase),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, phase, std::placeholders::_2));
      // MPI all reduce nodal momentum
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update stress first
    if (!usl_) {
      // Iterate over each particle to calculate strain
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, phase, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, phase, this->dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, phase));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        phase),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, phase, std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, phase));
      }

      // Iterate over each particle to compute stress
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, phase));

      // Iterate over each particle to compute strain energy
      if (strain_energy)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_strain_energy,
                      std::placeholders::_1, phase));
    }

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over each particle to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, phase, this->gravity_));

      // Iterate over each particle to map traction force to nodes
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                    std::placeholders::_1, phase));

      //! Apply nodal tractions
      if (nodal_tractions_) this->apply_nodal_tractions();
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over each particle to compute nodal internal force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                    std::placeholders::_1, phase));
    });
    task_group.wait();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
      // MPI all reduce internal force
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Use nodal velocity to update position
    if (velocity_update_)
      // Iterate over each particle to compute updated position
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position_velocity,
                    std::placeholders::_1, phase, this->dt_));
    else
      // Iterate over each particle to compute updated position
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                    std::placeholders::_1, phase, this->dt_));

    // Update Stress Last
    if (usl_ == true) {
      // Iterate over each particle to calculate strain
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, phase, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, phase, this->dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, phase));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        phase),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, phase, std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, phase));
      }

      // Iterate over each particle to compute stress
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, phase));

      // Iterate over each particle to compute strain energy
      if (strain_energy)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_strain_energy,
                      std::placeholders::_1, phase));
    }

    if (step_ != 0 && step_ == reference_step)
      // Iterate over each particle to assign reference coordinates
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::assign_reference_coordinates,
                    std::placeholders::_1));

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info("Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
                 (this->usl_ ? "USL" : "USF"),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}
