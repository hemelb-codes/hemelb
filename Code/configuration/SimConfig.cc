// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <string>

#include "configuration/SimConfig.h"
#include "reporting/BuildInfo.h"

namespace hemelb::configuration
{
    // Base IC
    ICConfigBase::ICConfigBase(std::optional<LatticeTimeStep> t) : t0(t) {
    }

    // Uniform equilibrium IC
    EquilibriumIC::EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p) : ICConfigBase(t), p_mmHg(p), v_ms(0.0) {
    }

    EquilibriumIC::EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v) : ICConfigBase(t), p_mmHg(p), v_ms(v) {
    }

    // checkpoint IC
    CheckpointIC::CheckpointIC(std::optional<LatticeTimeStep> t, std::filesystem::path cp, std::optional<std::filesystem::path> maybeOff)
            : ICConfigBase(t), cpFile(std::move(cp)), maybeOffFile(std::move(maybeOff)) {
    }


    std::unique_ptr<SimConfig> SimConfig::New(const path& path)
    {
      auto ans = std::unique_ptr<SimConfig>{new SimConfig{path}};
      ans->Init();
      return ans;
    }


    SimConfig::SimConfig(const path& path) :
        xmlFilePath(path)
    {
    }

    void SimConfig::Init()
    {
      if (!std::filesystem::exists(xmlFilePath))
      {
        throw Exception() << "Config file '" << xmlFilePath << "' does not exist";
      }
      auto rawXmlDoc = io::xml::Document(xmlFilePath);
      DoIO(rawXmlDoc.GetRoot());
    }

    // Turn an input XML-relative path into a full path
    std::filesystem::path SimConfig::RelPathToFullPath(std::string_view path) const {
        auto xml_dir = xmlFilePath.parent_path();
        return std::filesystem::absolute(xml_dir / path);
    }

    void SimConfig::DoIO(io::xml::Element topNode)
    {
      // Top element must be:
      // <hemelbsettings version="5" />
      if (topNode.GetName() != "hemelbsettings")
        throw Exception() << "Invalid root element: " << topNode.GetPath();

      auto version = topNode.GetAttributeOrThrow<unsigned>("version");
      if (version != 5U)
        throw Exception() << "Unrecognised XML version. Expected 5, got " << version;

      DoIOForSimulation(topNode.GetChildOrThrow("simulation"));

      DoIOForGeometry(topNode.GetChildOrThrow("geometry"));

      if (topNode.GetChildOrNull("colloids") != io::xml::Element::Missing())
      {
        hasColloidSection = true;
      }

      DoIOForInitialConditions(topNode.GetChildOrThrow("initialconditions"));

      inlets = DoIOForInOutlets(topNode.GetChildOrThrow("inlets"));
      outlets = DoIOForInOutlets(topNode.GetChildOrThrow("outlets"));

      // Optional element <properties>
      if (auto propertiesEl = topNode.GetChildOrNull("properties"))
        DoIOForProperties(propertiesEl);

      // Optional element <monitoring>
      if (auto monitoringEl = topNode.GetChildOrNull("monitoring"))
        DoIOForMonitoring(monitoringEl);

      // The RBC section must be parsed *after* the inlets and outlets have been
      // defined
      if (auto rbcEl = topNode.GetChildOrNull("redbloodcells")) {
#ifdef HEMELB_BUILD_RBC
	//rbcConf = new redblood::RBCConfig;
          rbcConf = DoIOForRedBloodCells(rbcEl);
#else
	      throw Exception() << "Input XML has redbloodcells section but HEMELB_BUILD_RBC=OFF";
#endif
      }
    }

    void SimConfig::DoIOForSimulation(const io::xml::Element simEl)
    {
        // Required element
        // <stresstype value="enum lb::StressTypes" />
        sim_info.stress_type = [](unsigned v) {
            switch (v) {
                case lb::IgnoreStress:
                    return lb::IgnoreStress;
                case lb::ShearStress:
                    return lb::IgnoreStress;
                case lb::VonMises:
                    return lb::IgnoreStress;
                default:
                    throw Exception() << "Invalid stresstype: " << v;
            }
        }(simEl.GetChildOrThrow("stresstype").GetAttributeOrThrow<unsigned>("value"));

        // Required element
        // <steps value="unsigned" units="lattice />
        const io::xml::Element stepsEl = simEl.GetChildOrThrow("steps");
        GetDimensionalValue(stepsEl, "lattice", sim_info.time.total_steps);

        // Required element
        // <step_length value="float" units="s" />
        const io::xml::Element tsEl = simEl.GetChildOrThrow("step_length");
        GetDimensionalValue(tsEl, "s", sim_info.time.step_s);

        // Optional element
        // <extra_warmup_steps value="unsigned" units="lattice" />
        if (auto wuEl = simEl.GetChildOrNull("extra_warmup_steps"))
        {
            GetDimensionalValue(wuEl, "lattice", sim_info.time.warmup_steps);
            sim_info.time.total_steps += sim_info.time.warmup_steps;
        } else {
            sim_info.time.warmup_steps = 0;
        }

        // Required element
        // <voxel_size value="float" units="m" />
        const io::xml::Element vsEl = simEl.GetChildOrThrow("voxel_size");
        GetDimensionalValue(vsEl, "m", sim_info.space.step_m);

        // Required element
        // <origin value="(x,y,z)" units="m" />
        const io::xml::Element originEl = simEl.GetChildOrThrow("origin");
        GetDimensionalValue(originEl, "m", sim_info.space.geometry_origin_m);

        // Optional element
        // <fluid_density value="float" units="kg/m3" />
        sim_info.fluid.density_kgm3 = simEl.GetChildOrNull("fluid_density").transform(
                [](io::xml::Element const& el) {
                    return GetDimensionalValue<PhysicalDensity>(el, "kg/m3");
                }).value_or(DEFAULT_FLUID_DENSITY_Kg_per_m3);

        // Optional element
        // <fluid_viscosity value="float" units="Pa.s" />
        sim_info.fluid.viscosity_Pas = simEl.GetChildOrNull("fluid_viscosity").transform(
                [](io::xml::Element const& el) {
                    return GetDimensionalValue<PhysicalDynamicViscosity>(el, "Pa.s");
                }).value_or(DEFAULT_FLUID_VISCOSITY_Pas);

        // Optional element (default = 0)
        // <reference_pressure value="float" units="mmHg" />
        sim_info.fluid.reference_pressure_mmHg = simEl.GetChildOrNull("reference_pressure").transform(
                [](io::xml::Element const& el) {
                    return GetDimensionalValue<PhysicalPressure>(el, "mmHg");
                }).value_or(0);
    }

    void SimConfig::DoIOForGeometry(const io::xml::Element geometryEl)
    {
      // Required element
      // <geometry>
      //  <datafile path="relative path to GMY" />
      // </geometry>
      dataFilePath = RelPathToFullPath(geometryEl.GetChildOrThrow("datafile").GetAttributeOrThrow("path"));
    }

    /**
     * Helper function to ensure that the iolet being created matches the compiled
     * iolet BC.
     * @param ioletEl
     * @param requiredBC
     */
    void SimConfig::CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                           const std::string& requiredBC) const
    {
      // Check that HEMELB_*LET_BOUNDARY is consistent with this
      const std::string& ioletTypeName = ioletEl.GetName();
      std::string hemeIoletBC;

      if (ioletTypeName == "inlet")
        hemeIoletBC = reporting::inlet_boundary_condition;
      else if (ioletTypeName == "outlet")
        hemeIoletBC = reporting::outlet_boundary_condition;
      else
        throw Exception() << "Unexpected element name '" << ioletTypeName
            << "'. Expected 'inlet' or 'outlet'";

      if (requiredBC != hemeIoletBC)
      {
        throw Exception() << "XML configuration for " << ioletTypeName << " (line "
            << ioletEl.GetLine()
            << ") not consistent with compile-time choice of boundary condition '" << hemeIoletBC
            << "'";
      }
    }

    auto SimConfig::DoIOForInOutlets(const io::xml::Element ioletsEl) const -> std::vector<IoletConfig>
    {
      const std::string& nodeName = ioletsEl.GetName();

      const std::string childNodeName = nodeName.substr(0, nodeName.size() - 1);
      std::vector<IoletConfig> ioletList;
      for (io::xml::Element currentIoletNode = ioletsEl.GetChildOrNull(childNodeName);
          currentIoletNode != io::xml::Element::Missing();
          currentIoletNode = currentIoletNode.NextSiblingOrNull(childNodeName))
      {
        // Determine which InOutlet to create
        io::xml::Element conditionEl = currentIoletNode.GetChildOrThrow("condition");
        auto conditionType = conditionEl.GetAttributeOrThrow("type");

        IoletConfig newIolet;

        if (conditionType == "pressure")
        {
          newIolet = DoIOForPressureInOutlet(currentIoletNode);
        }
        else if (conditionType == "velocity")
        {
          newIolet = DoIOForVelocityInOutlet(currentIoletNode);
        }
        else
        {
          throw Exception() << "Invalid boundary condition type '" << conditionType << "' in "
              << conditionEl.GetPath();
        }
        ioletList.push_back(std::move(newIolet));
      }
      return ioletList;
    }

    auto SimConfig::DoIOForPressureInOutlet(const io::xml::Element& ioletEl) const -> IoletConfig
    {
      CheckIoletMatchesCMake(ioletEl, "NASHZEROTHORDERPRESSUREIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      auto conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      if (conditionSubtype == "cosine")
      {
        return DoIOForCosinePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        return DoIOForFilePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "multiscale")
      {
        return DoIOForMultiscalePressureInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }
    }

    auto SimConfig::DoIOForVelocityInOutlet(const io::xml::Element& ioletEl) const -> IoletConfig
    {
      CheckIoletMatchesCMake(ioletEl, "LADDIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      auto conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      if (conditionSubtype == "parabolic")
      {
        return DoIOForParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "womersley")
      {
        return DoIOForWomersleyVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        return DoIOForFileVelocityInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }
    }

    void SimConfig::DoIOForProperties(const io::xml::Element& propertiesEl)
    {
      for (io::xml::ChildIterator poPtr = propertiesEl.IterChildren("propertyoutput");
          !poPtr.AtEnd(); ++poPtr)
      {
        propertyOutputs.push_back(DoIOForPropertyOutputFile(*poPtr));
      }

      if (auto cpEl = propertiesEl.GetChildOrNull("checkpoint")) {
	// Create a checkpoint property extractor.
	//
	// This is just a normal one, but fixed to be whole geometry,
	// only distributions, at double precision.
	//
	// Get stuff from XML
	auto file = extraction::PropertyOutputFile{};
	file.filename = cpEl.GetAttributeOrThrow("file");
	cpEl.GetAttributeOrThrow("period", file.frequency);
	// Configure the file
	file.geometry.reset(new extraction::WholeGeometrySelector());
	file.ts_mode = extraction::single_timestep_files{};

	extraction::OutputField field;
	field.name = "distributions";
	field.noffsets = 0;
	field.offset = {};
	field.typecode = distribn_t{0};
	field.src = extraction::source::Distributions{};
	file.fields.push_back(field);
	// Add to outputs
	propertyOutputs.push_back(std::move(file));
      }
    }

    extraction::PropertyOutputFile SimConfig::DoIOForPropertyOutputFile(
        const io::xml::Element& propertyoutputEl)
    {
      auto file = extraction::PropertyOutputFile{};

      auto&& ts_mode = propertyoutputEl.GetAttributeMaybe("timestep_mode").value_or("multi");
      if (ts_mode == "multi") {
	file.ts_mode = extraction::multi_timestep_file{};
      } else if (ts_mode == "single") {
	file.ts_mode = extraction::single_timestep_files{};
      } else {
	throw Exception()
	  << "Invalid value of timestep_mode attribute '" << ts_mode
	  << "' at: " << propertyoutputEl.GetPath();
      }

      file.filename = propertyoutputEl.GetAttributeOrThrow("file");
      if (std::holds_alternative<extraction::single_timestep_files>(file.ts_mode)) {
	auto const& p = file.filename.native();
	// Path must contain exactly one printf conversion specifier
	// for an integer.
	auto const errmsg = "For single timestep output files, "
	  "the path must contain exactly one '%d' and no other '%' characters";
	auto n_pc = std::count(p.begin(), p.end(), '%');
	auto i_pcd = p.find("%d", 0, 2);
	if (n_pc != 1 || i_pcd == std::string::npos) {
	  throw Exception() << errmsg;
	}
      }

      propertyoutputEl.GetAttributeOrThrow("period", file.frequency);

      io::xml::Element geometryEl = propertyoutputEl.GetChildOrThrow("geometry");
      auto type = geometryEl.GetAttributeOrThrow("type");

      if (type == "plane")
      {
        file.geometry.reset(DoIOForPlaneGeometry(geometryEl));
      }
      else if (type == "line")
      {
        file.geometry.reset(DoIOForLineGeometry(geometryEl));
      }
      else if (type == "whole")
      {
        file.geometry.reset(new extraction::WholeGeometrySelector());
      }
      else if (type == "surface")
      {
        file.geometry.reset(new extraction::GeometrySurfaceSelector());
      }
      else if (type == "surfacepoint")
      {
        file.geometry.reset(DoIOForSurfacePoint(geometryEl));
      }
      else
      {
        throw Exception() << "Unrecognised property output geometry selector '" << type
            << "' in element " << geometryEl.GetPath();
      }

      for (io::xml::ChildIterator fieldPtr = propertyoutputEl.IterChildren("field");
          !fieldPtr.AtEnd(); ++fieldPtr)
        file.fields.push_back(DoIOForPropertyField(*fieldPtr));

      return file;
    }

    extraction::StraightLineGeometrySelector* SimConfig::DoIOForLineGeometry(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element point1El = geometryEl.GetChildOrThrow("point");
      io::xml::Element point2El = point1El.NextSiblingOrThrow("point");

      PhysicalPosition point1;
      PhysicalPosition point2;

      GetDimensionalValue(point1El, "m", point1);
      GetDimensionalValue(point2El, "m", point2);

      return new extraction::StraightLineGeometrySelector(point1.as<float>(), point2.as<float>());
    }

    extraction::PlaneGeometrySelector* SimConfig::DoIOForPlaneGeometry(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");
      io::xml::Element normalEl = geometryEl.GetChildOrThrow("normal");

      PhysicalPosition point;
      util::Vector3D<float> normal;

      GetDimensionalValue(pointEl, "m", point);
      GetDimensionalValue(normalEl, "dimensionless", normal);

      auto radius = geometryEl.GetChildOrNull("radius").transform(
	[](io::xml::Element const& el) {
	  return GetDimensionalValue<PhysicalDistance>(el, "m");
	});
      if (radius)
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal, *radius);
      else	
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal);
    }

    extraction::SurfacePointSelector* SimConfig::DoIOForSurfacePoint(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");

      PhysicalPosition point;
      GetDimensionalValue(pointEl, "m", point);
      return new extraction::SurfacePointSelector(point.as<float>());
    }

    extraction::OutputField SimConfig::DoIOForPropertyField(const io::xml::Element& fieldEl)
    {
      extraction::OutputField field;
      auto type = fieldEl.GetAttributeOrThrow("type");
      // Default name is identical to type.
      field.name = fieldEl.GetAttributeMaybe("name").value_or(type);

      // Default offset is none
      field.noffsets = 0;
      field.offset = {};

      // The default type to be written is float
      field.typecode = float{0.0};

      // Check and assign the type.
      if (type == "pressure")
      {
        field.src = extraction::source::Pressure{};
	// Pressure has an offset of the reference pressure
	field.noffsets = 1;
	field.offset = {sim_info.fluid.reference_pressure_mmHg};
      }
      else if (type == "velocity")
      {
        field.src = extraction::source::Velocity{};
      }
      else if (type == "vonmisesstress")
      {
        field.src = extraction::source::VonMisesStress{};
      }
      else if (type == "shearstress")
      {
        field.src = extraction::source::ShearStress{};
      }
      else if (type == "shearrate")
      {
        field.src = extraction::source::ShearRate{};
      }
      else if (type == "stresstensor")
      {
        field.src = extraction::source::StressTensor{};
      }
      else if (type == "traction")
      {
        field.src = extraction::source::Traction{};
      }
      else if (type == "tangentialprojectiontraction")
      {
        field.src = extraction::source::TangentialProjectionTraction{};
      }
      else if (type == "distributions")
      {
        field.src = extraction::source::Distributions{};
      }
      else if (type == "mpirank")
      {
        field.src = extraction::source::MpiRank{};
	field.typecode = int{0};
      }
      else
      {
        throw Exception() << "Unrecognised field type '" << type << "' in " << fieldEl.GetPath();
      }
      return field;
    }

    void SimConfig::DoIOForBaseInOutlet(const io::xml::Element& ioletEl,
                                        IoletConfigBase& ioletConf) const
    {
        io::xml::Element positionEl = ioletEl.GetChildOrThrow("position");
        io::xml::Element normalEl = ioletEl.GetChildOrThrow("normal");

        GetDimensionalValue(positionEl, "m", ioletConf.position);
        GetDimensionalValue(normalEl, "dimensionless", ioletConf.normal);

        if (sim_info.time.warmup_steps) {
            ioletConf.warmup_steps = sim_info.time.warmup_steps;
        }

        // Optional element <flowextension>
        if (auto const flowEl = ioletEl.GetChildOrNull("flowextension")) {
            FlowExtensionConfig flowConf;
            GetDimensionalValue(flowEl.GetChildOrThrow("length"), "m", flowConf.length_m);
            GetDimensionalValue(flowEl.GetChildOrThrow("radius"), "m", flowConf.radius_m);

            // Optional element - default is length of flowext
            // <fadelength ioletConf="float" units="m" />
            flowConf.fadelength_m = flowEl.GetChildOrNull("fadelength").transform(
                    [](io::xml::Element const &el) {
                        return GetDimensionalValue<PhysicalDistance>(el, "m");
                    }).value_or(flowConf.length_m);

            // Infer normal and position from inlet
            // However, normals point in *opposite* direction, and, as a flowConf, origin are at opposite
            // end of the cylinder
            flowConf.normal = -ioletConf.normal;
            flowConf.normal.Normalise();

            flowConf.origin_m = ioletConf.position - flowConf.normal * flowConf.length_m;
            ioletConf.flow_extension = flowConf;
        }

        // Optional element(s) <insertcell>
        std::int64_t horrible_hack_seed = 0;
        for (auto insertEl = ioletEl.IterChildren("insertcell"); !insertEl.AtEnd(); ++insertEl) {
            CellInserterConfig inserterConf;
            auto maybeSeed = insertEl->GetChildOrNull("seed").transform(
                    [](io::xml::Element const& _) { return _.GetAttributeOrThrow<std::int64_t>("value"); }
            );
            if (maybeSeed) {
                inserterConf.seed = *maybeSeed;
            } else {
                // We need to seed each of the RBCInserterWithPerturbation objects consistently across MPI processes
                if (horrible_hack_seed == 0) {
                    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
                    auto comm_world = hemelb::net::MpiCommunicator::World();
                    comm_world.Broadcast(seed, 0);
                    std::stringstream message;
                    message << "RBC insertion random seed: " << std::hex << std::showbase << seed;
                    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>(message.str());
                    horrible_hack_seed = seed;
                } else {
                    horrible_hack_seed += 1;
                }
                inserterConf.seed = horrible_hack_seed;
            }
            inserterConf.template_name = insertEl->GetAttributeOrThrow("template");
//                    if (templateCells.count(templateName) == 0)
//                    {
//                        throw Exception() << "Template cell name does not match a known template cell";
//                    }

            // Rotate cell to align z axis with given position, and then z axis with flow
            // If phi == 0, then cell symmetry axis is aligned with the flow
            inserterConf.theta_rad = GetDimensionalValueWithDefault<Angle>(*insertEl, "theta", "rad", 0e0);
            inserterConf.phi_rad = GetDimensionalValueWithDefault<Angle>(*insertEl, "phi", "rad", 0e0);
            inserterConf.translation_m = {
                    GetDimensionalValueWithDefault<PhysicalDistance>(*insertEl, "x", "m", 0e0),
                    GetDimensionalValueWithDefault<PhysicalDistance>(*insertEl, "y", "m", 0e0),
                    GetDimensionalValueWithDefault<PhysicalDistance>(*insertEl, "z", "m", 0e0)
            };

            inserterConf.offset_s = GetDimensionalValueWithDefault<PhysicalTime>(*insertEl,
                                                                                 "offset",
                                                                                 "s",
                                                                                 0);
            inserterConf.drop_period_s = GetDimensionalValue<LatticeTime>(
                    insertEl->GetChildOrThrow("every"), "s");
            inserterConf.dt_s = GetDimensionalValueWithDefault<LatticeTime>(
                    *insertEl, "delta_t", "s", 0e0);

            inserterConf.dtheta_rad = GetDimensionalValueWithDefault<Angle>(*insertEl,
                                                                      "delta_theta",
                                                                      "rad",
                                                                      0e0);
            inserterConf.dphi_rad = GetDimensionalValueWithDefault<Angle>(*insertEl,
                                                                    "delta_phi",
                                                                    "rad",
                                                                    0e0);
            inserterConf.dx_m = GetDimensionalValueWithDefault<LatticeDistance>(*insertEl,
                                                                            "delta_x",
                                                                            "m",
                                                                            0e0);
            inserterConf.dy_m = GetDimensionalValueWithDefault<LatticeDistance>(*insertEl,
                                                                            "delta_y",
                                                                            "m",
                                                                            0e0);
            ioletConf.cell_inserters.push_back(inserterConf);
        }
    }

    template <typename T, typename F>
    auto opt_transform (std::optional<T>&& o, F&& f) {
        using U = std::remove_cv_t<std::invoke_result_t<F, T>>;
        if (o.has_value())
            return std::make_optional<U>(
                    std::invoke(
                            std::forward<F>(f),
                            std::move(o.value())
                    )
            );
        return std::optional<U>{};
    }

    void SimConfig::DoIOForInitialConditions(io::xml::Element initialconditionsEl)
    {
        // The <time> element may be present - if so, it will set the
        // initial timestep value
        auto t0 = initialconditionsEl.GetChildOrNull("time").transform([](auto el) {
            return GetDimensionalValue<LatticeTimeStep>(el, "lattice");
        });

        // Exactly one of {<pressure>, <checkpoint>} must be present
        // TODO: use something other than an if-tree
        auto pressureEl = initialconditionsEl.GetChildOrNull("pressure");
        auto checkpointEl = initialconditionsEl.GetChildOrNull("checkpoint");
        if (pressureEl) {
            if (checkpointEl) {
                // Both are present - this is an error
                throw Exception()
                        << "XML contains both <pressure> and <checkpoint> sub elements of <initialconditions>";
            } else {
                // Only pressure
                io::xml::Element uniformEl = pressureEl.GetChildOrThrow("uniform");
                PhysicalPressure p0_mmHg;
                GetDimensionalValue(uniformEl, "mmHg", p0_mmHg);
                initial_condition = EquilibriumIC(t0, p0_mmHg);
            }
        } else {
            if (checkpointEl) {
                // Only checkpoint
                initial_condition = CheckpointIC(
                        t0,
                        RelPathToFullPath(checkpointEl.GetAttributeOrThrow("file")),
                        opt_transform(
                                checkpointEl.GetAttributeMaybe("offsets"),
                                [&](std::string_view sv) { return RelPathToFullPath(sv); }
                        )
                );
            } else {
                // No IC!
                throw Exception() << "XML <initialconditions> element contains no known initial condition type";
            }
        }
    }

    auto SimConfig::DoIOForCosinePressureInOutlet(
        const io::xml::Element& ioletEl) const -> IoletConfig
    {
      //auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetCosine>();
      CosinePressureIoletConfig newIolet;
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      newIolet.amp_mmHg = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("amplitude"), "mmHg");
      newIolet.mean_mmHg = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("mean"), "mmHg");
      newIolet.phase_rad = GetDimensionalValue<Angle>(conditionEl.GetChildOrThrow("phase"), "rad");
      newIolet.period_s = GetDimensionalValue<LatticeTime>(conditionEl.GetChildOrThrow("period"), "s");
      return newIolet;
    }

    auto SimConfig::DoIOForFilePressureInOutlet(
        const io::xml::Element& ioletEl) const -> IoletConfig
    {
      FilePressureIoletConfig newIolet;
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element pathEl = ioletEl.GetChildOrThrow("condition").GetChildOrThrow("path");
      newIolet.file_path = RelPathToFullPath(pathEl.GetAttributeOrThrow("value"));

      return newIolet;
    }

    auto SimConfig::DoIOForMultiscalePressureInOutlet(
        const io::xml::Element& ioletEl) const -> IoletConfig
    {
      MultiscalePressureIoletConfig newIolet;
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element pressureEl = conditionEl.GetChildOrThrow("pressure");
      GetDimensionalValue(pressureEl, "mmHg", newIolet.pressure_reference_mmHg);

      const io::xml::Element velocityEl = conditionEl.GetChildOrThrow("velocity");
      GetDimensionalValue(velocityEl, "m/s", newIolet.velocity_reference_ms);

      newIolet.label = conditionEl.GetChildOrThrow("label").GetAttributeOrThrow("value");
      return newIolet;
    }

    auto SimConfig::DoIOForParabolicVelocityInOutlet(
        const io::xml::Element& ioletEl) const -> IoletConfig
    {
        ParabolicVelocityIoletConfig newIolet;
        DoIOForBaseInOutlet(ioletEl, newIolet);

        const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);
        GetDimensionalValue(conditionEl.GetChildOrThrow("maximum"), "m/s", newIolet.max_speed_ms);
        return newIolet;
    }

    auto SimConfig::DoIOForWomersleyVelocityInOutlet(
        const io::xml::Element& ioletEl) const -> IoletConfig
    {
        WomersleyVelocityIoletConfig newIolet;
        DoIOForBaseInOutlet(ioletEl, newIolet);

        const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);
        GetDimensionalValue(conditionEl.GetChildOrThrow("pressure_gradient_amplitude"),
                            "mmHg/m", newIolet.pgrad_amp_mmHgm);
        GetDimensionalValue(conditionEl.GetChildOrThrow("period"), "s", newIolet.period_s);
        GetDimensionalValue(conditionEl.GetChildOrThrow("womersley_number"), "dimensionless", newIolet.womersley);
        return newIolet;
    }

    auto SimConfig::DoIOForFileVelocityInOutlet(
            const io::xml::Element& ioletEl) const -> IoletConfig
    {
        FileVelocityIoletConfig newIolet;
        DoIOForBaseInOutlet(ioletEl, newIolet);

        const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

        newIolet.file_path = RelPathToFullPath(conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value"));
        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);

        return newIolet;
    }

    bool SimConfig::HasColloidSection() const
    {
      return hasColloidSection;
    }

    void SimConfig::DoIOForMonitoring(const io::xml::Element& monEl)
    {
      io::xml::Element convEl = monEl.GetChildOrNull("steady_flow_convergence");
      if (convEl != io::xml::Element::Missing())
      {
        DoIOForSteadyFlowConvergence(convEl);
      }

      monitoringConfig.doIncompressibilityCheck = (monEl.GetChildOrNull("incompressibility")
          != io::xml::Element::Missing());
    }

    void SimConfig::DoIOForSteadyFlowConvergence(const io::xml::Element& convEl)
    {
      monitoringConfig.doConvergenceCheck = true;
      monitoringConfig.convergenceRelativeTolerance = convEl.GetAttributeOrThrow<double>("tolerance");
      monitoringConfig.convergenceTerminate = (convEl.GetAttributeOrThrow("terminate") == "true");

      if (convEl.IterChildren("criterion").AtEnd())
      {
        throw Exception() << "At least one convergence criterion must be provided in "
            << convEl.GetPath();
      }

      for (io::xml::ChildIterator criteriaIt = convEl.IterChildren("criterion");
          !criteriaIt.AtEnd(); ++criteriaIt)
      {
        DoIOForConvergenceCriterion(*criteriaIt);
      }
    }

    void SimConfig::DoIOForConvergenceCriterion(const io::xml::Element& criterionEl)
    {
      auto criterionType = criterionEl.GetAttributeOrThrow("type");

      // We only allow velocity-based convergence check for the time being
      if (criterionType != "velocity")
      {
        throw Exception() << "Invalid convergence criterion type " << criterionType << " in "
            << criterionEl.GetPath();
      }
      monitoringConfig.convergenceVariable = extraction::source::Velocity{};
      monitoringConfig.convergenceReferenceValue = GetDimensionalValue<PhysicalSpeed>(criterionEl, "m/s");
    }

    TemplateCellConfig SimConfig::readCell(const io::xml::Element& cellNode) const {
        TemplateCellConfig ans;
        ans.name = cellNode.GetAttributeMaybe("name").value_or("default");
        auto const shape = cellNode.GetChildOrThrow("shape");
        ans.mesh_path = RelPathToFullPath(
                shape.GetAttributeOrThrow("mesh_path")
        );

        auto get_fmt = [](io::xml::Element const& shape, std::string attr_name) -> MeshFormat {
            auto fmt = shape.GetAttributeOrThrow(attr_name);
            if (fmt == "VTK") {
                return VTKMeshFormat{};
            } else if (fmt == "Krueger") {
                log::Logger::Log<log::Warning, log::Singleton>("Krueger format meshes are deprecated, move to VTK when you can.");
                return KruegerMeshFormat{};
            }
            throw Exception() << "Invalid " << attr_name << " '" << fmt << "' on element " << shape.GetPath();
        };

        ans.format = get_fmt(shape, "mesh_format");

        ans.scale_m = GetDimensionalValue<PhysicalDistance>(cellNode.GetChildOrThrow("scale"), "m");
        auto reference_mesh_path = shape.GetAttributeMaybe("reference_mesh_path");
        if (reference_mesh_path) {
            ans.reference_mesh_path = RelPathToFullPath(*reference_mesh_path);
            ans.reference_mesh_format = get_fmt(shape, "reference_mesh_format");
        }

        auto const moduliNode = cellNode.GetChildOrNull("moduli");
        ans.moduli.bending_Nm = GetDimensionalValueWithDefault<PhysicalModulus>(
                moduliNode, "bending", "Nm", 2e-19);
        ans.moduli.surface_lat = GetDimensionalValueWithDefault<LatticeModulus>(
                moduliNode, "surface", "lattice", 1e0);
        ans.moduli.volume_lat = GetDimensionalValueWithDefault<LatticeModulus>(
                moduliNode, "volume", "lattice", 1e0);
        ans.moduli.dilation_lat = GetDimensionalValueWithDefault<LatticeModulus>(
                moduliNode, "dilation", "lattice", 0.75);
        if (1e0 < ans.moduli.dilation_lat or ans.moduli.dilation_lat < 0.5)
        {
            log::Logger::Log<log::Critical, log::Singleton>("Dilation modulus is outside the recommended range 1e0 >= m >= 0.5");
        }
        ans.moduli.strain_Npm = GetDimensionalValueWithDefault<PhysicalModulus>(
                moduliNode, "strain", "N/m", 5e-6);
        return ans;
    }

    std::map<std::string, TemplateCellConfig> SimConfig::readTemplateCells(io::xml::Element const& cellsEl) const {
        std::map<std::string, TemplateCellConfig> ans;
        for (auto cellNode = cellsEl.GetChildOrNull("cell");
             cellNode;
             cellNode = cellNode.NextSiblingOrNull("cell"))
        {
            auto const key = std::string{cellNode.GetAttributeMaybe("name").value_or("default")};
            if (ans.contains(key))
                throw Exception() << "Multiple template mesh with same name: " << key;

            ans[key] = readCell(cellNode);
        }
        return ans;
    }

    NodeForceConfig readNode2NodeForce(const io::xml::Element& node) {
        auto intensity = GetDimensionalValueWithDefault<PhysicalEnergy>(
                node, "intensity", {"Nm", "lattice"}, {1.0, "Nm"});

        auto cutoffdist = GetDimensionalValueWithDefault<LatticeDistance>(
                node, "cutoffdistance", "lattice", 1.0);

        // exponent doesnt have units apparently
        auto exponent = node.and_then(
                [](io::xml::Element const& el) { return el.GetChildOrNull("exponent"); }
        ).transform(
                [](io::xml::Element const& el) { return el.GetAttributeOrThrow<std::size_t>("value"); }
        ).value_or(std::size_t(2));

        return {intensity.first, intensity.second, cutoffdist, exponent};
    }

    RBCConfig SimConfig::DoIOForRedBloodCells(const io::xml::Element &rbcEl) const {
        RBCConfig ans;

        const io::xml::Element controllerNode = rbcEl.GetChildOrThrow("controller");
        ans.boxSize = GetDimensionalValue<LatticeDistance>(controllerNode.GetChildOrThrow("boxsize"), "lattice");

        if (auto cellsEl = rbcEl.GetChildOrNull("cells"))
            ans.meshes = readTemplateCells(rbcEl.GetChildOrNull("cells"));

        // Now we can check that the cell inserters only refer to templates that exist
        auto cell_inserter_templates_exist = [&] (IoletConfig const& iolet_v) {
            bool ok = true;
            std::visit([&](auto const& iolet) {
                if constexpr(std::is_same_v<std::decay_t<decltype(iolet)>, std::monostate>) {
                    throw Exception() << "invalid iolet";
                } else {
                    for (auto const &ins: iolet.cell_inserters) {
                        if (!ans.meshes.contains(ins.template_name)) {
                            ok = false;
                            log::Logger::Log<log::Error, log::Singleton>(
                                    "Cell inserter for iolet refers to template that does not exist '%s'",
                                    ins.template_name.c_str());
                        }
                    }
                }
            }, iolet_v);
            return ok;
        };
        if (!(std::all_of(inlets.begin(), inlets.end(), cell_inserter_templates_exist)
              &&
              std::all_of(outlets.begin(), outlets.end(), cell_inserter_templates_exist))) {
            throw Exception() << "One or more cell inserters refer to a template that does not exist";
        }

        ans.cell2cell = readNode2NodeForce(rbcEl.GetChildOrNull("cell2Cell"));
        ans.cell2wall = readNode2NodeForce(rbcEl.GetChildOrNull("cell2Wall"));
        if (ans.boxSize < ans.cell2wall.cutoffdist)
            throw Exception() << "Box-size < cell-wall interaction size: "
                                 "cell-wall interactions cannot be all accounted for.";

        if (ans.boxSize < ans.cell2cell.cutoffdist)
            throw Exception() << "Box-size < cell-cell interaction size: "
                                 "cell-cell interactions cannot be all accounted for.";

        ans.output_period = GetDimensionalValue<LatticeTimeStep>(
                rbcEl.GetChildOrThrow("output").GetChildOrThrow("period"),
                "lattice"
        );

        return ans;
    }

    const MonitoringConfig& SimConfig::GetMonitoringConfiguration() const
    {
      return monitoringConfig;
    }

}
