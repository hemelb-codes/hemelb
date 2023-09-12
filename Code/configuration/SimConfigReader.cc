// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/SimConfigReader.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <memory>
#include <optional>
#include <utility>
#include <variant>

#include "build_info.h"
#include "constants.h"
#include "extraction/GeometrySelectors.h"
#include "io/TimePattern.h"
#include "log/Logger.h"
#include "net/MpiCommunicator.h"
#include "util/Vector3D.h"

namespace hemelb::configuration {
    using Element = io::xml::Element;

    // Given an element, get those child elements specified by name (repeats OK)
    // Any missing or extra elements are an error.
    //
    // Probably easier to use the helper template below
    template <std::size_t N>
    auto GetExactChildren(Element const& el, std::array<char const*, N> names) {
        using PAIR = std::pair<std::size_t, std::string_view>;
        std::vector<PAIR> idx_names(N); //= {{0, NAMES}...};
        for (std::size_t i = 0; i < N; ++i) {
            idx_names[i] = {i, names[i]};
        }

        std::array<Element, N> ans;

        // Must inspect all children of the node
        for (auto child: el.Children()) {
            // Check the found child against allowed names
            auto match = std::find_if(idx_names.begin(), idx_names.end(),
                                      [&](PAIR const& i_name) {
                                          return i_name.second == child.GetName();
                                      });
            if (match == idx_names.end())
                throw (Exception() << "Unexpected element found: " << child.GetPath());
            // It was one of the expected, add to output
            ans[match->first] = child;
            // Remove from allowed names
            idx_names.erase(match);
        }
        // Now check for required names not found
        for (auto& [_, name]: idx_names) {
            throw (Exception() << "Required child not found: " << name);
        }
        return ans;
    }
    // Helper for the above.
    template <typename... Ts>
    auto GetExactChildren(Element const& el, Ts&... names) {
        constexpr auto N = sizeof...(Ts);
        return GetExactChildren<N>(el, std::array<char const*, N>{std::forward<Ts>(names)...});
    }

    // Like std::foreach but for the child sublements of the given name, any others being an error.
    template <std::invocable<Element const&> F>
    void ForExactChildren(Element const& el, std::string_view sub_el_name, F&& func) {
        for (auto sub_el: el.Children()) {
            if (auto name = sub_el.GetName(); name != sub_el_name) {
                throw (Exception() << "Unexpected child element: " << sub_el.GetPath());
            }
            std::invoke(func, sub_el);
        }
    }

    SimConfigReader::SimConfigReader(path path) :
            xmlFilePath(std::move(path))
    {
    }

    SimConfig SimConfigReader::Read() const
    {
        if (!std::filesystem::exists(xmlFilePath))
        {
            throw Exception() << "Config file '" << xmlFilePath << "' does not exist";
        }
        auto rawXmlDoc = std::make_shared<io::xml::Document>(xmlFilePath);
        auto ans = DoIO(rawXmlDoc->GetRoot());
        ans.source_xml = rawXmlDoc;
        return ans;
    }

    // Turn an input XML-relative path into a full path
    std::filesystem::path SimConfigReader::RelPathToFullPath(std::string_view path) const {
        auto xml_dir = xmlFilePath.parent_path();
        return std::filesystem::weakly_canonical(xml_dir / path);
    }

    SimConfig SimConfigReader::DoIO(Element topNode) const
    {
        SimConfig ans;
        // Top element must be:
        // <hemelbsettings version="6" />
        constexpr unsigned VERSION = 6;
        if (topNode.GetName() != "hemelbsettings")
            throw Exception() << "Invalid root element: " << topNode.GetPath();

        auto version = topNode.GetAttributeOrThrow<unsigned>("version");
        if (version != VERSION)
            throw Exception() << "Unrecognised XML version. Expected " << VERSION << ", got " << version;

        ans.sim_info = DoIOForSimulation(topNode.GetChildOrThrow("simulation"));

        ans.dataFilePath = DoIOForGeometry(topNode.GetChildOrThrow("geometry"));

        if (topNode.GetChildOrNull("colloids") != Element::Missing())
        {
            ans.colloid_xml_path = xmlFilePath;
        }

        ans.initial_condition = DoIOForInitialConditions(topNode.GetChildOrThrow("initialconditions"));

        ans.inlets = DoIOForInOutlets(ans.sim_info, topNode.GetChildOrThrow("inlets"));
        ans.outlets = DoIOForInOutlets(ans.sim_info, topNode.GetChildOrThrow("outlets"));

        // Optional element <properties>
        if (auto propertiesEl = topNode.GetChildOrNull("properties"))
            ans.propertyOutputs = DoIOForProperties(ans.sim_info, propertiesEl);

        // Optional element <monitoring>
        if (auto monitoringEl = topNode.GetChildOrNull("monitoring"))
            ans.monitoringConfig = DoIOForMonitoring(monitoringEl);

        // The RBC section must be parsed *after* the inlets and outlets have been
        // defined
        if (auto rbcEl = topNode.GetChildOrNull("redbloodcells")) {
            if constexpr (build_info::BUILD_RBC) {
                ans.rbcConf = DoIOForRedBloodCells(ans, rbcEl);
            } else {
                throw Exception() << "Input XML has redbloodcells section but HEMELB_BUILD_RBC=OFF";
            }
        }
        return ans;
    }

    GlobalSimInfo SimConfigReader::DoIOForSimulation(const Element simEl) const
    {
        GlobalSimInfo ans;
        // Required element
        // <stresstype value="enum lb::StressTypes" />
        ans.stress_type = [](unsigned v) {
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
        const Element stepsEl = simEl.GetChildOrThrow("steps");
        GetDimensionalValue(stepsEl, "lattice", ans.time.total_steps);

        // Required element
        // <step_length value="float" units="s" />
        const Element tsEl = simEl.GetChildOrThrow("step_length");
        GetDimensionalValue(tsEl, "s", ans.time.step_s);

        // Optional element
        // <extra_warmup_steps value="unsigned" units="lattice" />
        if (auto wuEl = simEl.GetChildOrNull("extra_warmup_steps"))
        {
            GetDimensionalValue(wuEl, "lattice", ans.time.warmup_steps);
            ans.time.total_steps += ans.time.warmup_steps;
        } else {
            ans.time.warmup_steps = 0;
        }

        // Required element
        // <voxel_size value="float" units="m" />
        const Element vsEl = simEl.GetChildOrThrow("voxel_size");
        GetDimensionalValue(vsEl, "m", ans.space.step_m);

        // Required element
        // <origin value="(x,y,z)" units="m" />
        const Element originEl = simEl.GetChildOrThrow("origin");
        GetDimensionalValue(originEl, "m", ans.space.geometry_origin_m);

        // Optional element
        // <fluid_density value="float" units="kg/m3" />
        ans.fluid.density_kgm3 = simEl.GetChildOrNull("fluid_density").transform(
                [](Element const& el) {
                    return GetDimensionalValue<PhysicalDensity>(el, "kg/m3");
                }).value_or(DEFAULT_FLUID_DENSITY_Kg_per_m3);

        // Optional element
        // <fluid_viscosity value="float" units="Pa.s" />
        ans.fluid.viscosity_Pas = simEl.GetChildOrNull("fluid_viscosity").transform(
                [](Element const& el) {
                    return GetDimensionalValue<PhysicalDynamicViscosity>(el, "Pa.s");
                }).value_or(DEFAULT_FLUID_VISCOSITY_Pas);

        // Optional element (default = 0)
        // <reference_pressure value="float" units="mmHg" />
        ans.fluid.reference_pressure_mmHg = simEl.GetChildOrNull("reference_pressure").transform(
                [](Element const& el) {
                    return GetDimensionalValue<PhysicalPressure>(el, "mmHg");
                }).value_or(0);

        ans.checkpoint = simEl.GetChildOrNull("checkpoint").transform(
                [](Element const& el) {
                    CheckpointInfo ans;
                    el.GetAttributeOrThrow("period", ans.period);
                    return ans;
                }
        );

        return ans;
    }

    auto SimConfigReader::DoIOForGeometry(const Element geometryEl) const -> path
    {
        // Required element
        // <geometry>
        //  <datafile path="relative path to GMY" />
        // </geometry>
        return RelPathToFullPath(geometryEl.GetChildOrThrow("datafile").GetAttributeOrThrow("path"));
    }

    /**
     * Helper function to ensure that the iolet being created matches the compiled
     * iolet BC.
     * @param ioletEl
     * @param requiredBC
     */
    void SimConfigReader::CheckIoletMatchesCMake(const Element& ioletEl,
                                                 std::string_view requiredBC) const
    {
        // Check that HEMELB_*LET_BOUNDARY is consistent with this
        auto ioletTypeName = ioletEl.GetName();
        std::string_view hemeIoletBC;

        if (ioletTypeName == "inlet")
            hemeIoletBC = build_info::INLET_BOUNDARY.view();
        else if (ioletTypeName == "outlet")
            hemeIoletBC = build_info::OUTLET_BOUNDARY.view();
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

    auto SimConfigReader::DoIOForInOutlets(GlobalSimInfo const& sim_info, const Element ioletsEl) const -> std::vector<IoletConfig>
    {
        auto nodeName = ioletsEl.GetName();
        auto childNodeName = std::string{nodeName.substr(0, nodeName.size() - 1)};
        std::vector<IoletConfig> ioletList;
        for (auto currentIoletNode: ioletsEl.Children(std::string{nodeName.substr(0, nodeName.size() - 1)})) {
            // Determine which InOutlet to create
            Element conditionEl = currentIoletNode.GetChildOrThrow("condition");
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
            std::visit([&](auto& conf) {
                           if constexpr (std::is_same_v<decltype(conf), std::monostate&>) {
                               throw Exception();
                           } else {
                               DoIOForBaseInOutlet(sim_info, currentIoletNode, conf);
                           }
                       },
                       newIolet);
            ioletList.push_back(std::move(newIolet));
        }
        return ioletList;
    }

    auto SimConfigReader::DoIOForPressureInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        CheckIoletMatchesCMake(ioletEl, "NASHZEROTHORDERPRESSUREIOLET");
        Element conditionEl = ioletEl.GetChildOrThrow("condition");
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

    auto SimConfigReader::DoIOForVelocityInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        CheckIoletMatchesCMake(ioletEl, "LADDIOLET");
        Element conditionEl = ioletEl.GetChildOrThrow("condition");
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



    std::vector<extraction::PropertyOutputFile> SimConfigReader::DoIOForProperties(GlobalSimInfo const& sim_info, const Element& propertiesEl) const
    {
        std::vector<extraction::PropertyOutputFile> propertyOutputs;
        ForExactChildren(propertiesEl, "propertyoutput",
                         [&](Element const& _) {
                             propertyOutputs.push_back(DoIOForPropertyOutputFile(sim_info, _));
                         });
        return propertyOutputs;
    }

    extraction::PropertyOutputFile SimConfigReader::DoIOForPropertyOutputFile(GlobalSimInfo const& sim_info,
                                                                              const Element& propertyoutputEl) const
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
            if (!io::TimePattern::Check(file.filename.native()))
                throw (Exception() << "For single timestep output files, "
                                      "the path must contain exactly one '%d' and no other '%' characters");
        }

        propertyoutputEl.GetAttributeOrThrow("period", file.frequency);

        Element geometryEl = propertyoutputEl.GetChildOrThrow("geometry");
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

        for (auto fieldEl: propertyoutputEl.Children("field"))
            file.fields.push_back(DoIOForPropertyField(sim_info, fieldEl));

        return file;
    }

    extraction::StraightLineGeometrySelector*
    SimConfigReader::DoIOForLineGeometry(const Element& geometryEl) const
    {
        std::vector<PhysicalPosition> points;
        for (auto child: geometryEl.Children()) {
            if (auto name = child.GetName(); name != "point")
                throw (Exception() << "Unknown child element: " << name);
            points.push_back(GetDimensionalValue<PhysicalPosition>(child, "m"));
        }

        if (auto n = points.size(); n != 2)
            throw (Exception() << "Must have exactly 2 points for line, instead have " << n);

        return new extraction::StraightLineGeometrySelector(points[0].as<float>(), points[1].as<float>());
    }



    extraction::PlaneGeometrySelector*
    SimConfigReader::DoIOForPlaneGeometry(const Element& geometryEl) const
    {
        auto [pointEl, normalEl] = GetExactChildren(geometryEl, "point", "normal");
        PhysicalPosition point;
        util::Vector3D<float> normal;

        GetDimensionalValue(pointEl, "m", point);
        GetDimensionalValue(normalEl, "dimensionless", normal);

        auto radius = geometryEl.GetChildOrNull("radius").transform(
                [](Element const& el) {
                    return GetDimensionalValue<PhysicalDistance>(el, "m");
                });
        if (radius)
            return new extraction::PlaneGeometrySelector(point.as<float>(), normal, *radius);
        else
            return new extraction::PlaneGeometrySelector(point.as<float>(), normal);
    }

    extraction::SurfacePointSelector*
    SimConfigReader::DoIOForSurfacePoint(const Element& geometryEl) const
    {
        Element pointEl = geometryEl.GetChildOrThrow("point");

        PhysicalPosition point;
        GetDimensionalValue(pointEl, "m", point);
        return new extraction::SurfacePointSelector(point.as<float>());
    }

    extraction::OutputField
    SimConfigReader::DoIOForPropertyField(GlobalSimInfo const& sim_info, const Element& fieldEl) const
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

    void SimConfigReader::DoIOForBaseInOutlet(GlobalSimInfo const& sim_info, const Element& ioletEl, IoletConfigBase& ioletConf) const
    {
        Element positionEl = ioletEl.GetChildOrThrow("position");
        Element normalEl = ioletEl.GetChildOrThrow("normal");

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
                    [](Element const &el) {
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
        for (auto insertEl: ioletEl.Children("insertcell")) {
            CellInserterConfig inserterConf;
            inserterConf.seed = insertEl.GetChildOrThrow("seed").GetAttributeOrThrow<PrngSeedType>("value");
            inserterConf.template_name = insertEl.GetAttributeOrThrow("template");

            // Rotate cell to align z axis with given position, and then z axis with flow
            // If phi == 0, then cell symmetry axis is aligned with the flow
            inserterConf.theta_rad = GetDimensionalValueWithDefault<Angle>(insertEl, "theta", "rad", 0e0);
            inserterConf.phi_rad = GetDimensionalValueWithDefault<Angle>(insertEl, "phi", "rad", 0e0);
            inserterConf.translation_m = {
                    GetDimensionalValueWithDefault<PhysicalDistance>(insertEl, "x", "m", 0e0),
                    GetDimensionalValueWithDefault<PhysicalDistance>(insertEl, "y", "m", 0e0),
                    GetDimensionalValueWithDefault<PhysicalDistance>(insertEl, "z", "m", 0e0)
            };

            inserterConf.offset = GetDimensionalValueWithDefault<
                    quantity_union<double, "s", "lattice">
            >(insertEl, "offset", quantity<double, "s">(0));
            inserterConf.drop_period_s = GetDimensionalValue<LatticeTime>(
                    insertEl.GetChildOrThrow("every"), "s");
            inserterConf.dt_s = GetDimensionalValueWithDefault<LatticeTime>(
                    insertEl, "delta_t", "s", 0e0);

            inserterConf.dtheta_rad = GetDimensionalValueWithDefault<Angle>(insertEl,
                                                                            "delta_theta",
                                                                            "rad",
                                                                            0e0);
            inserterConf.dphi_rad = GetDimensionalValueWithDefault<Angle>(insertEl,
                                                                          "delta_phi",
                                                                          "rad",
                                                                          0e0);
            inserterConf.dx_m = GetDimensionalValueWithDefault<LatticeDistance>(insertEl,
                                                                                "delta_x",
                                                                                "m",
                                                                                0e0);
            inserterConf.dy_m = GetDimensionalValueWithDefault<LatticeDistance>(insertEl,
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

    ICConfig SimConfigReader::DoIOForInitialConditions(Element initialconditionsEl) const
    {
        ICConfig initial_condition;
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
                Element uniformEl = pressureEl.GetChildOrThrow("uniform");
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
        return initial_condition;
    }

    auto SimConfigReader::DoIOForCosinePressureInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        CosinePressureIoletConfig newIolet;
        const Element conditionEl = ioletEl.GetChildOrThrow("condition");

        newIolet.amp_mmHg = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("amplitude"), "mmHg");
        newIolet.mean_mmHg = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("mean"), "mmHg");
        newIolet.phase_rad = GetDimensionalValue<Angle>(conditionEl.GetChildOrThrow("phase"), "rad");
        newIolet.period_s = GetDimensionalValue<LatticeTime>(conditionEl.GetChildOrThrow("period"), "s");
        return newIolet;
    }

    auto SimConfigReader::DoIOForFilePressureInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        FilePressureIoletConfig newIolet;
        const Element pathEl = ioletEl.GetChildOrThrow("condition").GetChildOrThrow("path");
        newIolet.file_path = RelPathToFullPath(pathEl.GetAttributeOrThrow("value"));

        return newIolet;
    }

    auto SimConfigReader::DoIOForMultiscalePressureInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        MultiscalePressureIoletConfig newIolet;
        const Element conditionEl = ioletEl.GetChildOrThrow("condition");

        const Element pressureEl = conditionEl.GetChildOrThrow("pressure");
        GetDimensionalValue(pressureEl, "mmHg", newIolet.pressure_reference_mmHg);

        const Element velocityEl = conditionEl.GetChildOrThrow("velocity");
        GetDimensionalValue(velocityEl, "m/s", newIolet.velocity_reference_ms);

        newIolet.label = conditionEl.GetChildOrThrow("label").GetAttributeOrThrow("value");
        return newIolet;
    }

    auto SimConfigReader::DoIOForParabolicVelocityInOutlet(
            const Element& ioletEl) const -> IoletConfig
    {
        ParabolicVelocityIoletConfig newIolet;

        const Element conditionEl = ioletEl.GetChildOrThrow("condition");

        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);
        GetDimensionalValue(conditionEl.GetChildOrThrow("maximum"), "m/s", newIolet.max_speed_ms);
        return newIolet;
    }

    auto SimConfigReader::DoIOForWomersleyVelocityInOutlet(const Element& ioletEl) const -> IoletConfig
    {
        WomersleyVelocityIoletConfig newIolet;

        const Element conditionEl = ioletEl.GetChildOrThrow("condition");

        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);
        GetDimensionalValue(conditionEl.GetChildOrThrow("pressure_gradient_amplitude"),
                            "mmHg/m", newIolet.pgrad_amp_mmHgm);
        GetDimensionalValue(conditionEl.GetChildOrThrow("period"), "s", newIolet.period_s);
        GetDimensionalValue(conditionEl.GetChildOrThrow("womersley_number"), "dimensionless", newIolet.womersley);
        return newIolet;
    }

    auto SimConfigReader::DoIOForFileVelocityInOutlet(
            const Element& ioletEl) const -> IoletConfig
    {
        FileVelocityIoletConfig newIolet;

        const Element conditionEl = ioletEl.GetChildOrThrow("condition");

        newIolet.file_path = RelPathToFullPath(conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value"));
        GetDimensionalValue(conditionEl.GetChildOrThrow("radius"), "m", newIolet.radius_m);

        return newIolet;
    }


    MonitoringConfig SimConfigReader::DoIOForMonitoring(const Element& monEl) const
    {
        MonitoringConfig monitoringConfig;
        Element convEl = monEl.GetChildOrNull("steady_flow_convergence");
        if (convEl != Element::Missing())
        {
            DoIOForSteadyFlowConvergence(convEl, monitoringConfig);
        }

        monitoringConfig.doIncompressibilityCheck = (monEl.GetChildOrNull("incompressibility")
                                                     != Element::Missing());
        return monitoringConfig;
    }

    void SimConfigReader::DoIOForSteadyFlowConvergence(const Element& convEl, MonitoringConfig& monitoringConfig) const
    {
        monitoringConfig.doConvergenceCheck = true;
        monitoringConfig.convergenceRelativeTolerance = convEl.GetAttributeOrThrow<double>("tolerance");
        monitoringConfig.convergenceTerminate = (convEl.GetAttributeOrThrow("terminate") == "true");

        auto n = 0;
        for (auto criterionEl: convEl.Children("criterion")) {
            DoIOForConvergenceCriterion(criterionEl, monitoringConfig);
            n++;
        }
        if (n == 0) {
            throw Exception() << "At least one convergence criterion must be provided in "
                              << convEl.GetPath();
        }
    }

    void SimConfigReader::DoIOForConvergenceCriterion(const Element& criterionEl, MonitoringConfig& monitoringConfig) const
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

    TemplateCellConfig SimConfigReader::readCell(const Element& cellNode) const {
        TemplateCellConfig ans;
        ans.name = cellNode.GetAttributeMaybe("name").value_or("default");
        auto const shape = cellNode.GetChildOrThrow("shape");
        ans.mesh_path = RelPathToFullPath(
                shape.GetAttributeOrThrow("mesh_path")
        );

        auto get_fmt = [](Element const& shape, char const* attr_name) -> MeshFormat {
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

    std::map<std::string, TemplateCellConfig> SimConfigReader::readTemplateCells(Element const& cellsEl) const {
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

    NodeForceConfig readNode2NodeForce(const Element& node) {
        auto intensity = GetDimensionalValueWithDefault<
                quantity_union<double, "Nm", "lattice">
        >(node, "intensity", quantity<double, "Nm">(1.0));

        auto cutoffdist = GetDimensionalValueWithDefault<LatticeDistance>(
                node, "cutoffdistance", "lattice", 1.0);

        // exponent doesnt have units apparently
        auto exponent = node.and_then(
                [](Element const& el) { return el.GetChildOrNull("exponent"); }
        ).transform(
                [](Element const& el) { return el.GetAttributeOrThrow<std::size_t>("value"); }
        ).value_or(std::size_t(2));

        return {intensity, cutoffdist, exponent};
    }

    RBCConfig SimConfigReader::DoIOForRedBloodCells(SimConfig const& conf, const Element &rbcEl) const {
        RBCConfig ans;

        const Element controllerNode = rbcEl.GetChildOrThrow("controller");
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
        if (!(std::all_of(conf.inlets.begin(), conf.inlets.end(), cell_inserter_templates_exist)
              &&
              std::all_of(conf.outlets.begin(), conf.outlets.end(), cell_inserter_templates_exist))) {
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
}