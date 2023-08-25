// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIGREADER_H
#define HEMELB_CONFIGURATION_SIMCONFIGREADER_H

#include <filesystem>


#include "quantity.h"
#include "configuration/SimConfig.h"
#include "extraction/GeometrySelectors.h"
#include "io/xml.h"

namespace hemelb::configuration {

    // Note on specifying physical quantities.
    //
    // In the XML, these are encoded in attributes of an element (the name
    // of the element doesn't matter for these purposes). The two required
    // attributes are:
    // 1. "units" (e.g. "m" for a physical length or "lattice" for something specified in scaled units)
    // 2. "value" which is a string encoding of the representation (typically a double)
    //
    // Now some parts of the code (e.g. RBC) can specify units as multiple options (e.g. physical OR lattice)
    // for these, use the quantity union types. The consumer (i.e. SimBuilder etc) will have
    // to take care around this.
    inline void check_unit_spec(const io::xml::Element& elem, std::string_view actual, std::string_view const& expected) {
        if (actual != expected)
            throw Exception() << "Invalid units for element " << elem.GetPath()
                              << "."" Expected '" << expected
                              << "', got '" << actual << "'";
    }

    //! Check the units of the quantity and decode the value into @param value
    template<typename T>
    void GetDimensionalValue(const io::xml::Element& elem, std::string_view units, T& value)
    {
        auto got = elem.GetAttributeOrThrow("units");
        check_unit_spec(elem, got, units);
        elem.GetAttributeOrThrow("value", value);
    }
    //! Check the units of the quantity and return the value
    template<typename T>
    T GetDimensionalValue(const io::xml::Element& elem, std::string_view units)
    {
        auto got = elem.GetAttributeOrThrow("units");
        check_unit_spec(elem, got, units);

        return elem.GetAttributeOrThrow<T>("value");
    }

    //! Given an element (@param elem), check for a child with the given @param name.
    //! If it exists, return the unit-checked (against @param unit) value.
    //! If it doesn't exist, return @param default_value.
    template <typename T>
    T GetDimensionalValueWithDefault(const io::xml::Element& elem,
                                     char const* name, char const* unit, T default_value) {
        using Element = io::xml::Element;
        return elem
                .and_then([&](const Element& _) { return _.GetChildOrNull(name); })
                .transform([&](const Element& _) { return GetDimensionalValue<T>(_, unit); })
                .value_or(default_value);
    }

    template <QuantityUnion QUnion, IsVariantAlternative<QUnion> DefaultQ>
    QUnion GetDimensionalValueWithDefault(const io::xml::Element& elem, char const* name, DefaultQ const& default_q) {
        using RepT = typename quantity_union_traits<QUnion>::representation_type;
        using io::xml::Element;
        return elem
                .and_then([&](Element const& _) { return _.GetChildOrNull(name); })
                .transform([](Element const& _) {
                    auto xml_units = _.GetAttributeOrThrow("units");
                    auto xml_val = _.GetAttributeOrThrow<RepT>("value");
                    return quantity_union_factory<QUnion>()(xml_val, xml_units);
                })
                .value_or(default_q);
    }

    class SimConfigReader {
        using path = std::filesystem::path;
        using Element = io::xml::Element;
        path xmlFilePath;

    public:
        SimConfigReader(path);

        virtual ~SimConfigReader() = default;

        [[nodiscard]] virtual SimConfig Read() const;

        // Turn an input XML-relative path into a full path
        [[nodiscard]] path RelPathToFullPath(std::string_view path) const;

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const Element &ioletEl,
                                            std::string_view requiredBC) const;

        [[nodiscard]] virtual SimConfig DoIO(const Element xmlNode) const;

        [[nodiscard]] virtual GlobalSimInfo DoIOForSimulation(const Element simEl) const;

        [[nodiscard]] virtual path DoIOForGeometry(const Element geometryEl) const;

        [[nodiscard]] virtual std::vector <IoletConfig>
        DoIOForInOutlets(GlobalSimInfo const &sim_info, const Element xmlNode) const;

        void
        DoIOForBaseInOutlet(GlobalSimInfo const &sim_info, const Element &ioletEl, IoletConfigBase &ioletConf) const;

        [[nodiscard]] IoletConfig DoIOForPressureInOutlet(const Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForCosinePressureInOutlet(const Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForFilePressureInOutlet(const Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForMultiscalePressureInOutlet(
                const Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForVelocityInOutlet(const Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForParabolicVelocityInOutlet(
                const Element &ioletEl) const;

        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForWomersleyVelocityInOutlet(const Element &ioletEl) const;

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForFileVelocityInOutlet(const Element &ioletEl) const;

        [[nodiscard]] virtual std::vector <extraction::PropertyOutputFile>
        DoIOForProperties(GlobalSimInfo const &sim_info, const Element &xmlNode) const;

        [[nodiscard]] virtual extraction::OutputField
        DoIOForPropertyField(GlobalSimInfo const &sim_info, const Element &xmlNode) const;

        [[nodiscard]] virtual extraction::PropertyOutputFile
        DoIOForPropertyOutputFile(GlobalSimInfo const &sim_info, const Element &propertyoutputEl) const;

        [[nodiscard]] extraction::StraightLineGeometrySelector *DoIOForLineGeometry(
                const Element &xmlNode) const;

        [[nodiscard]] extraction::PlaneGeometrySelector *DoIOForPlaneGeometry(const Element &) const;

        [[nodiscard]] extraction::SurfacePointSelector *DoIOForSurfacePoint(const Element &) const;

        [[nodiscard]] virtual ICConfig DoIOForInitialConditions(Element parent) const;
        //virtual void DoIOForCheckpointFile(const Element& checkpointEl) const;

        /**
         * Reads monitoring configuration from XML file
         *
         * @param monEl in memory representation of <monitoring> xml element
         */
        [[nodiscard]] virtual MonitoringConfig DoIOForMonitoring(const Element &monEl) const;

        /**
         * Reads configuration of steady state flow convergence check from XML file
         *
         * @param convEl in memory representation of the <steady_flow_convergence> XML element
         */
        void DoIOForSteadyFlowConvergence(const Element &convEl, MonitoringConfig &monitoringConfig) const;

        /**
         * Reads the configuration of one of the possible several converge criteria provided
         *
         * @param criterionEl in memory representation of the <criterion> XML element
         */
        void DoIOForConvergenceCriterion(const Element &criterionEl, MonitoringConfig &monitoringConfig) const;

        [[nodiscard]] virtual TemplateCellConfig readCell(const Element &cellNode) const;

        [[nodiscard]] virtual std::map <std::string, TemplateCellConfig>
        readTemplateCells(Element const &cellsEl) const;

        [[nodiscard]] virtual RBCConfig DoIOForRedBloodCells(SimConfig const &, const Element &rbcEl) const;

    };
}
#endif