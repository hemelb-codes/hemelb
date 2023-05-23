// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_H
#define HEMELB_LB_STREAMERS_H

#include "build_info.h"

#include "lb/streamers/BulkStreamer.h"
#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/SimpleBounceBack.h"
#include "lb/streamers/BouzidiFirdaousLallemand.h"
#include "lb/streamers/GuoZhengShi.h"
#include "lb/streamers/JunkYang.h"
#include "lb/streamers/NashZerothOrderPressure.h"
#include "lb/streamers/LaddIolet.h"

namespace hemelb::lb {
    namespace detail {
        using namespace streamers;
        template <typename C>
        constexpr auto get_default_wall_streamer(InitParams& i) {
            constexpr auto WALL = build_info::WALL_BOUNDARY;
            if constexpr (WALL == "BFL") {
                return StreamerTypeFactory < BouzidiFirdaousLallemandLink < C >, NullLink < C >> {i};
            } else if constexpr (WALL == "GZS") {
                return StreamerTypeFactory < GuoZhengShiLink < C >, NullLink < C >> {i};
            } else if constexpr (WALL == "SIMPLEBOUNCEBACK") {
                return StreamerTypeFactory < BounceBackLink < C >, NullLink < C >> {i};
            } else if constexpr (WALL == "JUNKYANG") {
                return JunkYangFactory<NullLink<C> >{i};
            } else {
                throw (Exception() << "Configured with invalid WALL_BOUNDARY");
            }
        }

        template <ct_string NAME, typename C>
        constexpr auto get_default_iolet_streamer(InitParams& i) {
            if constexpr (NAME == "NASHZEROTHORDERPRESSUREIOLET") {
                return StreamerTypeFactory<
                        NullLink<C>,
                        NashZerothOrderPressureLink < C >
                >{i};
            } else if constexpr (NAME == "LADDIOLET") {
                return StreamerTypeFactory<
                        NullLink<C>,
                        LaddIoletLink < C >
                >{i};
            } else {
                throw (Exception() << "Configured with invalid IOLET boundary");
            }
        }

    }

    template <typename C>
    using DefaultStreamer = streamers::BulkStreamer<C>;

    template <typename C>
    using DefaultWallStreamer = decltype(detail::get_default_wall_streamer<C>(std::declval<InitParams&>()));

    template <typename C>
    using DefaultInletStreamer = decltype(detail::get_default_iolet_streamer<build_info::INLET_BOUNDARY, C>(std::declval<InitParams&>()));

    template <typename C>
    using DefaultOutletStreamer = decltype(detail::get_default_iolet_streamer<build_info::OUTLET_BOUNDARY, C>(std::declval<InitParams&>()));

    // Given wall and iolet streamers construct a combination in the `type` output member
    // Primary template
    template <typename WS, typename IS>
    struct CombineWallAndIoletStreamers;

    // Specialisation for "standard" wall streamers:
    // need to deduce the underlying wall and iolet
    // delegates out of the type factory.
    template<
            typename C, // collision
            template <typename> class WallT, // wall delegate template
            template <typename> class IoletT // iolet delegate template
    >
    struct CombineWallAndIoletStreamers<
            streamers::StreamerTypeFactory<WallT<C>, streamers::NullLink<C>>,
            streamers::StreamerTypeFactory<streamers::NullLink<C>, IoletT<C>>
    > {
        using type = streamers::StreamerTypeFactory<WallT<C>, IoletT<C>>;
    };

    // Junk Yang is different: the pure wall streamer has a special tag type for no-iolet
    template<
            typename C, // collision
            template <typename> class IoletT // iolet delegate template
    >
    struct CombineWallAndIoletStreamers<
            streamers::JunkYangFactory<streamers::NullLink<C> >,
            streamers::StreamerTypeFactory<streamers::NullLink<C>, IoletT<C>>
    > {
        using type = streamers::JunkYangFactory<IoletT<C>>;
    };
}
#endif /* HEMELB_LB_STREAMERS_H */
