// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_H
#define HEMELB_LB_STREAMERS_H

#include "build_info.h"

#include "lb/streamers/SimpleCollideAndStream.h"
#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiDelegate.h"
#include "lb/streamers/JunkYangFactory.h"
#include "lb/streamers/NashZerothOrderPressureDelegate.h"
#include "lb/streamers/LaddIoletDelegate.h"

namespace hemelb::lb {
    namespace detail {
        using namespace streamers;
        template <typename C>
        constexpr auto get_default_wall_streamer(InitParams& i) {
            constexpr auto WALL = build_info::WALL_BOUNDARY;
            if constexpr (WALL == "BFL") {
                return WallStreamerTypeFactory<C, BouzidiFirdaousLallemandDelegate<C>>{i};
            } else if constexpr (WALL == "GZS") {
                return WallStreamerTypeFactory<C, GuoZhengShiDelegate<C>>{i};
            } else if constexpr (WALL == "SIMPLEBOUNCEBACK") {
                return WallStreamerTypeFactory<C, SimpleBounceBackDelegate<C>>{i};
            } else if constexpr (WALL == "JUNKYANG") {
                return JunkYangFactory<C, NoIoletLink<C> >{i};
            } else {
                throw (Exception() << "Configured with invalid WALL_BOUNDARY");
            }
        }

        template <ct_string NAME, typename C>
        constexpr auto get_default_iolet_streamer(InitParams& i) {
            if constexpr (NAME == "NASHZEROTHORDERPRESSUREIOLET") {
                return IoletStreamerTypeFactory<
                        C,
                        NashZerothOrderPressureDelegate<C>
                >{i};
            } else if constexpr (NAME == "LADDIOLET") {
                return IoletStreamerTypeFactory<
                        C,
                        LaddIoletDelegate<C>
                >{i};
            } else {
                throw (Exception() << "Configured with invalid IOLET boundary");
            }
        }

    }

    template <typename C>
    using DefaultStreamer = streamers::SimpleCollideAndStream<C>;

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
            streamers::WallStreamerTypeFactory<C, WallT<C>>,
            streamers::IoletStreamerTypeFactory<C, IoletT<C>>
    > {
        using type = streamers::WallIoletStreamerTypeFactory<C, WallT<C>, IoletT<C>>;
    };

    // Junk Yang is different: the pure wall streamer has a special tag type for no-iolet
    template<
            typename C, // collision
            template <typename> class IoletT // iolet delegate template
    >
    struct CombineWallAndIoletStreamers<
            streamers::JunkYangFactory<C, streamers::NoIoletLink<C> >,
            streamers::IoletStreamerTypeFactory<C, IoletT<C>>
    > {
        using type = streamers::JunkYangFactory<C, IoletT<C>>;
    };
}
#endif /* HEMELB_LB_STREAMERS_H */
