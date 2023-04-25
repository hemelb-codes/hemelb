// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIWINDOW_H
#define HEMELB_NET_MPIWINDOW_H

#include <span>
#include "net/MpiCommunicator.h"

namespace hemelb::net
{
    // An MPI Window created with MPI_Win_allocate that gives access to the data, PGAS style
    template <typename T>
    class WinData {
        using value_type = T;
        using span = std::span<T>;
        using traits = net::MpiDataTypeTraits<value_type>;

        // We own the window (and data) so invoke rule of five
        MPI_Win window = MPI_WIN_NULL;
        span data;

    public:
        WinData() = default;

        // Construct window and allocate data - collective and `n`
        // must be the same across the communicator.
        WinData(MPI_Aint n, MpiCommunicator comm, value_type const& init) {
            MPI_Info info;
            MPI_Info_create(&info);
            MPI_Info_set(info, "same_size", "true");
            MPI_Info_set(info, "same_disp_unit", "true");
            T* tmp;
            HEMELB_MPI_CALL(
                    MPI_Win_allocate,
                    (n * sizeof(value_type), sizeof(value_type), info, comm, &tmp, &window)
            );
            MPI_Info_free(&info);
            data = span(tmp, n);
            // Windows do not inherit the handler from parent communicator
            HEMELB_MPI_CALL(MPI_Win_set_errhandler, (window, MPI_ERRORS_RETURN));
            std::fill(data.begin(), data.end(), init);
        }

        // No copying, uniquely owns the window
        WinData(WinData const&) = delete;
        WinData& operator=(WinData const&) = delete;

        // Moving is fine
        WinData(WinData&& other) : WinData() {
            // Default constructor sets us to null, so swap leaves other in null state :)
            std::swap(window, other.window);
            std::swap(data, other.data);
        }
        // Collective if we own a real window.
        WinData& operator=(WinData&& other) {
            Free();
            // Free sets us to null
            std::swap(window, other.window);
            std::swap(data, other.data);
            return *this;
        }

        // Destructor frees the window and data - collective
        ~WinData() noexcept(false) {
            Free();
        }

        void Free() {
            if (window != MPI_WIN_NULL) {
                HEMELB_MPI_CALL(MPI_Win_free, (&window));
                data = span{};
            }
        }

        // Get a span to our data
        span Span() const {
            return data;
        }

        // End an RMA epoch and ensure all communications within it are complete (collective)
        void Fence(int assertions = 0) const {
            HEMELB_MPI_CALL(MPI_Win_fence, (assertions, window));
        }

        // Possibly remote reference to a data item
        struct const_reference {
            int trank;
            MPI_Aint tdisp;
            WinData const * win;

            operator T() const {
                // A lock/unlock with MPI_MODE_NOCHECK should have zero overhead
                // on proper RDMA machines but does satisfy MPI RMA requirements.
                HEMELB_MPI_CALL(MPI_Win_lock,
                                (MPI_LOCK_SHARED, trank, MPI_MODE_NOCHECK, win->window)
                );
                int ans;
                HEMELB_MPI_CALL(MPI_Get, (
                        &ans, 1, traits::GetMpiDataType(),
                                trank, tdisp, 1, traits::GetMpiDataType(),
                                win->window
                ));
                HEMELB_MPI_CALL(MPI_Win_unlock,
                                (trank, win->window)
                );
                return ans;
            }
        };

        struct reference : public const_reference {
            // Up to the caller to deal with synchronisation!
            reference operator=(value_type const& val) {
                HEMELB_MPI_CALL(MPI_Put, (
                        &val, 1, traits::GetMpiDataType(),
                                this->trank, this->tdisp, 1, traits::GetMpiDataType(),
                                this->win->window
                ));
                return *this;
            }
        };

        // This type deals with synchronisation for a set of writes
        //
        // Requires that there NOT be an open RMA epoch on construction.
        // You can then do zero or more writes (MPI_Put - a pure one-sided operation),
        // then destructor calls fence (collective) and these are then guaranteed to have occurred.
        //
        // Does NOT start a new epoch so it puts you back where you started.
        class WriteSession {
            WinData* win;
        public:
            // Collective as it starts an epoch with fence.
            WriteSession(WinData* w) : win{w} {
                // Start the writing epoch
                win->Fence(MPI_MODE_NOPRECEDE);
            }
            // No copying: really we should ensure there is only one
            // active per window but leave that to the programmer.
            WriteSession(WriteSession const&) = delete;
            WriteSession& operator=(WriteSession const&) = delete;
            // Moving could be possible, but disallow to reduce the
            // chances of having multiple sessions going
            WriteSession(WriteSession &&) = delete;
            WriteSession& operator=(WriteSession&&) = delete;

            ~WriteSession() noexcept(false) {
                win->Fence(MPI_MODE_NOSUCCEED);
            }
            reference operator()(int rank, MPI_Aint i) {
                return {rank, i, win};
            }
        };

        // Start a load of writes
        WriteSession begin_writes() {
            return {this};
        }

        // Get reference to a possibly remote element
        // Read-only
        const_reference operator()(int rank, MPI_Aint i) const {
            return {rank, i, this};
        }
        // Read-write
        reference operator()(int rank, MPI_Aint i) {
            return {rank, i, this};
        }

    };
}

#endif
