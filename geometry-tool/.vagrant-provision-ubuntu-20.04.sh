# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# Ubuntu pacakages
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y \
    cmake \
    gcc-9 \
    jq \
    libboost-dev \
    libcgal-dev \
    libopengl0 \
    libosmesa6 \
    libxt6 \
    ninja-build \
    python-is-python3 \
    python3.8-venv \
    python3.8-dev \
    python3-pip \
    python3-wxgtk4.0 \
    xfce4 \
    virtualbox-guest-dkms \
    virtualbox-guest-utils \
    virtualbox-guest-x11

# Allow non-root to run GUI
echo "allowed_users=anybody" > /etc/X11/Xwrapper.config
# NOTE: VM requires reboot before running GUI
echo "Before using the GUI run \`vagrant reload gmy\`"
