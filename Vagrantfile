# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  config.vm.define "gmy" do |gmy|
    gmy.vm.box = "ubuntu/focal64"

    # Provider-specific configuration
    gmy.vm.provider "virtualbox" do |vb|
      # Display the VirtualBox GUI when booting the machine
      vb.gui = true
      # Make sure the graphics aren't tiny
      vb.customize ["modifyvm", :id, "--graphicscontroller", "vmsvga"]
      # Moar RAM and GPU
      vb.memory = "4096"
      vb.cpus = 4
    end

    # Install things
    gmy.vm.provision "shell", path: "geometry-tool/.vagrant-provision-ubuntu-20.04.sh"
  end
end
