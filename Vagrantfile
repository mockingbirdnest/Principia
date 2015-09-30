# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
  config.vm.box = "debian/jessie64"
  config.vm.synced_folder ".", "/home/vagrant/principia"

  script = <<SCRIPT
echo Provisioning Principia
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
apt-get update
apt-get install -y clang git unzip wget libc++-dev binutils make automake libtool curl cmake monodevelop
cd principia
if [ -d deps ]; then echo "Dependencies already installed, remove deps/ to reinstall." else ./install_deps.sh; fi
SCRIPT

  config.vm.provision "shell", inline: script
end
