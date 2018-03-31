# -*- mode: ruby -*-
# vi: set ft=ruby :

# Note: Use 'sudo su' in the ubuntu machine before make'ing.
Vagrant.configure("2") do |config|

  config.vm.define "ubuntu" do |ubuntu|
    ubuntu.vm.box = "ubuntu/bionic64"
    ubuntu.vm.box_url = "https://app.vagrantup.com/ubuntu/boxes/bionic64/versions/20180315.0.0/providers/virtualbox.box"
    ubuntu.vm.synced_folder "../KSP Assemblies", "/home/ubuntu/KSP Assemblies", id: "Assemblies"

    script = <<SCRIPT
echo Provisioning Principia
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
wget -O - http://apt.llvm.org/llvm-snapshot.gpg.key|sudo apt-key add -
add-apt-repository 'deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-5.0 main'
apt-get update
apt-get install -y clang-5.0 clang-tools-5.0 clang-format-5.0 clang-tidy-5.0 git unzip wget libc++-dev libc++abi-dev binutils make automake libtool curl cmake monodevelop
ln -s /usr/bin/clang-5.0 /usr/bin/clang
clang --version
ln -s /usr/bin/clang++-5.0 /usr/bin/clang++
clang++ --version
ln -s /usr/bin/clang-format-5.0 /usr/bin/clang-format
clang-format --version
ln -s /usr/bin/clang-format-diff-5.0 /usr/bin/clang-format-diff
clang-format-diff --help
ln -s /usr/bin/clang-tidy-5.0 /usr/bin/clang-tidy
clang-tidy --version
git clone https://github.com/mockingbirdnest/Principia.git principia
chown -R vagrant:vagrant principia KSP
cd principia
if [ -d deps ]; then echo "Dependencies already installed, remove deps/ to reinstall."; else ./install_deps.sh; fi
SCRIPT

    ubuntu.vm.provision "shell", inline: script

    ubuntu.vm.provider "virtualbox" do |v|
      v.memory = 2048
      v.cpus = 4
    end
  end
end
