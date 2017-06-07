# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/zesty64"
  config.vm.synced_folder "../KSP Assemblies", "/home/ubuntu/KSP Assemblies", id: "Assemblies"

  script = <<SCRIPT
echo Provisioning Principia
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
echo "deb http://download.mono-project.com/repo/debian wheezy main" | tee /etc/apt/sources.list.d/mono-xamarin.list
wget -O - http://apt.llvm.org/llvm-snapshot.gpg.key|sudo apt-key add -
add-apt-repository 'deb http://apt.llvm.org/wily/ llvm-toolchain-wily-4.0 main'
apt-get update
apt-get install -y clang-4.0 clang-format-4.0 clang-tidy-4.0 git unzip wget libc++-dev libc++abi-dev binutils make automake libtool curl cmake monodevelop
ln -s /usr/bin/clang-4.0 /usr/bin/clang
clang --version
ln -s /usr/bin/clang++-4.0 /usr/bin/clang++
clang++ --version
ln -s /usr/bin/clang-format-4.0 /usr/bin/clang-format
clang-format --version
ln -s /usr/bin/clang-format-diff-4.0 /usr/bin/clang-format-diff
clang-format-diff --help
ln -s /usr/bin/clang-tidy-4.0 /usr/bin/clang-tidy
clang-tidy --version
git clone https://github.com/mockingbirdnest/Principia.git principia
chown -R vagrant:vagrant principia KSP
cd principia
if [ -d deps ]; then echo "Dependencies already installed, remove deps/ to reinstall."; else ./install_deps.sh; fi
SCRIPT

  config.vm.provision "shell", inline: script

  config.vm.provider "virtualbox" do |v|
    v.memory = 2048
    v.cpus = 4
  end
end
