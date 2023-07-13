## Astra:

#### MKL and iFort installation:
    
``` bash
Ubuntu Installation (tested on 20.04 in AMD Ryzen 9)

wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-compiler-fortran
sudo apt install intel-oneapi-mkl-devel
source /opt/intel/oneapi/setvars.sh --force
```
    

#### Astra installation:
    TODO: Made some changes to some of the Makefiles
        1. Switched deprecated -mkl usage to -qmkl
        2. Changed MKL_PATH and MKL_INCLUDE to reflect the OneAPI installation

  

