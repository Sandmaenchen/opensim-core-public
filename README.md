<!-- OpenSim Logo -->
<p align=center>
    <a href="https://opensim.stanford.edu/">
        <img src="https://drive.google.com/uc?id=1urYfucgR4pCM5OeXySMBVc3i5oGfYRAf" alt="OpenSim Logo">
</p>

<!-- Badges -->
<p align=center>
    <a href="https://github.com/opensim-org/opensim-core/actions">
        <img src="https://github.com/opensim-org/opensim-core/workflows/continuous-integration/badge.svg" alt="Continuous Integration Badge">
    </a>
    <a href="https://github.com/opensim-org/opensim-core/releases">
        <img src="https://img.shields.io/github/v/release/opensim-org/opensim-core?include_prereleases" alt="Releases Badge">
    </a>
    <a href="https://github.com/opensim-org/opensim-core/blob/master/LICENSE.txt">
        <img src="https://img.shields.io/hexpm/l/apa" alt="License Badge">
    </a>
    <a href="https://github.com/opensim-org/opensim-core/wiki/Build-Instructions">
        <img src="https://img.shields.io/badge/platform-windows%20%7C%20macos%20%7C%20linux-lightgrey" alt="Supported Platforms Badge">
    </a>
    <a href="https://github.com/opensim-org/opensim-core/graphs/contributors">
        <img src="https://img.shields.io/github/contributors/opensim-org/opensim-core" alt="ZenHub Badge">
    </a>
    <a href="https://zenhub.com">
        <img src="https://img.shields.io/badge/Shipping%20faster%20with-ZenHub-blueviolet" alt="ZenHub Badge">
    </a>
</p>

---

**NOTE: This branch contains the source code for my customized version of OpenSim 4.5.**


## Building and installing on Ubuntu

Instructions partially follow [this reference](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53085346/Scripting+in+Python) [Accessed: 2024-05-13].

1. In terminal, run the command `sudo apt update`, followed by {`sudo apt upgrade`}.
2. Install Eigen 3.4.0 with the command `sudo apt install libeigen3-dev`.
3. Copy the Eigen directory with the command `sudo cp -r /usr/include/eigen3/Eigen/ /usr/local/include/Eigen`.
4. Download the script `opensim-core-ukf-linux-build-script.sh` from https://github.com/Sandmaenchen/opensim-core-public/tree/ukf-uks-tools/scripts/build. Ensure that *CORE_BRANCH* flag is set to *ukf-uks-tools*. 
5. Make the script runnable with `chmod +x opensim-core-ukf-linux-build-script.sh`
6. Run the script with the command `./opensim-core-ukf-linux-build-script.sh`. 
* This can take a lot of time.
7. Install Python setup tools with the command `sudo apt-get install python-setuptools`.
8. Navigate: `cd ~/opensim-core/sdk/Python`.
9. Run the command `sudo python3 setup.py install`.
* By default, this command installs to `/usr/` directory where your user has no writing permission by default. If, in addition, `root` cannot access the directory where OpenSim was built and installed (e.g., your user's `/home/` directory), you can set the *PYTHONUSERBASE* environment variable (add this also to your bash profile file) to point to the directory accessible by your user; afterwards, install the package with the command `python3 setup.py install --user`.
10. Edit your bash profile file with the commmand `nano ~/.bashrc`. Add the line `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: /home/<your_username>/opensim-core/sdk/Simbody/lib}` at the end of the bash profile file.
11. Run the command `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: /home/<your_username>/opensim-core/sdk/Simbody/lib` in order to be able to use OpenSim in the current session.

## How to run on Windows 11

We haven't been able to build OpenSim on Windows 11. If you wish to test our software on Windows machine, the easiest path would be to run a Linux virtual machine (VM) using Python API.

### Setting up and running WSL2

1. Select Start $\rightarrow$ Turn Windows features on or off. Ensure the following are enabled:
* Hyper-V
* Virtual Machine Platform
* Windows PowerShell 2.0
* Windows Subsystem for Linux
2. Start PowerShell (as admin), and run the command `wsl --install`. Set your username and password when asked. Restart.
3. Set WSL2 as default by running in PowerShell the command `wsl --set-default-version 2`.
4. Update with the command `wsl --update`.
5. Change to correct Linux distribution with the command `wsl --install -d Ubuntu-22.04`. 
6. Run the VM and login with the command `wsl --user <your_username>` (Quit with the command `exit`).

### Accessing files on virtual machine disk

Open File Explorer and navigate to `\\wsl\$`. The mounting point of VM (Ubuntu-22.04) is the root directory (i.e., `/`) of the Linux system.

## An example template for Python script

The following is a template for Python script that creates an object from our UKF-based tool class and runs it.

~~~~
import opensim as osim\
tool = osim.UKFIMUInverseKinematicsTool()
tool.set_model_file(modelFileName)
tool.set_orientations_file(orientationsFileName)
tool.set_sensor_to_opensim_rotations(osim.Vec3(0, 0, 0))
tool.set_results_directory(resultsDirectory)
tool.set_alpha(1.0)
tool.set_beta(2.0)
tool.set_kappa(-1.337) #value for (3-n)
tool.set_order(2)
tool.set_lag_length(5)
tool.set_num_threads(6)
tool.set_write_UKF(True)
tool.set_sgma2w(2.0**14)
tool.set_missing_data_scale(100.0)
tool.set_imu_RMS_in_deg(osim.Vec3(0.5, 1.0, 0.5))
tool.set_enable_resampling(True)
tool.set_process_covariance_method(0)
tool.set_time_range(0, startTime) 
tool.set_time_range(1, endTime)   
tool.run(False, osim.Vector((order+1), 1.0))
~~~~



## License [![License](https://img.shields.io/hexpm/l/apa)](https://github.com/opensim-org/opensim-core/blob/master/LICENSE.txt)

Licensed under the Apache License, Version 2.0.  See the full text of the [Apache License, Version 2.0](https://github.com/opensim-org/opensim-core/blob/master/LICENSE.txt) for more information. This license makes OpenSim suitable for commercial, government, academic, and personal use. 

Third-party components have their own licenses. see our [Notice](https://github.com/opensim-org/opensim-core/blob/master/NOTICE), and [Acknowledgements](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Acknowledgements) webpages for more information.

### How to acknowledge us

Acknowledging the OpenSim project helps us and helps you. It allows us to track our impact, which is essential for securing funding to improve the software and provide support to our users (you). If you use OpenSim, we would be extremely grateful if you acknowledge us by citing the following paper:

> Seth A, Hicks JL, Uchida TK, Habib A, Dembia CL, et al. (2018) **OpenSim: Simulating musculoskeletal dynamics and neuromuscular control to study human and animal movement.** _PLOS Computational Biology_ 14(7): e1006223. https://doi.org/10.1371/journal.pcbi.1006223

If you use plugins, models, or other components contributed by your fellow researchers, you must acknowledge their work as described in the license that accompanies each of these files.


## Funding

The OpenSim project is currently supported by the following:
 - United States National Institutes of Health (NIH)
    - [Mobilize Center](https://mobilize.stanford.edu/) (P41 EB027060)
    - [Restore Center](https://restore.stanford.edu/) (P2C HD101913)
 - [Wu Tsai Human Performance Alliance](https://humanperformancealliance.org/)

Past funding includes the following grants and contracts:

 - United States National Institutes of Health (NIH)
    - Simulation of Biological Structures (Simbios; U54 GM072970)
    - Simulation in Rehabilitation Research (NCSRR; R24 HD065690, P2C HD065690)
    - Mobilize Center (U54 EB020405)
 - United States Defense Advanced Research Projects Agency (DARPA)
    - Warrior Web (W911QX-12-C-0018)
