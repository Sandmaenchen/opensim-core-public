A CMake project template to test the AUKSMIKT code. Step-by-step instructions.

1. Copy the two files (CMakeLists.txt and runAUKSMIKT.cpp) in the directory of your choice (referred to as /root/ here).
2. In the CMakeLists.txt file, set the path where OpenSim core was installed in the find_package() command (line 15).
3. In the runAUKSMIKT.cpp file, set the paths and filenames.
4. At /root/, create subdirectory called 'build'
5. In terminal, navigate to /root/build/
6. run command ``cmake ..''
7. run command ``cmake --build .''
8. run  command ``./runAUKSMIKT''

If OpenSim had been built and the paths set successfully, you should see log messages from the execution of AUKSMIKT on your data.