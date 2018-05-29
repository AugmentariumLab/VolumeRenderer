# VolumeRenderer
This code compresses a volumetric dataset (specifically the Richtmeyer-Instability dataset) using a progressive kd-tree scheme, saves the data, decompresses, and renders it. 

#### Quick Start Tutorial:
1. Download and uncompress the timesteps of the Richtmeyer-Instability dataset you wish to visualize. There will be one binary file per brick per timestep.
2. Clone the VolumeRenderer repo
3. Change the 'findBinaryBrickFile' function on line #580 in main.cpp to return the Richtmeyer-Instability brick on your computer from the timestep and the brick number.
4. Run the project in Visual Studio.

##### To change the renderer:
You can switch from a color compositing renderer to a isosurface renderer by changing the shaders on lines 71-74 of main.cpp.

##### To change the bricks and timesteps to process:
Change the parameters in the LoadBricksToTexture function call on line 242 of main.cpp. 
> LoadBricksToTexture(number_of_bricks, size_x, size_y, size_z, timestep, send_to_gpu)

##### To set the number of epochs and error tolerance for compression:
Adjust parameters on lines 253 & 254 in main.cpp.

##### To save the compressed data:
Change the filename in line 267 in main.cpp.

