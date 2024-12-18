# SpheroPAK3D
A Plug-in for Generating Periodic Composite Cells with Random Sphere Inclusions in Abaqus/CAE

Developer: Youngbin LIM <br>
Contact: lyb0684@naver.com

![Github_Fig](https://github.com/user-attachments/assets/bb6609c8-fb16-44e2-9507-a87ced5df065)

# Supported versions
Abaqus 2022 and higher

# Installation
Download the contents in SpheroPAK3D and put those files in one folder. Move the folder to the abaqus_plugins directory (ex: C:\Users\User_Name\abaqus_plugins). Restart the Abaqus/CAE and the plug-in will be available whenever the Abaqus/CAE is running. The plug-in is activated for all modules and can be found in Plug-ins drop down menu.

# Geometry and Mesh
Spheres with random sizes are generated based on user-defined parameter (avaerage radius and standard deviation of radius) <br>
Periodic images of spheres are generated along with random sphere generation, meaning that the geometry is completely periodic <br> However, free meshing algorithm in Abaqus is used and the mesh is not periodic. <br>
Thus, surface elements with regular mesh is generated and tied with the cell surface to enforce the periodic boundary condition <br>

![Github_Fig_2](https://github.com/user-attachments/assets/280ba2e4-bf56-49a5-bc3c-5ab5d0d27ed3)


# Periodic boundary condition
Rectangular shell part is generated and then meshed with surface elements (SFM3D4). The surface elements are tied with RVE cell surface to enforce the periodic boudnary conditions. Face, edge and corner node sets are first generated and equation constraints are applied. <br> Hide the surface element part in Abaqus/Viewer when analyzing the results

![Github_Fig_3](https://github.com/user-attachments/assets/b6bd69f2-9500-46f9-a3ea-978e2561a142)


# Disclaimer
The author do not represent DASSAULT SYSTEMES SIMULIA CORP (Collectively "DS") and the plug-in is not part of official product provided by DS. Your use of this download shall be at your sole risk. No support of any kind of the download is provided by DS. Regarding the inquiries on the plug-in, please contact developer.
