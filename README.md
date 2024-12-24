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
Spheres with random sizes are generated based on user-defined parameters (average radius and standard deviation of the radius). Periodic images of spheres are created simultaneously with the random sphere generation, ensuring that the geometry is completely periodic. However, the free meshing algorithm in Abaqus is used, and as a result, the mesh is not periodic. <br>

![Github_Fig_2](https://github.com/user-attachments/assets/280ba2e4-bf56-49a5-bc3c-5ab5d0d27ed3)


# Periodic boundary condition
A rectangular shell part is created and meshed with surface elements (SFM3D4). These surface elements are tied to the RVE cell surface to enforce periodic boundary conditions. Face, edge, and corner node sets are first generated, and equation constraints are applied. The surface elements have no stiffness and do not affect the overall stiffness of the system. When analyzing the results in Abaqus/Viewer, hide the surface element part to avoid visual clutter.

![Shell](https://github.com/user-attachments/assets/5c654b2f-ac04-45b0-b2ab-1c46f8b7756a)

# Disclaimer
The author do not represent DASSAULT SYSTEMES SIMULIA CORP (Collectively "DS") and the plug-in is not part of official product provided by DS. Your use of this download shall be at your sole risk. No support of any kind of the download is provided by DS. Regarding the inquiries on the plug-in, please contact developer.
