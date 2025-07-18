# Spatium
Plug-in for Computational Homogenization of Particle-Reinforced Composite with Spherical Inclusions in Abaqus/CAE

Developer: Youngbin LIM <br>
Contact: lyb0684@naver.com

![image](https://github.com/user-attachments/assets/435cd454-bc40-49c0-a8b5-3d880ff02343)


Supported versions
--------------------------
Abaqus 2024 and higher <br>

Installation
--------------------------
Download the contents in Spatium and put those files in one folder. Move the folder to the abaqus_plugins directory (ex: C:\Users\User_Name\abaqus_plugins). Restart the Abaqus/CAE and the plug-in will be available whenever the Abaqus/CAE is running. The plug-in is activated for all modules and can be found in Plug-ins drop down menu.

Geometry and Mesh
--------------------------
Spheres with random sizes are generated based on user-defined parameters (average radius and standard deviation of the radius). Periodic images of spheres are created simultaneously with the random sphere generation, ensuring that the geometry is completely periodic. However, the free meshing algorithm in Abaqus is used, and as a result, the mesh is not periodic. <br>

![Github_Fig_2](https://github.com/user-attachments/assets/280ba2e4-bf56-49a5-bc3c-5ab5d0d27ed3)


Periodic boundary condition
--------------------------
A rectangular shell part is created and meshed with surface elements (SFM3D4). These surface elements are tied to the RVE cell surface to enforce periodic boundary conditions. Face, edge, and corner node sets are first generated, and equation constraints are applied. The surface elements have no stiffness and do not affect the overall stiffness of the system. When analyzing the results in Abaqus/Viewer, hide the surface element part.

![Fig](https://github.com/user-attachments/assets/92e38222-745b-4eed-966c-726d71a88ea4)

Volume average
--------------------------
The homogenized stress-strain curve is obtained by using the following equation.

![image](https://github.com/user-attachments/assets/bd732954-9da2-414b-90e5-88b9b1a449ed)


Disclaimer
--------------------------
The author do not represent DASSAULT SYSTEMES SIMULIA CORP (Collectively "DS") and the plug-in is not part of official product provided by DS. Your use of this download shall be at your sole risk. No support of any kind of the download is provided by DS. Regarding the inquiries on the plug-in, please contact developer.

Citation:
--------------------------
Y. Lim, M. Song, and S. Ha, "Spatium: An open-source Abaqus/CAE plug-in for computational homogenization of particle-reinforced composites" SoftwareX, vol. 31, p. 102219, 2025. DOI: 10.1016/j.softx.2025.102219 <br><br>

LaTex
```bibtex
@article{lim2025spatium,
  title={Spatium: An open-source Abaqus/CAE plug-in for computational homogenization of particle-reinforced composites},
  author={Lim, Youngbin; Song Minseo and Ha, Sangyul},
  journal={SoftwareX},
  volume={31},
  pages={102219},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.softx.2025.102219},
  url={https://doi.org/10.1016/j.softx.2025.102219}
}
