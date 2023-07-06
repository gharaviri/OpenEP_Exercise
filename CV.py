import sys
import os
import math
import pyvista as pv
import numpy as np

class CV:
    
    def Triangulation(self,egmX,LAT,data_type=''):
        ''' This function calculates local atrial conduction velocity from recorded 
                electro-anatomical data, using triangulation method. 
                (C.D. Cantwell et.al. 2015,''Techniques for automated local activation 
                time annotation and conduction velocity estimation in cardiac mapping")
            
            egmX: Recorded Electrode locations (x,y,z). electric --> bipolar_egm --> points
            LAT : Recorded local activation time. electric --> annotations --> local_activation_time - reference Time
            '''
        ##### initiation ######
        cv = []
        min_theta = 30
        max_elec_distance = 10
        min_elec_distance = 1.5
        min_lat_difference = 2
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        cv_centers_id = []
        alpha = 5
        
        ######## Create a triangulation mesh from the recording electrodes
        egmX_point_cloud = pv.PolyData(egmX)
        print('---> Triangulation process')
        surf = egmX_point_cloud.delaunay_3d(alpha=alpha)
        print('---> Triangulation process')
        print('---> CV calculation (Triangulation method started ...)')
        for num_point in range(len(surf.cells_dict[5])):
            
            vtx_id = []
            lat = []
            
            for i in range (3):
                vtx_id.append(surf.cells_dict[5][num_point][i])
                lat.append(LAT[surf.cells_dict[5][num_point][i]])
            
            id_lat_sorted  = np.argsort(lat)

            O = [egmX[int(vtx_id[id_lat_sorted[0]])],lat[id_lat_sorted[0]]]
            A = [egmX[int(vtx_id[id_lat_sorted[1]])],lat[id_lat_sorted[1]]]
            B = [egmX[int(vtx_id[id_lat_sorted[2]])],lat[id_lat_sorted[2]]]
            
            OA = np.sqrt(sum(np.power(np.subtract(O[0],A[0]),2)))
            OB = np.sqrt(sum(np.power(np.subtract(O[0],B[0]),2)))
            AB = np.sqrt(sum(np.power(np.subtract(A[0],B[0]),2)))
            
            tOA = A[1] - O[1]
            tOB = B[1] - O[1]
            
            theta = np.arccos((np.power(OA,2)+ np.power(OB,2) - np.power(AB,2))/(2 * OA * OB))
            # check if the conditions are meet to accept the triangle set as a valid one
            if (math.degrees(theta) >= min_theta and OA >= min_elec_distance and OA <= max_elec_distance and
                OB >= min_elec_distance and OB <= max_elec_distance and tOA >= min_lat_difference and tOB >= min_lat_difference):

                alpha = np.arctan((tOB * OA - tOA * OB * np.cos(theta)) / (tOA * OB * np.sin(theta)))
                cv_temp = (OA/tOA) * np.cos(alpha)
                
                if cv_temp >= 0.2 and cv_temp <=2:
                    cv.append(cv_temp)
                    cv_centers_x.append(O[0][0])
                    cv_centers_y.append(O[0][1])
                    cv_centers_z.append(O[0][2])
                    cv_centers_id.append(vtx_id[id_lat_sorted[0]])
                else:
                    
                    cv_centers_x.append(O[0][0])
                    cv_centers_y.append(O[0][1])
                    cv_centers_z.append(O[0][2])
                    cv_centers_id.append(vtx_id[id_lat_sorted[0]])
                    cv.append(cv_temp)
            ####### Create a triangulation mesh from egmX pointcloud #####
        cv_centers = [cv_centers_x, cv_centers_y, cv_centers_z]
           
        print('---> CV calculation ended')
        return cv, cv_centers, cv_centers_id
    
    
    def visualization(self, data = '',mesh = '', egmX = ''):
        ''' This function prepare and visualize calculated data
            data: Calculated data (e.g. conduction velocity or voltage)
            egmX: Recorded Electrode locations (x,y,z). electric --> bipolar_egm --> points
            mesh: Surface mesh provided by electro-anatomical mapping device. you can access it 
                  using this command "surface = case.create_mesh()" openep-py library
        '''
        
        boring_cmap = plt.cm.get_cmap("jet", 256)
        
        pl = pv.Plotter()
        
        egmX_point_cloud = pv.PolyData(egmX)
        egmX_point_cloud['scalars'] = data
        interpolated_mesh_cv = mesh.interpolate(egmX_point_cloud,radius=10,strategy='closest_point')
        pl.add_mesh(egmX_point_cloud,scalars='scalars',render_points_as_spheres=True, 
                    point_size=20, cmap = boring_cmap)
        pl.add_mesh(interpolated_mesh_cv,scalars='scalars', clim = [0, 1], cmap = boring_cmap)
        pl.show()
            
            