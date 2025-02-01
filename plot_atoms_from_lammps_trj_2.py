
import matplotlib.pyplot as plt
import math
import numpy as np


def plot_energy_vs_distance():
    reference_atoms = ["1", "2"] #this are atom types for now
    time_steps = []
    distance_data = []    
    trajectory_file = "AIREBO_C2H6_variable_pull_test.lammpstrj"
    path = "C:\\LAMMPS\\topo\\C2H6\\AIREBO RUN\\"
    path_and_file = path + trajectory_file
    log_file = path  + "log.lammps"
    count1 = 0
    count2 = 0
    atom1= np.zeros(3)
    atom2= np.zeros(3)   
    
    with open(path_and_file) as f:
        igot = f.readlines()         
        
        for line in igot:
                atom_data = line.split() 
                if len(atom_data) >2 :
                    
                        
                    #print(line)
                    countzero_in1 = np.count_nonzero(atom1)
                    countzero_in2 = np.count_nonzero(atom2)

                    
                    if atom_data[1] == "1":
                        #print("found reference atom", line)  
                        #print(atom_data[1])
                        atom1_x, atom1_y, atom1_z = float(atom_data[2]), float(atom_data[3]) , float(atom_data[4])
                        atom1 = np.array([atom1_x, atom1_y, atom1_z])
                        count1 =count1 +1 
                        #print(atom1_x, atom1_y, atom1_z)
                    
                    if atom_data[1] == "2":
                        #print("found reference atom", line)  
                        #print(atom_data[1])
                        atom2_x, atom2_y, atom2_z = float(atom_data[2]), float(atom_data[3]) ,float(atom_data[4])
                        atom2 = np.array([atom2_x, atom2_y, atom2_z])
                        count2 =count2 +1
                        #print(atom2_x, atom2_y, atom2_z)
                    
                   
                    if countzero_in1 != 0 and countzero_in2 != 0 :
                        #print(atom1, atom2)
                        #distance = np.linalg.norm(atom2 - atom1) ##############################################try the manual method first to double check                    
                        distance = math.sqrt((atom2_x - atom1_x)**2 +  (atom2_y -atom1_y)**2 +  (atom2_z- atom1_z)**2)
                        
                        if distance < 2:
                            #print("this is distance", distance)
                            distance_data.append(distance)
                            
                        #testing is all the steps are measured only once:
                        #distance_data.append(distance)
                        atom1.fill(0)
                        atom2.fill(0)
                        #print(atom1, atom2)
                                                  

        try:                        
            
    
            print("steps to break", len(distance_data), ": bond distance at breaking step", distance_data[-1], "Angstroms")
            print(count1, count2)  
        except:
            pass
    
    #extract the forces at each step from the log file
    forces = []    
    with open(log_file) as f:
        igot = f.readlines()
        for count, line in enumerate(igot):
            if line.find("999") > -1 and "dataflag" not in line:  #this needs to skip over the dataflag line in the input sections
                force_line= line.split()
                forces.append(force_line[1])
                #print(line)        
        print("force data points", len(forces))
    
    #calculate the energy by multiplying the force * the distance up to the breaking step
    energies = []
    eq_length = 1.54
    for point in range(len(distance_data)):
        
        if point == 0:
            delta_d = round(float(distance_data[point]) -  eq_length , 4 )
            
        else:
            delta_d = round(float(distance_data[point]), 4) - round(float(distance_data[point-1]), 4)
        
        energy = delta_d  * float(forces[point])
        #print(delta_d, forces[point], energy)
        energies.append(energy)
    
    #sum up all the enerygy to supplement the graph
    total_e =  0
    for delta_e in range(len(energies)):
        total_e = round(float(energies[delta_e]) + total_e, 2)
    
    BDE = "Bond dissociation energy: " +  str(total_e) + " eV"    
    print("Bond dissociation energy: ", total_e)
    bond_length = "Bond length at breaking step: " + str( round(float(distance_data[-1]),2)) +  " Angstroms"
    
    last_data_point = len(distance_data)
                   
    plt.figure()                
    plt.plot([float(i) for i in distance_data[0:last_data_point]], [float(j) for j in energies[0:last_data_point]], color='r', label="energy")
    plt.text(1.57, 0.1, BDE, fontsize=10, color='black', zorder=20)
    plt.text(1.57, 0.08, bond_length, fontsize=10, color='black', zorder=20)
    plt.legend(loc='upper left')
    plt.xlabel("distance A")
    plt.ylabel("energy eV")
    
    plt.title(" AIREBO bond dissociation energy")
    
    plt.show()        
        




def simple_plot():
    plot_all = "no"
    stopping_step = 1444
    
    atoms_to_plot = [1668]
    steps_start = 0 
    time_steps = []
    x_cor = []
    y_cor= []
    z_cor = []


    
    trajectory_file = "AIREBO_C2H6_variable_pull_test.lammpstrj"
    path_and_file = "C:\\LAMMPS\\topo\\C2H6\\AIREBO RUN\\" + trajectory_file
    
    with open(path_and_file) as f:
        igot = f.readlines()  
        
        #lets get how many steps are in the file by seeing how many times the atoms groups are written out
        
        for count, line in enumerate(igot):
            
            if line.find("ITEM: TIMESTEP") > -1:
                
                #print(count) #this is a timestep line
                new_step = (int(igot[count+1]))
                time_steps.append(new_step)
                
              
        #somewhere around here, if plot all = yes then stopping step equals last value in time_steps list else it stopping step
        if plot_all == "yes":
            last_data_point =  len(time_steps)
        else:
            last_data_point = stopping_step 
        
        print("there are this many data points in total", len(time_steps)) #this is how many data points we'll have for the atoms in atoms_to_plot
        
       
        for atoms in range(len(atoms_to_plot)):
            for lines in igot:
                atom_data = lines.split()
            #if line.find("ITEM: ATOMS id type xs ys zs") > -1:
                #print("found set of coordinates for a time step")
               
                if atom_data[0] == str(atoms_to_plot[atoms]):
                    #print(atom_data[0])
                    x_cor2append = atom_data[2]
                    x_cor.append(x_cor2append)
                    
                    y_cor2append = atom_data[3]
                    y_cor.append(y_cor2append)         
                    
                    z_cor2append = atom_data[4]
                    z_cor.append(z_cor2append)
                    
                    
        plt.plot([float(i) for i in time_steps[0:last_data_point]], [float(j) for j in x_cor[0:last_data_point]], color='r', label="X")
        #plt.plot([float(i) for i in time_steps], [float(j) for j in y_cor], color='g', label="Y")
        #plt.plot([float(i) for i in time_steps], [float(k) for k in z_cor], color='b', label="Z")
        
        plt.legend(loc='upper left')
        plt.xlabel("step number")
        plt.ylabel("position")
        plt.title(trajectory_file)
        
        plt.show()        









#########################################plots and atom from a single trajecotry file##############################
def single_plot_defined_steps():  #this always plots the first and last step, so it is not used to zoom in on a section of the graph

    atoms_to_plot = [5]
    how_many_timesteps = 1200  
    
    steps_start = 0 #leave this as a variable for graphing data starting in the middle
    
    time_steps = []
    time_steps_cut = []
    x_cor = []
    y_cor = []
    z_cor = []
    x_cor_cut = []
    y_cor_cut = []
    z_cor_cut = []
    
    trajectory_file = "ethane_variable_pull_test_with_FF4_run2break.lammpstrj"
    path_and_file = "C:\\LAMMPS\\topo\\C2H6\\" + trajectory_file
    
    with open(path_and_file) as f:
        igot = f.readlines()
    
        #lets get how many steps are in the file by seeing how many times the atoms groups are written out
        numberSteps = 0
        for count, line in enumerate(igot):
            
            if line.find("ITEM: TIMESTEP") > -1:
                numberSteps = numberSteps + 1
                #print(count) #this is a timestep line
                new_step = (int(igot[count+1]))
                time_steps.append(new_step)
            
        for atoms in range(len(atoms_to_plot)):
            for lines in igot:
                atom_data = lines.split()
            #if line.find("ITEM: ATOMS id type xs ys zs") > -1:
                #print("found set of coordinates for a time step")
               
                if atom_data[0] == str(atoms_to_plot[atoms]):
                    #print(atom_data[0])
                    x_cor2append = atom_data[2]
                    x_cor.append(x_cor2append)
                    
                    y_cor2append = atom_data[3]
                    y_cor.append(y_cor2append)         
                    
                    z_cor2append = atom_data[4]
                    z_cor.append(z_cor2append)                
                    #print(line)
                
        #print(numberSteps, time_steps)
        #print("X", x_cor)
        #print("Y", y_cor)
        #print("Z", z_cor)
        
    ##get number of time steps
    num_steps = len(time_steps)
    #print(num_steps)
    
    #extracting a subset of the full data set based on how many point I actually want to have on the graph. 
    ##always append the starting step and last step
    time_steps_cut.append(time_steps[steps_start])
    time_steps_cut.append(time_steps[num_steps-1])
    
    x_cor_cut.append(x_cor[steps_start]) 
    x_cor_cut.append(x_cor[num_steps-1])
    
    y_cor_cut.append(y_cor[steps_start]) 
    y_cor_cut.append(y_cor[num_steps-1])
    
    z_cor_cut.append(z_cor[steps_start]) 
    z_cor_cut.append(z_cor[num_steps-1])
    
    ##get step increment
    step_increment = int(num_steps / how_many_timesteps)
    data_index = step_increment
    
    for steps in range(1, how_many_timesteps):
        time_steps_cut.insert(steps, time_steps[data_index]) 
        x_cor_cut.insert(steps, x_cor[data_index])
        y_cor_cut.insert(steps, y_cor[data_index])
        z_cor_cut.insert(steps, z_cor[data_index])
        ##print(steps)
        #print(x_cor_cut, time_steps_cut)
        data_index = data_index + step_increment
        #print(data_index)    
    
    plt.plot([float(i) for i in time_steps_cut], [float(j) for j in x_cor_cut], color='r', label="X")
    #plt.plot([float(i) for i in time_steps_cut], [float(j) for j in y_cor_cut], color='g', label="Y")
    #plt.plot([float(i) for i in time_steps_cut], [float(k) for k in z_cor_cut], color='b', label="Z")
    
    plt.legend(loc='upper left')
    plt.xlabel("step number")
    plt.ylabel("position")
    plt.title(trajectory_file)
    
    plt.show()
#####################################################################################

def multiple_plot():
    import glob
    files_list = (glob.glob("C:/LAMMPS/topo/hinge files//Data/*.lammpstrj")) 
    
    for trajs in files_list:
        #these values get reset for each traj value, but hopefully the graph can be appended to the graph and i persistent
        atoms_to_plot = [1]
        how_many_timesteps = 2
        
        steps_start = 0 #leave this as a variable for graphing data starting in the middle
        
        time_steps = []
        time_steps_cut = []
        x_cor = []
        y_cor = []
        z_cor = []
        x_cor_cut = []
        y_cor_cut = []
        z_cor_cut = [] 
        
        with open(trajs) as f:
            print(trajs)
            igot = f.readlines()
            numberSteps = 0
            for count, line in enumerate(igot):
                
                if line.find("ITEM: TIMESTEP") > -1:
                    numberSteps = numberSteps + 1
                    #print(count) #this is a timestep line
                    new_step = (int(igot[count+1]))
                    time_steps.append(new_step)
                
            for atoms in range(len(atoms_to_plot)):
                for lines in igot:
                    atom_data = lines.split()
                #if line.find("ITEM: ATOMS id type xs ys zs") > -1:
                    #print("found set of coordinates for a time step")
                   
                    if atom_data[0] == str(atoms_to_plot[atoms]):
                        #print(atom_data[0])
                        x_cor2append = atom_data[2]
                        x_cor.append(x_cor2append)
                        
                        y_cor2append = atom_data[3]
                        y_cor.append(y_cor2append)         
                        
                        z_cor2append = atom_data[4]
                        z_cor.append(z_cor2append)                        
                        
            ##get number of time steps
            num_steps = len(time_steps)
            print(num_steps)
            
            #extracting a subset of the full data set based on how many point I actually want to have on the graph. 
            ##always append the starting step and last step
            time_steps_cut.append(time_steps[steps_start])
            time_steps_cut.append(time_steps[num_steps-1])
            
            x_cor_cut.append(x_cor[steps_start]) 
            x_cor_cut.append(x_cor[num_steps-1])
            
            y_cor_cut.append(y_cor[steps_start]) 
            y_cor_cut.append(y_cor[num_steps-1])
            
            z_cor_cut.append(z_cor[steps_start]) 
            z_cor_cut.append(z_cor[num_steps-1])
            
            ##get step increment
            step_increment = int(num_steps / how_many_timesteps)
            data_index = step_increment
            
            for steps in range(1, how_many_timesteps):
                time_steps_cut.insert(steps, time_steps[data_index]) 
                x_cor_cut.insert(steps, x_cor[data_index])
                y_cor_cut.insert(steps, y_cor[data_index])
                z_cor_cut.insert(steps, z_cor[data_index])
                ##print(steps)
                #print(x_cor_cut, time_steps_cut)
                data_index = data_index + step_increment
                #print(data_index)       
                
            plt.plot([float(i) for i in time_steps_cut], [float(j) for j in x_cor_cut], color='r', label="X")
            plt.plot([float(i) for i in time_steps_cut], [float(j) for j in y_cor_cut], color='g', label="Y")
            plt.plot([float(i) for i in time_steps_cut], [float(k) for k in z_cor_cut], color='b', label="Z")
            
            plt.legend(loc='upper left')
            plt.xlabel("step number")
            plt.ylabel("position")
            plt.title(str(trajs))
            
            plt.show()        



#################################################################compare trajectories in directory####################################
def compare_plots():    
    import glob
    
    files_list = (glob.glob("C:/LAMMPS/DATA/*.lammpstrj")) 
    
    for dimensions in range(0,3):
    
        for count, trajs in enumerate(files_list):
            if count == 0: 
                line_color = 'r'
            elif count == 1:
                line_color = 'g'
            else:
                line_color = 'b'
            #these values get reset for each traj value, but hopefully the graph can be appended to the graph and i persistent
            atoms_to_plot = [1]
            how_many_timesteps = 2
            
            steps_start = 0 #leave this as a variable for graphing data starting in the middle
            
            time_steps = []
            time_steps_cut = []
            x_cor = []
            y_cor = []
            z_cor = []
            dim_cor_cut = []
    
            
            with open(trajs) as f:
                print(trajs)
                label_name = str(trajs).split('\\') 
                
                
                igot = f.readlines()
                numberSteps = 0
                for count, line in enumerate(igot):
                    
                    if line.find("ITEM: TIMESTEP") > -1:
                        numberSteps = numberSteps + 1
                        #print(count) #this is a timestep line
                        new_step = (int(igot[count+1]))
                        time_steps.append(new_step)
                    
                for atoms in range(len(atoms_to_plot)):
                    for lines in igot:
                        atom_data = lines.split()
                    #if line.find("ITEM: ATOMS id type xs ys zs") > -1:
                        #print("found set of coordinates for a time step")
                       
                        if atom_data[0] == str(atoms_to_plot[atoms]):
                            #print(atom_data[0])
                            x_cor2append = atom_data[2]
                            x_cor.append(x_cor2append)
                            
                            y_cor2append = atom_data[3]
                            y_cor.append(y_cor2append)         
                            
                            z_cor2append = atom_data[4]
                            z_cor.append(z_cor2append)                        
                            
                ##get number of time steps
                num_steps = len(time_steps)
                print(num_steps)
                
                #extracting a subset of the full data set based on how many point I actually want to have on the graph. 
                ##always append the starting step and last step
                time_steps_cut.append(time_steps[steps_start])
                time_steps_cut.append(time_steps[num_steps-1])
                
                
                ###get step increment
                step_increment = int(num_steps / how_many_timesteps)
                data_index = step_increment            
                
                
                #for the three xyz dimensions
                
                if dimensions == 0:
                    plt.title("X")
                    dim_cor_cut.append(x_cor[steps_start]) 
                    dim_cor_cut.append(x_cor[num_steps-1])
                    
                    for steps in range(1, how_many_timesteps):
                        time_steps_cut.insert(steps, time_steps[data_index]) 
                        dim_cor_cut.insert(steps, x_cor[data_index])    
                        data_index = data_index + step_increment
                    
                elif dimensions == 1:
                    plt.title("Y")
                    dim_cor_cut.append(y_cor[steps_start]) 
                    dim_cor_cut.append(y_cor[num_steps-1])
                    
                    for steps in range(1, how_many_timesteps):
                        time_steps_cut.insert(steps, time_steps[data_index]) 
                        dim_cor_cut.insert(steps, y_cor[data_index])    
                        data_index = data_index + step_increment                
                        
                else:
                    plt.title("Z")
                    dim_cor_cut.append(z_cor[steps_start]) 
                    dim_cor_cut.append(z_cor[num_steps-1])  
                    
                    for steps in range(1, how_many_timesteps):
                        time_steps_cut.insert(steps, time_steps[data_index]) 
                        dim_cor_cut.insert(steps, z_cor[data_index])    
                        data_index = data_index + step_increment                
                    
                      
                plt.plot([float(i) for i in time_steps_cut], [float(j) for j in dim_cor_cut], color=line_color, label= label_name[1])
    
                
        plt.legend(loc='upper left')
        plt.xlabel("step number")
        plt.ylabel("position")
        
        picture_name = "DIMENSION_"+str(dimensions)
        plt.savefig("C:/LAMMPS/DATA/"+picture_name+".png")
 
        plt.show()
        plt.clf()
def make_scatter_plot():
    plot_all = "yes"
    stopping_step = 1444
    
    atoms_to_plot = [19174]
    steps_start = 0 
    time_steps = []
    x_cor = []
    y_cor= []
    z_cor = []
    
    
    
    trajectory_file = "12_0_revolute_joint.lammpstrj"
    path_and_file = "C:\\LAMMPS\\rovelute joint plot compare\\" + trajectory_file
    
    with open(path_and_file) as f:
        igot = f.readlines()  
        
        #lets get how many steps are in the file by seeing how many times the atoms groups are written out
        
        for count, line in enumerate(igot):
            
            if line.find("ITEM: TIMESTEP") > -1:
                
                #print(count) #this is a timestep line
                new_step = (int(igot[count+1]))
                time_steps.append(new_step)
                
              
        #somewhere around here, if plot all = yes then stopping step equals last value in time_steps list else it stopping step
        if plot_all == "yes":
            last_data_point =  len(time_steps)
        else:
            last_data_point = stopping_step 
        
        print("there are this many data points in total", len(time_steps)) #this is how many data points we'll have for the atoms in atoms_to_plot
        
       
        for atoms in range(len(atoms_to_plot)):
            for lines in igot:
                atom_data = lines.split()
            #if line.find("ITEM: ATOMS id type xs ys zs") > -1:
                #print("found set of coordinates for a time step")
               
                if atom_data[0] == str(atoms_to_plot[atoms]):
                    #print(atom_data[0])
                    x_cor2append = atom_data[2]
                    x_cor.append(float(x_cor2append))
                    
                    y_cor2append = atom_data[3]
                    y_cor.append(float(y_cor2append))       
                    
                    z_cor2append = atom_data[4]
                    z_cor.append(float(z_cor2append))
                    
        x_var = statistics.variance(x_cor)
        y_var = statistics.variance(y_cor)
        z_var = statistics.variance(z_cor)
        
        print(round(x_var, 5) , round(y_var,5), round(z_var,5))
        
        xmax_value = max(x_cor)
        xmin_value = min(x_cor)
        ymax_value = max(y_cor)
        ymin_value = min(y_cor) 
        zmax_value = max(z_cor)
        zmin_value = min(z_cor)      
                    
        # Create a figure and 3D axes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Create the scatter plot
        ax.scatter(x_cor, y_cor, z_cor, c='r', marker='o')
        ax.set_xticks([xmin_value, xmax_value ])
        ax.set_yticks([ymin_value, ymax_value])
        ax.set_zticks([zmin_value, zmax_value])
        
        # Set labels
        ax.set_xlabel('X ')
        ax.set_ylabel('Y ')
        ax.set_zlabel('Z ')
        
    
        
        # Update the axis view and title
        ax.view_init(60, 60, 0)
        #plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))
    
        plt.draw()
        plt.show()
        plt.clf()



#plot_energy_vs_distance()
#simple_plot()
#single_plot_defined_steps()
#multiple_plot()
#compare_plots()
make_scatter_plot()
