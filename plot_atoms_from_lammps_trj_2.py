
import matplotlib.pyplot as plt
import numpy as np


def simple_plot():
    plot_all = "no"
    stopping_step = 1444
    
    atoms_to_plot = [1668]
    steps_start = 0 
    time_steps = []
    x_cor = []
    y_cor= []
    z_cor = []


    
    trajectory_file = "30_variable_pull_test.lammpstrj"
    path_and_file = "C:\\LAMMPS\\tensile_DNW\\30_nonH\\" + trajectory_file
    
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


simple_plot()
#single_plot_defined_steps()
#multiple_plot()
#compare_plots()