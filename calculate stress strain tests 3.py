import re
import matplotlib.pyplot as plt


project_path = "C:\\LAMMPS\\tensile_DNW\\30_nonH\\"  
mmp_file = project_path + "30_pull_group_tensile_DNW_nonH.mmp" #just use the trajectory file, include stress force checks
trajectory_file = project_path +  "30_variable_pull_test.lammpstrj"
log_file = project_path  + "30_pull_group_log.lammps"

#get a list of the atoms in the anchor and pull groups of the mmp


#get a list of the atoms in the anchor and pull groups
atomsets = {}
 

with open(mmp_file) as f:
    igot = f.readlines()

    for line_raw in igot:
        if line_raw.find("(Clipboard)") > -1:
            break    
        if line_raw.find("atomset") > -1:
                #This section creates a dictionary for all the atomsets. It is possible to just use atomsets and not jigs. 
                jig_name_raw = re.findall(r'\(.*?\)', line_raw)
                atom_set_name = jig_name_raw[0].replace('(','').replace(')','')
                #atomset_name_list.append(atom_set_name)
                print("name list", atom_set_name)
                clean_line_1 = line_raw.replace('(0, 0, 0)', '')
                clean_line = clean_line_1.strip('\n')
                
                s1 = [int(i) for i in clean_line.split() if i.isdigit()]
                atomsets[atom_set_name] = s1
    
        line = line_raw.replace('(','').replace(')','')              
        line = line.replace(',','')

#print(atomsets)

#get number of pull group atoms


num_pull_group_atoms = 0

for outer_key in atomsets:
    #print(outer_key)
    if outer_key == "Pull_X":
        for inner_key in atomsets[outer_key]: 
            num_pull_group_atoms = num_pull_group_atoms + 1
            
print("number of atoms in pull group", num_pull_group_atoms)


#get number of steps
#this doesn't seem necessary as I search again down below
#step_num = 0
#with open(trajectory_file) as dumpfile:
        #igot = dumpfile.readlines()
        #for line in igot:
            #if line.find("TIMESTEP") > -1: 
                #step_num = step_num + 1
        #print("total number of steps in dump file", step_num)
    

#####calculate strain

Fx_sum = 0
Fy_sum = 0
Fz_sum = 0
X_centroid = []
anchor_centroid = -33.74 #from mmptolammps_atomsets
X_strain_list = []


numberSteps = 0

with open(trajectory_file) as dumpfile:
        igot = dumpfile.readlines()

    #lets get how many steps are in the file by seeing how many times the atoms groups are written out
    
        #extract the number of atoms
        for count, line in enumerate(igot):
                if line.find("ITEM: NUMBER OF ATOMS") > -1:
                        NumAtoms = int(igot[count+1])
                        print('number of atoms ', NumAtoms)
                        break
        
        for count, line in enumerate(igot):
                
                                
                if line.find("ITEM: ATOMS id type x y z fx fy fz") > -1 :
                        numberSteps = numberSteps + 1
                        groupnum = 0
                        Lx = 0
                        Ly = 0
                        Lz = 0
                        #print("atoms start on line ", count+1)
                        start_lines = count+1
                        #print(igot[start_lines:start_lines+NumAtoms])
            
                        for atoms in range(start_lines,start_lines+NumAtoms):
                                atom_info =  igot[atoms].split()
                                
                                atom_id = atom_info[0]
                                atom_type = atom_info[1]
                                atomLoc_x = atom_info[2]
                                atomLoc_y = atom_info[3]
                                atomLoc_z = atom_info[4]
                                atom_Fx = atom_info[5]
                                atom_Fy = atom_info[6]
                                atom_Fz = atom_info[7]
                                atom_name = "atom_"+atom_id+"_"+"type_"+atom_type
                                
                                #pull group
                                if atom_type == "3" or atom_type == "4":
                                    groupnum = groupnum +  1
                                    Fx_sum = Fx_sum + float(atom_Fx)
                                    Lx =  Lx + float(atomLoc_x)
                                    Ly =  Ly + float(atomLoc_y)
                                    Lz =  Lz + float(atomLoc_z)
                                
                            
                        Cx = round(Lx/groupnum,4)
                        Cy = round(Ly/groupnum, 4)
                        Cz = round(Lz/groupnum, 4)
                        X_centroid.append(Cx)
                        #print(Cx, groupnum, numberSteps)
                                    
                                
  
#print("Fx SUM", Fx_sum)
original_len  = abs(X_centroid[0] - anchor_centroid)

for steps in range(len(X_centroid)):
    l = abs(X_centroid[steps] - anchor_centroid)
    
    X_strain =  abs((l-original_len)/original_len)
    #print("X strain", X_strain)
    X_strain_list.append(X_strain)


print(X_strain_list)    
# tries to find the bearking step by setting dif to a reasonable value
m = 0
n= 1 
print("there are these number of points", len(X_strain_list), " in the strain list")

plot_all = 0

if plot_all == 1 :
    breaking_step = len(X_strain_list)-1
else:

    for index, value in enumerate(X_strain_list):
        A = float(X_strain_list[m])
        
        B = float(X_strain_list[n])
        dif = B - A
        if dif > .2:
            print("possible break", index)
            breaking_step =  index
            break
        else:
            #print(dif)
            m = m + 1 
            n = n + 1  
            
#plot strain up to the breaking step
X_strain_plot = []
  
for points in range(breaking_step):
    X_strain_plot.append(X_strain_list[points])


plt.plot([float(i) for i in X_strain_plot])
#plt.legend(loc='upper left')
plt.xlabel("step number")
plt.ylabel("X_strain_raw")

title = "X Strain for " + str(num_pull_group_atoms) + " Pull Group Atoms " + "breaking at step "+ str(len(X_strain_plot))

plt.title(title) #add pull num to title strain then save file
plt.savefig(project_path + title + '.png')
plt.show()    
    

all_force_data = []
all_stress_data = []
X_force_plot_data = []
X_stress_plot_data = []


#set units flag
metal = 0
real = 1
if metal == 1:
    force_convert = 0.6242 #eV/A to 1 nN
    print("USING METAL UNITS FOR AIREBO")
else:
    force_convert = 0.0694 #kcal/mol/A to 1nN
    print("USING REAL UNITS FOR FF")

with open(log_file) as logfile:
        igot_force = logfile.readlines() 
        for line in igot_force:
            if line.find(" 999 ") > -1 and "dataflag" not in line:
                    force_data = line.split()
                    force_data_nN = float(force_data[1])/force_convert * num_pull_group_atoms #1nN = 0.6242 ev/Ang Metal units in LAMMPS, and LAMMPS documentation says this is per atom so * the number of atoms in the pull group....
                    cross_sec_area = 1.5 #######################SET THIS FOR THE STRUCTURE BEING TESTED in nm^2: 1nm * 1.5nm for the DNW
                    stress = force_data_nN / cross_sec_area
                    all_force_data.append(force_data_nN)
                    all_stress_data.append(stress)
            
        
        for force_points in range(breaking_step):            
            X_force_plot_data.append(all_force_data[force_points])
            X_stress_plot_data.append(all_stress_data[force_points])
            
        #Plot stress
        length_max_stress = len(X_stress_plot_data)
        max_stress = round(float(X_stress_plot_data[length_max_stress-1]),1)
        plt.plot([float(i) for i in X_stress_plot_data])
        #plt.legend(loc='upper left')
        plt.xlabel("step number")
        plt.ylabel("X Stress nano-Pa")
        stress_plot_name = "X Stress for "+ str(num_pull_group_atoms)+ " Pull Group Atoms"
        plt.title(stress_plot_name)
        plt.figtext(.3, 0.8, "Maxium stress = "+ str(max_stress)+" nano-Pa")
        plt.savefig(project_path + stress_plot_name + '.png')
        plt.show() 
        
        #plot LOG force
        length_max_force = len(X_force_plot_data)
        max_force = round(float(X_force_plot_data[length_max_force-1]),1)            
        plt.plot([float(i) for i in X_force_plot_data])
        #plt.legend(loc='upper left')
        plt.figtext(.5, 0.8, "Maxium force = "+ str(max_force)+" nN")
        plt.xlabel("step number")
        plt.ylabel("X Pull Force (nN) ")
        force_plot_name = ("X Pull Force "+ str(num_pull_group_atoms)+ " Pull Group Atoms")
        plt.title(force_plot_name)
        plt.savefig(project_path + force_plot_name + '.png')
        plt.show()             
    
    
#print(len(X_strain_plot))
#print(len(X_stress_plot_data))

#plot stress strain

plt.plot([float(i) for i in X_strain_plot], [float(j) for j in X_stress_plot_data], color='r')

plt.xlabel("Strain")
plt.ylabel("Stress (nN/nm^2)")
title = "Stress-Strain " + str(num_pull_group_atoms)+ " Pull Group Atoms"
plt.title(title)
plt.savefig(project_path + title + '.png')

plt.show()



breaking_step() 
create_stress_strain_plot()








#graph_bond_lenghts()
#def graph_bond_lenghts():
    ##this is for the last step, could do this for an inputed step or mutiple steps and try an animated histogram    
    #for line in reversed(list(open("C:\\LAMMPS\\tensile_DNW\\variable_pull_test.lammpstrj"))):
        #print(line)
        #if line.find("ITEM:") > -1:
            #break    
