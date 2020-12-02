
import numpy as np
import re
import pymol
from pymol import cmd
pymol.finish_launching()

########## set the sphere scales to corresponding value ##########
value_to_show = "Residence Time"

###### reading data from csv file ##########
binding_site_id = 10

fn_data = "./Interactions_DOPS.csv"
with open(fn_data, "r") as f:
    data_lines = f.readlines()

column_names = data_lines[0].strip().split(",")
for column_idx, column_name in enumerate(column_names):
    if column_name == "Residue":
        column_id_residue_set = column_idx
    elif column_name == "Residue idx":
        column_id_residue_index = column_idx
    elif column_name == "Binding site":
        column_id_BS = column_idx
    elif column_name == value_to_show:
        column_id_value_to_show = column_idx

residue_set = []
residue_rank_set = []
binding_site_identifiers = []
values_to_show = []
for line in data_lines[1:]:
    data_list = line.strip().split(",")
    residue_set.append(data_list[column_id_residue_set])
    residue_rank_set.append(data_list[column_id_residue_index])
    binding_site_identifiers.append(data_list[column_id_BS])
    values_to_show.append(data_list[column_id_value_to_show])

############## read information from pdb coordinates ##########
pdb_file = "./A2a_active_renumbered.pdb"
with open(pdb_file, "r") as f:
    pdb_lines = f.readlines()
residue_identifiers = []
for line in pdb_lines:
    if line.strip()[:4] == "ATOM":
        identifier = (line.strip()[22:26].strip(), line.strip()[17:20].strip(), line.strip()[21].strip())
##                           residue index,              resname,                     chain id
        if len(residue_identifiers) == 0:
            residue_identifiers.append(identifier)
        elif identifier != residue_identifiers[-1]:
            residue_identifiers.append(identifier)

######### calculate scale ###############
values_to_show = np.array(values_to_show, dtype=float)
MIN = np.percentile(np.unique(values_to_show), 5)
MAX = np.percentile(np.unique(values_to_show), 100)
print("MAX: {}".format(MAX))
print("MIN: {}".format(MIN))
#MIN = 0.05
#MAX = 17
X = (values_to_show - np.percentile(np.unique(values_to_show), 50))/(MAX - MIN)
SCALES = 1/(0.5 + np.exp(-X * 5))

######## some pymol settings #########
cmd.set("retain_order", 1)
cmd.set("cartoon_oval_length", 1.0)
cmd.set("cartoon_oval_width", 0.3)
cmd.set("cartoon_color", "white")
cmd.set("stick_radius", 0.35)

##################################
cmd.load("./A2a_active_renumbered.pdb", "Prot_DOPS")
prefix = "Prot_DOPS"
cmd.hide("everything")
cmd.show("cartoon", prefix)
cmd.center(prefix)
cmd.orient(prefix)
colors = np.array([np.random.choice(np.arange(256, dtype=float), size=3) for dummy in range(binding_site_id)])
colors /= 255.0

residue_set = np.array(residue_set, dtype=str)
residue_rank_set = np.array(residue_rank_set, dtype=int)
binding_site_identifiers = np.array(binding_site_identifiers, dtype=int)
residue_identifiers = list(residue_identifiers)
for bs_id in np.arange(binding_site_id):
    cmd.set_color("tmp_{}".format(bs_id), list(colors[bs_id]))
    for entry_id in np.where(binding_site_identifiers == bs_id)[0]:
        selected_residue = residue_set[entry_id]
        selected_residue_rank = residue_rank_set[entry_id]
        identifier_from_pdb = residue_identifiers[selected_residue_rank]
        if re.findall("[a-zA-Z]+$", selected_residue)[0] != identifier_from_pdb[1]:
            raise IndexError("The {}th residue in the provided pdb file ({}{}) is different from that in the simulations ({})!".format(entry_id+1,
                                                                                                                                     identifier_from_pdb[0],
                                                                                                                                     identifier_from_pdb[1],
                                                                                                                                     selected_residue))
        if identifier_from_pdb[2] != " ":
            cmd.select("BSid{}_{}".format(bs_id, selected_residue), "chain {} and resid {} and (not name C+O+N)".format(identifier_from_pdb[2],
                                                                                                                      identifier_from_pdb[0]))
        else:
            cmd.select("BSid{}_{}".format(bs_id, selected_residue), "resid {} and (not name C+O+N)".format(identifier_from_pdb[0]))
        cmd.show("spheres", "BSid{}_{}".format(bs_id, selected_residue))
        cmd.set("sphere_scale", SCALES[entry_id], selection="BSid{}_{}".format(bs_id, selected_residue))
        cmd.color("tmp_{}".format(bs_id), "BSid{}_{}".format(bs_id, selected_residue))
    cmd.group("BSid{}".format(bs_id), "BSid{}_*".format(bs_id))


