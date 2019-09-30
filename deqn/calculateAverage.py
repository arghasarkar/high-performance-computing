import numpy as np

filename = "./2000.txt"

f = open(filename, "r")
data = f.read()

# Keys of the values in the text file
diffusion = "diff1"
explicit_scheme_diffusion = "esdiff"
explicit_scheme_reset = "esres"
mesh_get_total_temp = "mgtt"

diff = 0
esdiff_arr = []
esres_arr = []
mgtt_arr = []

data = data.split("\n")
for i in range(len(data)):
    d_split = data[i].split(":")

    if d_split[0] == diffusion:
        diff = d_split[1]

    if d_split[0] == explicit_scheme_diffusion:
        esdiff_arr.append(float(d_split[1]))

    if d_split[0] == explicit_scheme_reset:
        esres_arr.append(float(d_split[1]))

    if d_split[0] == mesh_get_total_temp:
        mgtt_arr.append(float(d_split[1]))


# Time to print out all the means
print("diff: {}".format(str(diff)))
print("esdiff: {}".format(str(np.mean(esdiff_arr))))
print("esres: {}".format(str(np.mean(esres_arr))))
print("mgtt: {}".format(str(np.mean(mgtt_arr))))
